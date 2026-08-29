#include "match_spill.h"

#include "cache_manifest.h"
#include "runtime_resources.h"
#include "spdlog/spdlog.h"

#include <array>
#include <chrono>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#if defined(__linux__)
#include <unistd.h>
#endif

namespace RaMAxSpill {
namespace {

constexpr std::array<char, 8> kMagic{'R', 'M', 'X', 'M', 'A', 'T', 'C', 'H'};
constexpr uint32_t kFormatVersion = 1;
constexpr uint32_t kEndianMarker = 0x01020304U;
constexpr uint64_t kFixedHeaderBytes =
    kMagic.size() + 4 + 4 + 8 + 8 + 4 + 4 + 8 + 8 + 8 + 8 + 8 + 8;
constexpr uint64_t kMatchBytes = 5ULL * sizeof(uint32_t);

uint64_t checkedAdd(uint64_t left, uint64_t right) {
  if (right > std::numeric_limits<uint64_t>::max() - left) {
    throw std::overflow_error("Match spill size exceeds uint64_t");
  }
  return left + right;
}

uint64_t checkedMultiply(uint64_t left, uint64_t right) {
  if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left) {
    throw std::overflow_error("Match spill size exceeds uint64_t");
  }
  return left * right;
}

template <typename UInt>
void writeUnsigned(std::ostream& output, UInt value) {
  static_assert(std::is_unsigned_v<UInt>);
  std::array<unsigned char, sizeof(UInt)> bytes{};
  for (size_t index = 0; index < bytes.size(); ++index) {
    bytes[index] = static_cast<unsigned char>(value & 0xffU);
    value >>= 8U;
  }
  output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

template <typename UInt>
UInt readUnsigned(std::istream& input, std::string_view field) {
  static_assert(std::is_unsigned_v<UInt>);
  std::array<unsigned char, sizeof(UInt)> bytes{};
  input.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
  if (!input) {
    throw std::runtime_error("Truncated Match spill while reading " +
                             std::string(field));
  }
  UInt value = 0;
  for (size_t index = bytes.size(); index != 0; --index) {
    value = static_cast<UInt>((value << 8U) | bytes[index - 1]);
  }
  return value;
}

void writeString(std::ostream& output, std::string_view value) {
  output.write(value.data(), static_cast<std::streamsize>(value.size()));
}

std::string readString(std::istream& input, uint64_t length,
                       std::string_view field) {
  if (length > std::numeric_limits<size_t>::max()) {
    throw std::runtime_error("Match spill string is too large: " +
                             std::string(field));
  }
  std::string result(static_cast<size_t>(length), '\0');
  input.read(result.data(), static_cast<std::streamsize>(result.size()));
  if (!input) {
    throw std::runtime_error("Truncated Match spill while reading " +
                             std::string(field));
  }
  return result;
}

std::string kindName(MatchSpillKind kind) {
  switch (kind) {
    case MatchSpillKind::PRIMARY: return "primary";
    case MatchSpillKind::REPEAT_MASKED: return "repeat-masked";
    case MatchSpillKind::REPEAT_FULL: return "repeat-full";
    case MatchSpillKind::SECONDARY_COMBINED: return "secondary-combined";
  }
  throw std::invalid_argument("Unknown Match spill kind");
}

uint64_t processId() noexcept {
#if defined(__linux__)
  return static_cast<uint64_t>(::getpid());
#else
  return 0;
#endif
}

void validateIdentity(const MatchSpillIdentity& actual,
                      const MatchSpillIdentity& expected) {
  if (actual.round != expected.round ||
      actual.query_ordinal != expected.query_ordinal ||
      actual.kind != expected.kind ||
      actual.reference != expected.reference ||
      actual.query != expected.query) {
    throw std::runtime_error(
        "Match spill identity does not match the requested search result");
  }
}

}  // namespace

MatchSpillStore::MatchSpillStore(std::filesystem::path directory)
    : directory_(std::move(directory)) {
  std::error_code error;
  std::filesystem::create_directories(directory_, error);
  if (error || !std::filesystem::is_directory(directory_)) {
    throw std::runtime_error("Cannot create Match spill directory " +
                             directory_.string() + ": " + error.message());
  }
}

MatchSpillStats MatchSpillStore::inspect(const MatchVec3D& matches) const {
  MatchSpillStats stats;
  stats.outer_vectors = matches.size();
  stats.bytes = kFixedHeaderBytes;
  for (const auto& outer : matches) {
    stats.inner_vectors = checkedAdd(stats.inner_vectors, outer.size());
    stats.bytes = checkedAdd(stats.bytes, sizeof(uint64_t));
    for (const auto& inner : outer) {
      stats.matches = checkedAdd(stats.matches, inner.size());
      stats.bytes = checkedAdd(stats.bytes, sizeof(uint64_t));
    }
  }
  stats.bytes = checkedAdd(
      stats.bytes, checkedMultiply(stats.matches, kMatchBytes));
  return stats;
}

std::filesystem::path MatchSpillStore::pathFor(
    const MatchSpillIdentity& identity) const {
  return directory_ /
      ("match-p" + std::to_string(processId()) + "-r" +
       std::to_string(identity.round) + "-q" +
       std::to_string(identity.query_ordinal) + "-" +
       kindName(identity.kind) + ".bin");
}

std::filesystem::path MatchSpillStore::write(
    const MatchSpillIdentity& identity, const MatchVec3D& matches) const {
  const auto start = std::chrono::steady_clock::now();
  MatchSpillStats stats = inspect(matches);
  stats.bytes = checkedAdd(stats.bytes, identity.reference.size());
  stats.bytes = checkedAdd(stats.bytes, identity.query.size());
  RaMAxResources::RuntimeResourceManager::instance().requireTempSpace(
      stats.bytes, "match-spill-write");

  const auto final_path = pathFor(identity);
  auto partial_path = final_path;
  partial_path += ".partial";
  auto incomplete_path = final_path;
  incomplete_path += ".incomplete";
  RaMAxCache::removeIfPresent(partial_path);

  {
    std::ofstream marker(incomplete_path, std::ios::trunc);
    marker << "RAMAX_MATCH_SPILL_INCOMPLETE 1\n"
           << identity.round << ' ' << identity.query_ordinal << ' '
           << static_cast<uint32_t>(identity.kind) << '\n';
  }

  try {
    std::ofstream output(partial_path, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error("Cannot create Match spill " +
                               partial_path.string());
    }
    output.write(kMagic.data(), kMagic.size());
    writeUnsigned<uint32_t>(output, kFormatVersion);
    writeUnsigned<uint32_t>(output, kEndianMarker);
    writeUnsigned<uint64_t>(output, identity.round);
    writeUnsigned<uint64_t>(output, identity.query_ordinal);
    writeUnsigned<uint32_t>(output, static_cast<uint32_t>(identity.kind));
    writeUnsigned<uint32_t>(output, 0);
    writeUnsigned<uint64_t>(output, stats.outer_vectors);
    writeUnsigned<uint64_t>(output, stats.inner_vectors);
    writeUnsigned<uint64_t>(output, stats.matches);
    writeUnsigned<uint64_t>(output, stats.bytes);
    writeUnsigned<uint64_t>(output, identity.reference.size());
    writeUnsigned<uint64_t>(output, identity.query.size());
    writeString(output, identity.reference);
    writeString(output, identity.query);
    for (const auto& outer : matches) {
      writeUnsigned<uint64_t>(output, outer.size());
      for (const auto& inner : outer) {
        writeUnsigned<uint64_t>(output, inner.size());
        for (const Match& match : inner) {
          writeUnsigned<uint32_t>(output, match.ref_chr_index);
          writeUnsigned<uint32_t>(output, match.ref_start);
          writeUnsigned<uint32_t>(output, match.qry_chr_index);
          writeUnsigned<uint32_t>(output, match.qry_start);
          writeUnsigned<uint32_t>(output, match.len_and_strand);
        }
      }
    }
    output.flush();
    if (!output) {
      throw std::runtime_error("Failed to write Match spill " +
                               partial_path.string());
    }
    output.close();
    const uint64_t actual_size = std::filesystem::file_size(partial_path);
    if (actual_size != stats.bytes) {
      throw std::runtime_error("Match spill length mismatch before publish");
    }
    RaMAxCache::publishFile(partial_path, final_path);
    RaMAxCache::removeIfPresent(incomplete_path);
  } catch (...) {
    // Keep both the partial payload (when present) and diagnostic marker.
    throw;
  }

  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  spdlog::info(
      "[spill] event=write species={} reference={} round={} ordinal={} "
      "kind={} outer={} inner={} records={} bytes={} seconds={:.6f} path={}",
      identity.query, identity.reference, identity.round,
      identity.query_ordinal, kindName(identity.kind), stats.outer_vectors,
      stats.inner_vectors, stats.matches, stats.bytes, seconds,
      final_path.string());
  return final_path;
}

MatchVec3DPtr MatchSpillStore::read(
    const std::filesystem::path& path,
    const MatchSpillIdentity& expected) const {
  const auto start = std::chrono::steady_clock::now();
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("Cannot open Match spill " + path.string());
  }
  std::array<char, kMagic.size()> magic{};
  input.read(magic.data(), magic.size());
  if (!input || magic != kMagic) {
    throw std::runtime_error("Invalid Match spill magic: " + path.string());
  }
  const uint32_t version = readUnsigned<uint32_t>(input, "version");
  const uint32_t endian = readUnsigned<uint32_t>(input, "endian marker");
  if (version != kFormatVersion || endian != kEndianMarker) {
    throw std::runtime_error("Unsupported Match spill format: " +
                             path.string());
  }
  MatchSpillIdentity actual;
  actual.round = readUnsigned<uint64_t>(input, "round");
  actual.query_ordinal = readUnsigned<uint64_t>(input, "query ordinal");
  actual.kind = static_cast<MatchSpillKind>(
      readUnsigned<uint32_t>(input, "kind"));
  (void)readUnsigned<uint32_t>(input, "reserved");
  const uint64_t outer_count = readUnsigned<uint64_t>(input, "outer count");
  const uint64_t expected_inner = readUnsigned<uint64_t>(input, "inner count");
  const uint64_t expected_matches = readUnsigned<uint64_t>(input, "match count");
  const uint64_t expected_bytes = readUnsigned<uint64_t>(input, "file length");
  const uint64_t reference_length =
      readUnsigned<uint64_t>(input, "reference name length");
  const uint64_t query_length =
      readUnsigned<uint64_t>(input, "query name length");
  actual.reference = readString(input, reference_length, "reference name");
  actual.query = readString(input, query_length, "query name");
  validateIdentity(actual, expected);

  const uint64_t actual_bytes = std::filesystem::file_size(path);
  if (actual_bytes != expected_bytes) {
    throw std::runtime_error("Match spill file length is inconsistent: " +
                             path.string());
  }
  if (outer_count > std::numeric_limits<size_t>::max()) {
    throw std::runtime_error("Match spill outer vector count is too large");
  }
  auto result = std::make_shared<MatchVec3D>();
  result->resize(static_cast<size_t>(outer_count));
  uint64_t inner_count = 0;
  uint64_t match_count = 0;
  for (auto& outer : *result) {
    const uint64_t count = readUnsigned<uint64_t>(input, "inner vector count");
    inner_count = checkedAdd(inner_count, count);
    if (count > std::numeric_limits<size_t>::max()) {
      throw std::runtime_error("Match spill inner vector count is too large");
    }
    outer.resize(static_cast<size_t>(count));
    for (auto& inner : outer) {
      const uint64_t records = readUnsigned<uint64_t>(input, "record count");
      match_count = checkedAdd(match_count, records);
      if (records > std::numeric_limits<size_t>::max()) {
        throw std::runtime_error("Match spill record count is too large");
      }
      inner.resize(static_cast<size_t>(records));
      for (Match& match : inner) {
        match.ref_chr_index = readUnsigned<uint32_t>(input, "ref chromosome");
        match.ref_start = readUnsigned<uint32_t>(input, "ref start");
        match.qry_chr_index = readUnsigned<uint32_t>(input, "query chromosome");
        match.qry_start = readUnsigned<uint32_t>(input, "query start");
        match.len_and_strand = readUnsigned<uint32_t>(input, "length/strand");
      }
    }
  }
  if (inner_count != expected_inner || match_count != expected_matches ||
      input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("Match spill counts or trailing data are invalid: " +
                             path.string());
  }
  const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  spdlog::info(
      "[spill] event=read species={} reference={} round={} ordinal={} "
      "kind={} outer={} inner={} records={} bytes={} seconds={:.6f} path={}",
      actual.query, actual.reference, actual.round, actual.query_ordinal,
      kindName(actual.kind), outer_count, inner_count, match_count,
      actual_bytes, seconds, path.string());
  return result;
}

void MatchSpillStore::consume(const std::filesystem::path& path) const noexcept {
  RaMAxCache::removeIfPresent(path);
  auto incomplete_path = path;
  incomplete_path += ".incomplete";
  RaMAxCache::removeIfPresent(incomplete_path);
}

}  // namespace RaMAxSpill
