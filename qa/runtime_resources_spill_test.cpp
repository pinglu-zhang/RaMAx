#ifdef NDEBUG
#undef NDEBUG
#endif

#include "match_spill.h"
#include "mapped_file_region.h"
#include "ramesh.h"
#include "runtime_resources.h"

#include <cassert>
#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

bool sameMatches(const MatchVec3D& left, const MatchVec3D& right) {
  if (left.size() != right.size()) return false;
  for (size_t outer = 0; outer < left.size(); ++outer) {
    if (left[outer].size() != right[outer].size()) return false;
    for (size_t inner = 0; inner < left[outer].size(); ++inner) {
      if (left[outer][inner].size() != right[outer][inner].size()) {
        return false;
      }
      for (size_t index = 0; index < left[outer][inner].size(); ++index) {
        const Match& a = left[outer][inner][index];
        const Match& b = right[outer][inner][index];
        if (a.ref_chr_index != b.ref_chr_index ||
            a.ref_start != b.ref_start ||
            a.qry_chr_index != b.qry_chr_index ||
            a.qry_start != b.qry_start ||
            a.len_and_strand != b.len_and_strand) {
          return false;
        }
      }
    }
  }
  return true;
}

template <typename Function>
bool throws(Function&& function) {
  try {
    function();
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

}  // namespace

int main() {
  using namespace RaMAxResources;
  static_assert(sizeof(RaMesh::StageSharedLink) == sizeof(RaMesh::SegPtr));
  static_assert(sizeof(RaMesh::StageSharedMutex) == 1);
  // This is the compatibility backend, not the final stable-ID arena.  Keep a
  // hard ceiling so a future field/lock regression cannot silently restore the
  // former per-occurrence footprint.
  static_assert(sizeof(RaMesh::Segment) <= 104);
  static_assert(sizeof(RaMesh::Block) <= 96);
  static_assert(sizeof(RaMesh::AnchorPathKey) == 2 * sizeof(void*));

  {
    const RaMesh::InternedGraphName first(std::string("species-17"));
    const RaMesh::InternedGraphName second(std::string_view("species-17"));
    assert(&first.str() == &second.str());

    using LegacyMap = std::unordered_multimap<
        RaMesh::SpeciesChrPair, RaMesh::SegPtr, RaMesh::SpeciesChrPairHash>;
    LegacyMap legacy;
    RaMesh::ChrHeadMap compact;
    legacy.reserve(256);
    compact.reserve(256);

    std::vector<RaMesh::SpeciesChrPair> inserted;
    inserted.reserve(512);
    for (size_t index = 0; index < 512; ++index) {
      inserted.emplace_back(
          "species-" + std::to_string((index * 37) % 41),
          "chromosome-" + std::to_string((index * 53) % 23));
      legacy.emplace(inserted.back(), nullptr);
      compact.emplace(RaMesh::AnchorPathKey(inserted.back()), nullptr);
    }

    assert(legacy.bucket_count() == compact.bucket_count());
    assert(legacy.size() == compact.size());
    const RaMesh::SpeciesChrPairHash legacy_hash;
    const RaMesh::AnchorPathKeyHash compact_hash;
    for (const auto& key : inserted) {
      const RaMesh::AnchorPathKey interned(key);
      assert(legacy_hash(key) == compact_hash(interned));
      assert(legacy.bucket(key) == compact.bucket(interned));
    }

    auto legacy_it = legacy.begin();
    auto compact_it = compact.begin();
    for (; legacy_it != legacy.end(); ++legacy_it, ++compact_it) {
      assert(compact_it != compact.end());
      assert(legacy_it->first.first == compact_it->first.first.str());
      assert(legacy_it->first.second == compact_it->first.second.str());
    }
    assert(compact_it == compact.end());
  }

  assert(parseMemorySize("1KiB") == 1024ULL);
  assert(parseMemorySize("2MiB") == 2ULL * 1024ULL * 1024ULL);
  assert(parseMemorySize("3GiB") == 3ULL * 1024ULL * 1024ULL * 1024ULL);
  assert(parseMemorySize("1TiB") == 1024ULL * 1024ULL * 1024ULL * 1024ULL);
  assert(parseMemorySize("auto") == 0);
  assert(throws([] { (void)parseMemorySize("0GiB"); }));
  assert(throws([] { (void)parseMemorySize("1GB"); }));
  assert(throws([] { (void)parseMemorySize("1.5GiB"); }));

  const auto stamp = std::chrono::steady_clock::now()
      .time_since_epoch().count();
  const auto root = std::filesystem::current_path() /
      ("ramax-spill-test-" + std::to_string(stamp));
  std::filesystem::create_directories(root);
  RuntimeResourceManager::instance().configure(
      RuntimeResourceConfig{"auto", root / "spill", root, 4});

  MappedFileRegion mapped;
  mapped.allocate(1024 * 1024, "resource-test");
#if defined(__linux__)
  assert(mapped.fileBacked());
#endif
  auto* mapped_bytes = static_cast<unsigned char*>(mapped.data());
  mapped_bytes[0] = 17;
  mapped_bytes[mapped.size() - 1] = 29;
  mapped.adviseSequential();
  mapped.adviseRandom();
  MappedFileRegion moved = std::move(mapped);
  assert(mapped.data() == nullptr);
  assert(moved.size() == 1024 * 1024);
  assert(static_cast<unsigned char*>(moved.data())[0] == 17);
  assert(static_cast<unsigned char*>(moved.data())[moved.size() - 1] == 29);
  moved.clear();

  MatchVec3D source(4);
  source[1].resize(2);
  source[1][0].push_back(Match(7, 11, 13, 17, 19, FORWARD));
  source[1][0].push_back(Match(23, 29, 31, 37, 41, REVERSE));
  source[2].resize(1);  // Preserve an explicitly empty MatchVec.
  source[3].resize(3);
  source[3][2].push_back(Match(
      UINT32_MAX, UINT32_MAX - 1, UINT32_MAX - 2,
      UINT32_MAX - 3, 0x7fffffffU, REVERSE));

  RaMAxSpill::MatchSpillStore store(root / "spill");
  const RaMAxSpill::MatchSpillIdentity identity{
      9, 42, RaMAxSpill::MatchSpillKind::PRIMARY, "reference", "query"};
  const auto stats = store.inspect(source);
  assert(stats.outer_vectors == 4);
  assert(stats.inner_vectors == 6);
  assert(stats.matches == 3);
  const auto path = store.write(identity, source);
  assert(std::filesystem::is_regular_file(path));
  auto restored = store.read(path, identity);
  assert(restored);
  assert(sameMatches(source, *restored));
  assert(throws([&] {
    auto wrong = identity;
    ++wrong.query_ordinal;
    (void)store.read(path, wrong);
  }));
  store.consume(path);
  assert(!std::filesystem::exists(path));
  std::filesystem::remove_all(root);
  return 0;
}
