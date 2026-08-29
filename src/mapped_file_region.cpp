#include "mapped_file_region.h"

#include "runtime_resources.h"

#include <cstdlib>
#include <filesystem>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

#if defined(__linux__)
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace RaMAxResources {
namespace {

std::string safeLabel(std::string_view label) {
  std::string result;
  result.reserve(label.size());
  for (const char character : label) {
    if ((character >= 'a' && character <= 'z') ||
        (character >= 'A' && character <= 'Z') ||
        (character >= '0' && character <= '9') || character == '-') {
      result.push_back(character);
    } else {
      result.push_back('-');
    }
  }
  if (result.empty()) result = "arena";
  return result;
}

#if defined(__linux__)
[[noreturn]] void throwSystemError(const std::string& operation) {
  throw std::runtime_error(
      operation + " failed: " + std::string(std::strerror(errno)));
}
#endif

}  // namespace

MappedFileRegion::MappedFileRegion(MappedFileRegion&& other) noexcept
    : data_(std::exchange(other.data_, nullptr)),
      size_(std::exchange(other.size_, 0)),
      file_backed_(std::exchange(other.file_backed_, false)) {}

MappedFileRegion& MappedFileRegion::operator=(
    MappedFileRegion&& other) noexcept {
  if (this != &other) {
    clear();
    data_ = std::exchange(other.data_, nullptr);
    size_ = std::exchange(other.size_, 0);
    file_backed_ = std::exchange(other.file_backed_, false);
  }
  return *this;
}

MappedFileRegion::~MappedFileRegion() { clear(); }

void MappedFileRegion::allocate(size_t bytes, std::string_view label) {
  clear();
  if (bytes == 0) return;

  auto& resources = RuntimeResourceManager::instance();
#if defined(__linux__)
  if (resources.configured()) {
    resources.requireTempSpace(bytes, label);
    const std::filesystem::path directory = resources.tempDirectory();
    std::string path =
        (directory / ("ramax-" + safeLabel(label) + "-XXXXXX")).string();
    const int descriptor = ::mkstemp(path.data());
    if (descriptor < 0) throwSystemError("mkstemp for " + std::string(label));

    struct FileGuard {
      int descriptor;
      std::string path;
      ~FileGuard() {
        if (descriptor >= 0) ::close(descriptor);
        if (!path.empty()) ::unlink(path.c_str());
      }
    } guard{descriptor, path};

    if (bytes > static_cast<size_t>(std::numeric_limits<off_t>::max())) {
      throw std::bad_array_new_length();
    }
    if (::ftruncate(descriptor, static_cast<off_t>(bytes)) != 0) {
      throwSystemError("ftruncate for " + std::string(label));
    }
    void* mapping = ::mmap(
        nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, descriptor, 0);
    if (mapping == MAP_FAILED) {
      throwSystemError("mmap for " + std::string(label));
    }
    data_ = mapping;
    size_ = bytes;
    file_backed_ = true;
    return;
  }
#else
  (void)label;
#endif

  data_ = std::malloc(bytes);
  if (data_ == nullptr) throw std::bad_alloc();
  size_ = bytes;
  file_backed_ = false;
}

void MappedFileRegion::clear() noexcept {
  if (data_ == nullptr) return;
#if defined(__linux__)
  if (file_backed_) {
    (void)::munmap(data_, size_);
  } else {
    std::free(data_);
  }
#else
  std::free(data_);
#endif
  data_ = nullptr;
  size_ = 0;
  file_backed_ = false;
}

void MappedFileRegion::adviseSequential() noexcept {
#if defined(__linux__)
  if (file_backed_ && data_ != nullptr) {
    (void)::madvise(data_, size_, MADV_SEQUENTIAL);
  }
#endif
}

void MappedFileRegion::adviseRandom() noexcept {
#if defined(__linux__)
  if (file_backed_ && data_ != nullptr) {
    (void)::madvise(data_, size_, MADV_RANDOM);
  }
#endif
}

void MappedFileRegion::adviseDontNeed() noexcept {
#if defined(__linux__)
  if (file_backed_ && data_ != nullptr) {
    (void)::madvise(data_, size_, MADV_DONTNEED);
  }
#endif
}

}  // namespace RaMAxResources
