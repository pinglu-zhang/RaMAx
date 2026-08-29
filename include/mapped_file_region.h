#ifndef RAMAX_MAPPED_FILE_REGION_H
#define RAMAX_MAPPED_FILE_REGION_H

#include <cstddef>
#include <string_view>

namespace RaMAxResources {

// An uninitialized byte region.  Once the runtime resource manager is
// configured on Linux the storage is backed by an unlinked temporary file,
// allowing clean pages to be reclaimed by the kernel without changing the
// array interface used by the suffix-array builders.  Other environments use
// malloc with identical lifetime and alignment semantics.
class MappedFileRegion {
 public:
  MappedFileRegion() = default;
  MappedFileRegion(const MappedFileRegion&) = delete;
  MappedFileRegion& operator=(const MappedFileRegion&) = delete;
  MappedFileRegion(MappedFileRegion&& other) noexcept;
  MappedFileRegion& operator=(MappedFileRegion&& other) noexcept;
  ~MappedFileRegion();

  void allocate(size_t bytes, std::string_view label);
  void clear() noexcept;

  [[nodiscard]] void* data() noexcept { return data_; }
  [[nodiscard]] const void* data() const noexcept { return data_; }
  [[nodiscard]] size_t size() const noexcept { return size_; }
  [[nodiscard]] bool fileBacked() const noexcept { return file_backed_; }

  void adviseSequential() noexcept;
  void adviseRandom() noexcept;
  void adviseDontNeed() noexcept;

 private:
  void* data_{nullptr};
  size_t size_{0};
  bool file_backed_{false};
};

}  // namespace RaMAxResources

#endif  // RAMAX_MAPPED_FILE_REGION_H
