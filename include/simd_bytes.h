#ifndef RAMAX_SIMD_BYTES_H
#define RAMAX_SIMD_BYTES_H

#include <cstddef>
#include <cstdint>

namespace RaMAxSimd {

enum class ByteKernel : uint8_t {
  SCALAR = 0,
  SSE41 = 1,
  AVX2 = 2
};

struct ByteComparison {
  int order{0};
  size_t lcp{0};

  bool operator==(const ByteComparison&) const = default;
};

ByteComparison compareBytes(const uint8_t* text, size_t text_size,
                            const uint8_t* pattern, size_t pattern_size,
                            size_t known_lcp = 0) noexcept;

ByteComparison compareBytesWithKernel(
    ByteKernel kernel, const uint8_t* text, size_t text_size,
    const uint8_t* pattern, size_t pattern_size,
    size_t known_lcp = 0) noexcept;

bool byteKernelSupported(ByteKernel kernel) noexcept;
ByteKernel selectedByteKernel() noexcept;
const char* byteKernelName(ByteKernel kernel) noexcept;

}  // namespace RaMAxSimd

#endif  // RAMAX_SIMD_BYTES_H
