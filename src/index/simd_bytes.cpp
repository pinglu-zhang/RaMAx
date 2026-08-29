#include "simd_bytes.h"

#include <algorithm>

#if (defined(__GNUC__) || defined(__clang__)) && \
    (defined(__x86_64__) || defined(__i386__))
#define RAMAX_X86_TARGET_DISPATCH 1
#include <immintrin.h>
#else
#define RAMAX_X86_TARGET_DISPATCH 0
#endif

namespace RaMAxSimd {
namespace {

using CompareFunction = ByteComparison (*)(
    const uint8_t*, size_t, const uint8_t*, size_t, size_t) noexcept;

ByteComparison compareScalar(const uint8_t* text, size_t text_size,
                             const uint8_t* pattern, size_t pattern_size,
                             size_t known_lcp) noexcept {
  size_t index = std::min(known_lcp, pattern_size);
  while (index < pattern_size) {
    if (index >= text_size) return {-1, index};
    if (text[index] != pattern[index]) {
      return {text[index] < pattern[index] ? -1 : 1, index};
    }
    ++index;
  }
  return {0, index};
}

#if RAMAX_X86_TARGET_DISPATCH
__attribute__((target("sse4.1")))
ByteComparison compareSse41(const uint8_t* text, size_t text_size,
                            const uint8_t* pattern, size_t pattern_size,
                            size_t known_lcp) noexcept {
  size_t index = std::min(known_lcp, pattern_size);
  const size_t scalar_end =
      std::min(pattern_size, index + static_cast<size_t>(8));
  while (index < scalar_end) {
    if (index >= text_size) return {-1, index};
    if (text[index] != pattern[index]) {
      return {text[index] < pattern[index] ? -1 : 1, index};
    }
    ++index;
  }
  while (index + 16 <= pattern_size && index + 16 <= text_size) {
    const __m128i left = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(text + index));
    const __m128i right = _mm_loadu_si128(
        reinterpret_cast<const __m128i*>(pattern + index));
    const __m128i equal = _mm_cmpeq_epi8(left, right);
    if (!_mm_test_all_ones(equal)) {
      const unsigned mask =
          static_cast<unsigned>(_mm_movemask_epi8(equal));
      const size_t mismatch = static_cast<size_t>(
          __builtin_ctz((~mask) & 0xffffU));
      const size_t position = index + mismatch;
      return {text[position] < pattern[position] ? -1 : 1, position};
    }
    index += 16;
  }
  while (index < pattern_size) {
    if (index >= text_size) return {-1, index};
    if (text[index] != pattern[index]) {
      return {text[index] < pattern[index] ? -1 : 1, index};
    }
    ++index;
  }
  return {0, index};
}

__attribute__((target("avx2")))
ByteComparison compareAvx2(const uint8_t* text, size_t text_size,
                           const uint8_t* pattern, size_t pattern_size,
                           size_t known_lcp) noexcept {
  size_t index = std::min(known_lcp, pattern_size);
  const size_t scalar_end =
      std::min(pattern_size, index + static_cast<size_t>(8));
  while (index < scalar_end) {
    if (index >= text_size) return {-1, index};
    if (text[index] != pattern[index]) {
      return {text[index] < pattern[index] ? -1 : 1, index};
    }
    ++index;
  }
  while (index + 32 <= pattern_size && index + 32 <= text_size) {
    const __m256i left = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(text + index));
    const __m256i right = _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(pattern + index));
    const uint32_t mask = static_cast<uint32_t>(
        _mm256_movemask_epi8(_mm256_cmpeq_epi8(left, right)));
    if (mask != 0xffffffffU) {
      const size_t mismatch = static_cast<size_t>(__builtin_ctz(~mask));
      const size_t position = index + mismatch;
      return {text[position] < pattern[position] ? -1 : 1, position};
    }
    index += 32;
  }
  while (index < pattern_size) {
    if (index >= text_size) return {-1, index};
    if (text[index] != pattern[index]) {
      return {text[index] < pattern[index] ? -1 : 1, index};
    }
    ++index;
  }
  return {0, index};
}
#endif

ByteKernel detectKernel() noexcept {
#if RAMAX_X86_TARGET_DISPATCH
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx2")) return ByteKernel::AVX2;
  if (__builtin_cpu_supports("sse4.1")) return ByteKernel::SSE41;
#endif
  return ByteKernel::SCALAR;
}

CompareFunction functionFor(ByteKernel kernel) noexcept {
#if RAMAX_X86_TARGET_DISPATCH
  if (kernel == ByteKernel::AVX2 && byteKernelSupported(kernel)) {
    return compareAvx2;
  }
  if (kernel == ByteKernel::SSE41 && byteKernelSupported(kernel)) {
    return compareSse41;
  }
#endif
  return compareScalar;
}

}  // namespace

bool byteKernelSupported(ByteKernel kernel) noexcept {
  if (kernel == ByteKernel::SCALAR) return true;
#if RAMAX_X86_TARGET_DISPATCH
  __builtin_cpu_init();
  if (kernel == ByteKernel::SSE41) return __builtin_cpu_supports("sse4.1");
  if (kernel == ByteKernel::AVX2) return __builtin_cpu_supports("avx2");
#endif
  return false;
}

ByteKernel selectedByteKernel() noexcept {
  static const ByteKernel kernel = detectKernel();
  return kernel;
}

const char* byteKernelName(ByteKernel kernel) noexcept {
  switch (kernel) {
    case ByteKernel::SCALAR: return "scalar";
    case ByteKernel::SSE41: return "SSE4.1";
    case ByteKernel::AVX2: return "AVX2";
  }
  return "unknown";
}

ByteComparison compareBytesWithKernel(
    ByteKernel kernel, const uint8_t* text, size_t text_size,
    const uint8_t* pattern, size_t pattern_size,
    size_t known_lcp) noexcept {
  return functionFor(kernel)(
      text, text_size, pattern, pattern_size, known_lcp);
}

ByteComparison compareBytes(const uint8_t* text, size_t text_size,
                            const uint8_t* pattern, size_t pattern_size,
                            size_t known_lcp) noexcept {
  static const CompareFunction function = functionFor(selectedByteKernel());
  return function(text, text_size, pattern, pattern_size, known_lcp);
}

}  // namespace RaMAxSimd
