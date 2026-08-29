#ifdef NDEBUG
#undef NDEBUG
#endif

#include "align.h"
#include "ksw2_dispatch.h"
#include "simd_bytes.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main() {
  std::mt19937 generator(0x5eed41u);
  std::uniform_int_distribution<int> byte_distribution(0, 255);
  for (size_t fixture = 0; fixture < 100000; ++fixture) {
    const size_t text_size = generator() % 513;
    const size_t pattern_size = generator() % 513;
    std::vector<uint8_t> text(text_size);
    std::vector<uint8_t> pattern(pattern_size);
    for (uint8_t& byte : text) byte = byte_distribution(generator);
    for (uint8_t& byte : pattern) byte = byte_distribution(generator);
    if (fixture % 3 == 0) {
      const size_t copied = std::min(text_size, pattern_size);
      std::copy_n(text.begin(), copied, pattern.begin());
      if (copied != 0 && fixture % 2 == 0) {
        pattern[copied - 1] ^= 1U;
      }
    }
    size_t actual_lcp = 0;
    while (actual_lcp < text_size && actual_lcp < pattern_size &&
           text[actual_lcp] == pattern[actual_lcp]) {
      ++actual_lcp;
    }
    const size_t known = actual_lcp == 0
        ? 0
        : generator() % (actual_lcp + 1);
    const auto scalar = RaMAxSimd::compareBytesWithKernel(
        RaMAxSimd::ByteKernel::SCALAR,
        text.data(), text.size(), pattern.data(), pattern.size(), known);
    for (const auto kernel : {
             RaMAxSimd::ByteKernel::SSE41,
             RaMAxSimd::ByteKernel::AVX2}) {
      if (!RaMAxSimd::byteKernelSupported(kernel)) continue;
      assert(RaMAxSimd::compareBytesWithKernel(
                 kernel, text.data(), text.size(), pattern.data(),
                 pattern.size(), known) == scalar);
    }
    assert(RaMAxSimd::compareBytes(
               text.data(), text.size(), pattern.data(), pattern.size(),
               known) == scalar);
  }

  if (ramax_ksw_sse41_supported()) {
    constexpr std::string_view alphabet = "ACGTN";
    for (size_t fixture = 0; fixture < 10000; ++fixture) {
      std::string reference(1 + generator() % 128, 'A');
      std::string query(1 + generator() % 128, 'A');
      for (char& base : reference) base = alphabet[generator() % 5];
      for (char& base : query) base = alphabet[generator() % 5];
      const bool reverse = (generator() & 1U) != 0;
      ramax_ksw_set_kernel_for_testing(0);
      const AlignmentResult sse2 = ramaxExtendAlignKSW2RawForTesting(
          {reference, false}, {query, reverse}, 400);
      ramax_ksw_set_kernel_for_testing(1);
      const AlignmentResult sse41 = ramaxExtendAlignKSW2RawForTesting(
          {reference, false}, {query, reverse}, 400);
      assert(sse2.cigar == sse41.cigar);
      assert(sse2.summary.reference_length == sse41.summary.reference_length);
      assert(sse2.summary.query_length == sse41.summary.query_length);
      assert(sse2.summary.alignment_length == sse41.summary.alignment_length);
      assert(sse2.summary.match_length == sse41.summary.match_length);
    }
  }
  ramax_ksw_set_kernel_for_testing(-1);
  std::cout << "byte-kernel="
            << RaMAxSimd::byteKernelName(RaMAxSimd::selectedByteKernel())
            << " ksw-kernel=" << ramax_ksw_selected_kernel() << '\n';
  return 0;
}
