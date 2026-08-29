#include "ksw2.h"
#include "ksw2_dispatch.h"

#include <atomic>
#include <cstdint>

extern "C" {
void ksw_extz2_sse2(void*, int, const uint8_t*, int, const uint8_t*, int8_t,
                    const int8_t*, int8_t, int8_t, int, int, int, int,
                    ksw_extz_t*);
void ksw_extz2_sse41(void*, int, const uint8_t*, int, const uint8_t*, int8_t,
                     const int8_t*, int8_t, int8_t, int, int, int, int,
                     ksw_extz_t*);
}

namespace {

using KswFunction = void (*)(
    void*, int, const uint8_t*, int, const uint8_t*, int8_t,
    const int8_t*, int8_t, int8_t, int, int, int, int, ksw_extz_t*);

bool supportsSse41() noexcept {
#if (defined(__GNUC__) || defined(__clang__)) && \
    (defined(__x86_64__) || defined(__i386__))
  __builtin_cpu_init();
  return __builtin_cpu_supports("sse4.1");
#else
  return false;
#endif
}

KswFunction automaticKernel() noexcept {
  static KswFunction function =
      supportsSse41() ? ksw_extz2_sse41 : ksw_extz2_sse2;
  return function;
}

std::atomic<int> test_override{-1};

KswFunction selectedKernel() noexcept {
  const int override = test_override.load(std::memory_order_relaxed);
  if (override == 0) return ksw_extz2_sse2;
  if (override == 1 && supportsSse41()) return ksw_extz2_sse41;
  return automaticKernel();
}

}  // namespace

extern "C" void ksw_extz2_sse(
    void* km, int qlen, const uint8_t* query, int tlen,
    const uint8_t* target, int8_t m, const int8_t* mat, int8_t q,
    int8_t e, int w, int zdrop, int end_bonus, int flag,
    ksw_extz_t* ez) {
  selectedKernel()(km, qlen, query, tlen, target, m, mat, q, e, w,
                   zdrop, end_bonus, flag, ez);
}

extern "C" const char* ramax_ksw_selected_kernel(void) {
  return automaticKernel() == ksw_extz2_sse41 ? "SSE4.1" : "SSE2";
}

extern "C" int ramax_ksw_sse41_supported(void) {
  return supportsSse41() ? 1 : 0;
}

extern "C" void ramax_ksw_set_kernel_for_testing(int kernel) {
  test_override.store(kernel, std::memory_order_relaxed);
}
