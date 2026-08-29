#ifndef RAMAX_KSW2_DISPATCH_H
#define RAMAX_KSW2_DISPATCH_H

#ifdef __cplusplus
extern "C" {
#endif

const char* ramax_ksw_selected_kernel(void);
int ramax_ksw_sse41_supported(void);
void ramax_ksw_set_kernel_for_testing(int kernel);

#ifdef __cplusplus
}
#endif

#endif  // RAMAX_KSW2_DISPATCH_H
