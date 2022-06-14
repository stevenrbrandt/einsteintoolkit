#ifndef _VECTORS_TESTS_H_
#define _VECTORS_TESTS_H_
#include <cctk.h>

namespace Vectors {
bool Test_4_AVX(CCTK_ARGUMENTS);
bool Test_4_AVX_AVX2(CCTK_ARGUMENTS);
bool Test_4_AVX_FMA4(CCTK_ARGUMENTS);
bool Test_4_Altivec(CCTK_ARGUMENTS);
bool Test_4_SSE(CCTK_ARGUMENTS);
bool Test_4_SSE_SSE41(CCTK_ARGUMENTS);
bool Test_4_SSE_SSE4A(CCTK_ARGUMENTS);
bool Test_4_SSE_FMA4(CCTK_ARGUMENTS);
bool Test_4_SSE(CCTK_ARGUMENTS);
bool Test_4_default(CCTK_ARGUMENTS);
bool Test_8_AVX(CCTK_ARGUMENTS);
bool Test_8_AVX_AVX2(CCTK_ARGUMENTS);
bool Test_8_AVX_FMA4(CCTK_ARGUMENTS);
bool Test_8_AVX512(CCTK_ARGUMENTS);
bool Test_8_AVX512_AVX512ER(CCTK_ARGUMENTS);
bool Test_8_MIC(CCTK_ARGUMENTS);
bool Test_8_SSE2(CCTK_ARGUMENTS);
bool Test_8_SSE2_SSE41(CCTK_ARGUMENTS);
bool Test_8_SSE2_SSE4A(CCTK_ARGUMENTS);
bool Test_8_SSE2_FMA4(CCTK_ARGUMENTS);
bool Test_8_VSX(CCTK_ARGUMENTS);
bool Test_8_default(CCTK_ARGUMENTS);
}

#endif // _VECTORS_TESTS_H_
