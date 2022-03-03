#undef __AVX__
//#undef __knl__
//#undef __MIC__
#undef __AVX512F__
#undef __AVX512ER__
#undef __SSE2__
#undef __SSE__
#undef __bgq__
#undef __VECTOR4DOUBLE__
#undef __ALTIVEC__
#undef _ARCH_PWR7
#undef _ARCH_450D

#define VECTOR_REAL_PRECISION 8

// includes cctk.h of correct type
#include "test.hcc"

#include "test-vectors.h"

namespace Vectors {

bool Test_8_MIC(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

#if (defined __MIC__ || defined __knl__) && !defined DISABLE_AVX512 // Intel MIC
  return Test(CCTK_PASS_CTOC, vec_architecture);
#else
  return true;
#endif
}

}
