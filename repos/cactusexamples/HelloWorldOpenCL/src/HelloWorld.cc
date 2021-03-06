#include <cctk.h>
#include <cctk_Arguments.h>

#include <OpenCLRunTime.h>

using namespace std;



// Autogenerated strings containing the OpenCL source code for our
// kernels
// Note: These names of these variables contain the thorn name
extern char const *const OpenCL_source_HelloWorldOpenCL_HelloWorld;



extern "C"
void HelloWorldOpenCL_evol(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  // The variable groups accessed by this kernel (defining the
  // kernel's CCTK_ARGUMENTS list)
  const char* const groups[] = {NULL};
  
  const int imin[] = {0, 0, 0};
  const int imax[] = {cctk_lsh[0], cctk_lsh[1], cctk_lsh[2]};
  
  static struct OpenCLKernel *kernel = NULL;
  const char* const sources[] =
    {"", OpenCL_source_HelloWorldOpenCL_HelloWorld, NULL};
  OpenCLRunTime_CallKernel(cctkGH, CCTK_THORNSTRING, "HelloWorld",
                           sources, groups, NULL, NULL, NULL, -1,
                           imin, imax, &kernel);
}
