#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_Table.h"

extern "C" void CT_MultiLevel(CCTK_ARGUMENTS);
extern "C" void CT_V(CCTK_ARGUMENTS, const CCTK_INT toplevel, const CCTK_INT init_psi);
extern "C" void CT_FillFlagArrays(CCTK_INT *levels, CCTK_INT *downward, const CCTK_INT toplevel);
extern "C" void CT_InitializePsi(CCTK_ARGUMENTS);
extern "C" void CT_InitializeError(CCTK_ARGUMENTS);
extern "C" void CT_InitializeCoefficients(const int sindex, const int nequation, struct coeffptr *cptr, const struct stencil stnc);
extern "C" void CT_InitializeConstants(CCTK_ARGUMENTS);
extern "C" void CT_SolvePsiEquation(CCTK_ARGUMENTS, CCTK_REAL *norm, CCTK_INT *step);
extern "C" void CT_SolveErrorEquation(CCTK_ARGUMENTS, CCTK_REAL *norm, CCTK_INT *step);
extern "C" void CT_CalcPsiResidual(CCTK_ARGUMENTS, const CCTK_INT step, const CCTK_INT output, CCTK_REAL *norm, CCTK_INT outdated_psi);
extern "C" void CT_CalcErrResidual(CCTK_ARGUMENTS, const CCTK_INT step, const CCTK_INT output, CCTK_REAL *norm);
extern "C" void CT_CompareToExact(CCTK_ARGUMENTS, const CCTK_INT outdated_psi);
extern "C" void CT_AddErrorToPsi(CCTK_ARGUMENTS);
extern "C" void CT_RestoreError(CCTK_ARGUMENTS);
extern "C" void CT_PrintError(CCTK_ARGUMENTS, const CCTK_INT nequation);
extern "C" void CT_PrintPsi(CCTK_ARGUMENTS, const CCTK_INT nequation);
extern "C" void CT_CopyResidual(CCTK_ARGUMENTS, const CCTK_INT nequation);
extern "C" void CT_RelaxPsi(CCTK_ARGUMENTS, 
                            const CCTK_INT level,
                            const CCTK_INT nsteps,
                            const CCTK_INT init_psi,
                            const CCTK_INT prol_psi,
                            const CCTK_INT rest_psi,
                            const CCTK_INT add_err,
                            const CCTK_INT rset_psi,
                            const CCTK_INT init_coeff,
                            const CCTK_INT enf_int);
extern "C" void CT_RelaxError(CCTK_ARGUMENTS, 
                              const CCTK_INT level,
                              const CCTK_INT nsteps,
                              const CCTK_INT init_err,
                              const CCTK_INT prol_psi,
                              const CCTK_INT rest_psi,
                              const CCTK_INT rest_res,
                              const CCTK_INT add_err_psi,
                              const CCTK_INT rset_psi,
                              const CCTK_INT init_coeff);
extern "C" void CT_UpdateBoundaries(CCTK_ARGUMENTS, const char *varname);
extern "C" void CT_Reset(CCTK_ARGUMENTS, const char *varName, const CCTK_REAL *resetvalue, const CCTK_INT relative);
extern "C" void CT_InitializeResidual(CCTK_ARGUMENTS);
extern "C" void CT_EnforceInt(CCTK_ARGUMENTS, CCTK_REAL *value);

int CT_ProcessOwnsData();

extern "C" void CT_Restrict(CCTK_ARGUMENTS, const char *varname);
extern "C" void CT_Prolongate(CCTK_ARGUMENTS, const char *varname);
extern "C" void CT_ProlongateBndrs(CCTK_ARGUMENTS, const char *varname);
extern "C" void CT_FD(const CCTK_REAL *gfunc, const struct stencil stnc, const int nequation,
                      CCTK_REAL *derx, CCTK_REAL *dery, CCTK_REAL *derz,
                      CCTK_REAL *derxx, CCTK_REAL *deryy, CCTK_REAL *derzz,
                      CCTK_REAL *derxy, CCTK_REAL *derxz, CCTK_REAL *deryz);
extern "C" void CT_Norm(CCTK_ARGUMENTS, const char *varName, CCTK_REAL *norm, const CCTK_INT nequation);
extern "C" void CT_Integral(CCTK_ARGUMENTS, const char *varName, CCTK_REAL *norm, CCTK_INT refinement, CCTK_INT nequation);
CCTK_REAL CT_GetValue(CCTK_ARGUMENTS, CCTK_REAL xcoord, CCTK_REAL ycoord, CCTK_REAL zcoord, const char *varName);
extern "C" void CT_WriteTimeSeries(const CCTK_INT iteration, const CCTK_REAL value, const char *filename);
extern "C" void CT_WritePairs(const CCTK_REAL value1, const CCTK_REAL value2, const char *filename);
extern "C" void CT_OutputWalk(CCTK_ARGUMENTS);
extern "C" void CT_Copy(CCTK_ARGUMENTS, const char *varName1, const char *varName2, const CCTK_INT nequation);
extern "C" void CT_ClearFiles(CCTK_ARGUMENTS);
extern "C" void CT_InitializeADM(CCTK_ARGUMENTS);
extern "C" void CT_SetAuxiliaries(const int sindex, const struct stencil stnc, const struct coeffptr cptr);
extern "C" void CT_PopulatePointerStruct(CCTK_ARGUMENTS, struct coeffptr *cptr);
extern "C" void CT_Boundaries(CCTK_ARGUMENTS, const char *varname);

struct coeffptr
{
  CCTK_REAL *cxx;
  CCTK_REAL *cxy;
  CCTK_REAL *cxz;
  CCTK_REAL *cyy;
  CCTK_REAL *cyz;
  CCTK_REAL *czz;
  CCTK_REAL *cx;
  CCTK_REAL *cy;
  CCTK_REAL *cz;
  CCTK_REAL *c0;
  CCTK_REAL *c1;
  CCTK_REAL *c2;
  CCTK_REAL *c3;
  CCTK_REAL *c4;
  CCTK_REAL *a0;
  CCTK_REAL *a1;
  CCTK_REAL *a2;
  CCTK_REAL *a3;
  CCTK_REAL *a4;
  CCTK_REAL *ct_cxx;
  CCTK_REAL *ct_cxy;
  CCTK_REAL *ct_cxz;
  CCTK_REAL *ct_cyy;
  CCTK_REAL *ct_cyz;
  CCTK_REAL *ct_czz;
  CCTK_REAL *ct_cx;
  CCTK_REAL *ct_cy;
  CCTK_REAL *ct_cz;
  CCTK_REAL *ct_c0;
  CCTK_REAL *ct_c1;
  CCTK_REAL *ct_c2;
  CCTK_REAL *ct_c3;
  CCTK_REAL *ct_c4;
  CCTK_REAL *ct_a0;
  CCTK_REAL *ct_a1;
  CCTK_REAL *ct_a2;
  CCTK_REAL *ct_a3;
  CCTK_REAL *ct_a4;
  CCTK_REAL *ct_psi;
  CCTK_REAL *ct_auxiliary;
  CCTK_REAL **extras;
};

struct stencil
{
  int c;
  int mx;
  int px;
  int my;
  int py;
  int mz;
  int pz;
  int m2x;
  int p2x;
  int m2y;
  int p2y;
  int m2z;
  int p2z;
  int mxmy;
  int mxpy;
  int pxmy;
  int pxpy;
  int mxmz;
  int mxpz;
  int pxmz;
  int pxpz;
  int mymz;
  int mypz;
  int pymz;
  int pypz;
  CCTK_REAL dx;
  CCTK_REAL dy;
  CCTK_REAL dz;
  int order;
  int npoints;
};
