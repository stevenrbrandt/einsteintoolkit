#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "cctk_Schedule.h"

#include "carpet.hh"
#include "loopcontrol.h"
#include "Carpet/Carpet/src/modes.hh"

#include "CT_MultiLevel.hh"

using namespace Carpet;

extern "C" void CT_MultiLevel(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CT_ClearFiles(CCTK_PASS_CTOC);

  CT_InitializeResidual(CCTK_PASS_CTOC);
  CT_InitializeConstants(CCTK_PASS_CTOC);

  CCTK_INT rset_psi = CCTK_Equals(reset_psi, "no")? 0 : CCTK_Equals(reset_psi, "to value")? 1 : -1;
  CCTK_INT topMGl = topMGlevel > Carpet::reflevels - 1? Carpet::reflevels - 1: topMGlevel;

  if (CCTK_Equals(cycle_type, "V cycle"))
  {
    BEGIN_REFLEVEL_LOOP(cctkGH) {
      CT_InitializePsi(CCTK_PASS_CTOC);
      CT_InitializeError(CCTK_PASS_CTOC);
    } END_REFLEVEL_LOOP;
    
    CT_V(CCTK_PASS_CTOC, topMGl, 1);
  }
  else if (CCTK_Equals(cycle_type, "FMG cycle"))
  {
    CT_RelaxPsi(CCTK_PASS_CTOC, 0, nrelsteps_bottom, 1, 0, 0, 0, rset_psi, 1, enforce_int);

    for (int toplevel=1; toplevel<=topMGl; toplevel++)
      CT_V(CCTK_PASS_CTOC, toplevel, 0);

    for (int fmgit=1; fmgit < fmg_niter; fmgit++)
    {
      CCTK_Info(CCTK_THORNSTRING, "\t\t");
      CCTK_VInfo(CCTK_THORNSTRING, "\t\t=== FMG CYCLE #%d ===", fmgit+1);
      CCTK_Info(CCTK_THORNSTRING, "\t\t");

      for (int ilev=topMGl-1; ilev>=0; ilev--)
      {
        ENTER_LEVEL_MODE(cctkGH,ilev) {
          CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::psi");
          CT_Restrict(CCTK_PASS_CTOC, "CT_MultiLevel::coeffs");
          CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testK");
          CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdxK");
          CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdyK");
          CT_Restrict(CCTK_PASS_CTOC, "CT_Analytic::CT_testdzK");
        } LEAVE_LEVEL_MODE;
      }

      CT_RelaxPsi(CCTK_PASS_CTOC, 0, nrelsteps_bottom, 0, 0, 0, 0, rset_psi, 0, enforce_int);

      for (int toplevel=1; toplevel<=topMGl; toplevel++)
      {
        CT_V(CCTK_PASS_CTOC, toplevel, 0);
      }
    }
  }

  for (int rl=topMGl+1; rl<Carpet::reflevels; rl++)
    CT_RelaxPsi(CCTK_PASS_CTOC, rl, nrelsteps_top, 0, 1, 0, 0, 0, 1, 0);

  if (CCTK_Equals(fill_ADM, "yes")) CT_InitializeADM(CCTK_PASS_CTOC);

  return;
}
