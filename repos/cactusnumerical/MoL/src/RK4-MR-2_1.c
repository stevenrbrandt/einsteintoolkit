 /*@@
   @file      RK4-MR-2_1.c
   @date      2012-03-25
   @author    Christian Reisswig
   @desc 
   A routine to perform 3rd order 2:1 multirate RK evolution. Mostly copied
   from genericRK.c
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>
#include "ExternalVariables.h"

//#define MOLDEBUG

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_RK4_MR_2_1_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK4_MR_2_1_Add(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    MoL_RK4_MR_2_1_Add
   @date       
   @author     
   @desc 
   Performs a single step of a RK43 2:1 type time
   integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK4_MR_2_1_Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS_CHECKED(MoL_RK4_MR_2_1_Add)
  DECLARE_CCTK_PARAMETERS;
    
  CCTK_INT arraydim;

  static CCTK_INT scratchspace_firstindex = -99;
  static CCTK_INT scratchspace_firstindex_slow = -99;
  CCTK_INT index, var, scratchstep;
  CCTK_INT totalsize;
  CCTK_REAL alpha[9], beta[10];
  CCTK_REAL alpha_slow[9], beta_slow[10];
  CCTK_REAL * restrict UpdateVar;
  CCTK_REAL * restrict OldVar;
  CCTK_REAL const * restrict RHSVar;
  CCTK_REAL * restrict ScratchVar;

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
     totalsize *= cctk_ash[arraydim];
  }

  if (scratchspace_firstindex == -99)
  {
    scratchspace_firstindex = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
  }
  
  if (scratchspace_firstindex_slow == -99)
  {
    scratchspace_firstindex_slow = CCTK_FirstVarIndex("MOL::SCRATCHSPACESLOW");
  }

  switch (MoL_Intermediate_Steps - (*MoL_Intermediate_Step))
  {
    case 0:
      alpha[0] = 1.0/4.0; 
      alpha[1] = 0;
      alpha[2] = 0;
      alpha[3] = 0;
      alpha[4] = 0;
      alpha[5] = 0;
      alpha[6] = 0;
      alpha[7] = 0;
      alpha[8] = 0;
      alpha_slow[0] = 1.0/4.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 0.0; 
      alpha_slow[5] = 0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[8] = 0;
      break;
    case 1:
      alpha[0] = -1.0/12.0;
      alpha[1] = 1.0/3.0;
      alpha[2] = 0;
      alpha[3] = 0;
      alpha[4] = 0.0;
      alpha[5] = 0.0;
      alpha[6] = 0;
      alpha[7] = 0;
      alpha[8] = 0;
      alpha_slow[0] = 1.0/4.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 0.0; 
      alpha_slow[5] = 0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[8] = 0;
      break;
    case 2:
      alpha[0] = 1.0/6.0; 
      alpha[1] = -1.0/6.0;
      alpha[2] = 1.0/2.0;
      alpha[3] = 0.0;
      alpha[4] = 0.0; 
      alpha[5] = 0.0;
      alpha[6] = 0.0;
      alpha[7] = 0.0;
      alpha[8] = 0.0;
      alpha_slow[0] = 1.0/2.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 0.0; 
      alpha_slow[5] = 0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[8] = 0;
      break;
    case 3:
      alpha[0] = 1.0/12.0;
      alpha[1] = 1.0/6.0; 
      alpha[2] = 1.0/6.0; 
      alpha[3] = 1.0/12.0;
      alpha[4] = 0.0; 
      alpha[5] = 0.0;
      alpha[6] = 0.0;
      alpha[7] = 0.0;
      alpha[8] = 0.0;
      alpha_slow[0] = 1.0/2.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 0.0; 
      alpha_slow[5] = 0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[8] = 0;
      break;
    case 4:
      alpha[0] = 1.0/12.0;
      alpha[1] = 1.0/6.0; 
      alpha[2] = 1.0/6.0; 
      alpha[3] = 1.0/6.0;
      alpha[4] = 1.0/12.0; 
      alpha[5] = 0.0;
      alpha[6] = 0.0;
      alpha[7] = 0.0;
      alpha[8] = 0.0;
      alpha_slow[0] = -1.0/6.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 2.0/3.0; 
      alpha_slow[5] = 0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[7] = 0;
      break;
    case 5:
      alpha[0] = 1.0/12.0;
      alpha[1] = 1.0/6.0; 
      alpha[2] = 1.0/6.0; 
      alpha[3] = 1.0/12.0;
      alpha[4] = 0.0; 
      alpha[5] = 1.0/4.0;
      alpha[6] = 0.0;
      alpha[7] = 0.0;
      alpha[8] = 0.0;
      alpha_slow[0] = 1.0/12.0;
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 1.0/6.0; 
      alpha_slow[5] = 1.0/2.0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[7] = 0;
      break;
    case 6:
      alpha[0] = 1.0/12.0;
      alpha[1] = 1.0/6.0; 
      alpha[2] = 1.0/6.0; 
      alpha[3] = 1.0/12.0;
      alpha[4] = 0.0; 
      alpha[5] = -1.0/12.0;
      alpha[6] = 1.0/3.0;
      alpha[7] = 0.0;
      alpha[8] = 0.0;
      alpha_slow[0] = 1.0/12.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = 1.0/6.0; 
      alpha_slow[5] = 1.0/2.0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[7] = 0;
      break;
    case 7:
      alpha[0] = 1.0/12.0;
      alpha[1] = 1.0/6.0; 
      alpha[2] = 1.0/6.0; 
      alpha[3] = 1.0/12.0;
      alpha[4] = 0.0; 
      alpha[5] = 1.0/6.0;
      alpha[6] = -1.0/6.0;
      alpha[7] = 1.0/2.0;
      alpha[8] = 0.0;
      alpha_slow[0] = 1.0/3.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = -1.0/3.0; 
      alpha_slow[5] = 1.0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[7] = 0;
      break;
    case 8:
      alpha[0] = 1.0/12.0;
      alpha[1] = 1.0/6.0; 
      alpha[2] = 1.0/6.0; 
      alpha[3] = 1.0/12.0;
      alpha[4] = 0.0; 
      alpha[5] = 1.0/12.0;
      alpha[6] = 1.0/6.0;
      alpha[7] = 1.0/6.0;
      alpha[8] = 1.0/12.0;
      alpha_slow[0] = 1.0/3.0; 
      alpha_slow[1] = 0;
      alpha_slow[2] = 0;
      alpha_slow[3] = 0;
      alpha_slow[4] = -1.0/3.0; 
      alpha_slow[5] = 1.0;
      alpha_slow[6] = 0;
      alpha_slow[7] = 0;
      alpha_slow[7] = 0;
      break;
    case 9:
      /* alpha is not used in this case */
      break;
    default:
      CCTK_VERROR("Internal error. Unexpected substep %d", (int)*MoL_Intermediate_Step);
      break;
  }

  beta[0] = 1.0/12.0;
  beta[1] = 1.0/6.0; 
  beta[2] = 1.0/6.0; 
  beta[3] = 1.0/12.0;
  beta[4] = 0.0;
  beta[5] = 1.0/12.0;
  beta[6] = 1.0/6.0; 
  beta[7] = 1.0/6.0; 
  beta[8] = 1.0/12.0;
  beta[9] = 0.0;

  beta_slow[0] = 1.0/6.0; 
  beta_slow[1] = 0.0; 
  beta_slow[2] = 0.0; 
  beta_slow[3] = 0.0;
  beta_slow[4] = 1.0/3.0;
  beta_slow[5] = 1.0/3.0; 
  beta_slow[6] = 0.0; 
  beta_slow[7] = 0.0; 
  beta_slow[8] = 0.0;
  beta_slow[9] = 1.0/6.0;

  /* FIXME */


  /* Real GFs, the "fast" part */

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndex[var]);
    OldVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 1, 
                                              EvolvedVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSVariableIndex[var]);
/* #define MOLDEBUG 1 */
#ifdef MOLDEBUG
    printf("In multirate RK. Variable %d (%s). RHS %d (%s). beta %g.\n",
           EvolvedVariableIndex[var],
           CCTK_VarName(EvolvedVariableIndex[var]),
           RHSVariableIndex[var],
           CCTK_VarName(RHSVariableIndex[var]),
           beta[MoL_Intermediate_Steps - (*MoL_Intermediate_Step)]);
#endif

    ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                      scratchspace_firstindex
                                      + (MoL_Intermediate_Steps - (*MoL_Intermediate_Step)));

#pragma omp /*parallel for*/ simd
    for (index = 0; index < totalsize; index++)
    {
      ScratchVar[index] = (*Original_Delta_Time) / cctkGH->cctk_timefac * RHSVar[index];
    
#ifdef MOLDEBUG
      if (CCTK_EQUALS(verbose,"extreme"))
      {
        printf("Variable: %d. Index: %d. dt: %f. beta %f. RHS: %f. q: %f.\n",
               var, index, (*Original_Delta_Time) / cctkGH->cctk_timefac, beta, RHSVar[index], 
               UpdateVar[index]);
      }
#endif
    }
    
    
#pragma omp /*parallel for*/ simd
    for (index = 0; index < totalsize; index++)
    {
      UpdateVar[index] = OldVar[index];
    }

    if ((*MoL_Intermediate_Step)>1)
    {
      //printf("Step %d \n", MoL_Intermediate_Steps - (*MoL_Intermediate_Step));
      for (scratchstep = 0; scratchstep <= MoL_Intermediate_Steps - (*MoL_Intermediate_Step); scratchstep++)
      {
         
         //printf("Scratch Step %d, alpha %g \n", scratchstep, alpha[scratchstep]);
         
         ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                       scratchspace_firstindex
                                       + scratchstep);
      
#pragma omp /*parallel for*/ simd
         for (index = 0; index < totalsize; index++)
         {
            UpdateVar[index] += alpha[scratchstep] * ScratchVar[index];
         }
      }
    }
    else
    {
      //printf("Final Step!\n");
    
      for (scratchstep = 0; scratchstep < MoL_Intermediate_Steps; scratchstep++)
      {
         
         //printf("Scratch Step %d, beta %g \n", scratchstep, beta[scratchstep]);
         
         ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                       scratchspace_firstindex
                                       + scratchstep);
      
#pragma omp /*parallel for*/ simd
         for (index = 0; index < totalsize; index++)
         {
            UpdateVar[index] += beta[scratchstep] * ScratchVar[index];
         }
      }
    }

  }


  for (var = 0; var < MoLNumEvolvedVariablesSlow; var++)
  {
    
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndexSlow[var]);
    OldVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 1, 
                                              EvolvedVariableIndexSlow[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
                                                 RHSVariableIndexSlow[var]);
/* #define MOLDEBUG 1 */
#ifdef MOLDEBUG
    printf("In multirate RK. SLOW Variable %d (%s). RHS %d (%s). beta %g.\n",
           EvolvedVariableIndexSlow[var],
           CCTK_VarName(EvolvedVariableIndexSlow[var]),
           RHSVariableIndexSlow[var],
           CCTK_VarName(RHSVariableIndexSlow[var]),
           beta[MoL_Intermediate_Steps - (*MoL_Intermediate_Step)]);
#endif

    ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                      scratchspace_firstindex_slow
                                      + (MoL_Intermediate_Steps - (*MoL_Intermediate_Step)));

#pragma omp /*parallel for*/ simd
    for (index = 0; index < totalsize; index++)
    {
      ScratchVar[index] = (*Original_Delta_Time) / cctkGH->cctk_timefac * RHSVar[index];
    
#ifdef MOLDEBUG
      if (CCTK_EQUALS(verbose,"extreme"))
      {
        printf("SLOW Variable: %d. Index: %d. dt: %f. beta %f. RHS: %f. q: %f.\n",
               var, index, (*Original_Delta_Time) / cctkGH->cctk_timefac, beta, RHSVar[index], 
               UpdateVar[index]);
      }
#endif
    }
    
    
#pragma omp /*parallel for*/ simd
    for (index = 0; index < totalsize; index++)
    {
      UpdateVar[index] = OldVar[index];
    }

    if ((*MoL_Intermediate_Step)>1)
    {
      //printf("Step %d \n", MoL_Intermediate_Steps - (*MoL_Intermediate_Step));
      for (scratchstep = 0; scratchstep <= MoL_Intermediate_Steps - (*MoL_Intermediate_Step); scratchstep++)
      {
         
         //printf("Scratch Step %d, alpha %g \n", scratchstep, alpha[scratchstep]);
         
         ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                       scratchspace_firstindex_slow
                                       + scratchstep);
      
#pragma omp /*parallel for*/ simd
         for (index = 0; index < totalsize; index++)
         {
            UpdateVar[index] += alpha_slow[scratchstep] * ScratchVar[index];
         }
      }
    }
    else
    {
      //printf("Final Step!\n");
    
      for (scratchstep = 0; scratchstep < MoL_Intermediate_Steps; scratchstep++)
      {
         
         //printf("Scratch Step %d, beta %g \n", scratchstep, beta[scratchstep]);
         
         ScratchVar = CCTK_VarDataPtrI(cctkGH, var, 
                                       scratchspace_firstindex_slow
                                       + scratchstep);
      
#pragma omp /*parallel for*/ simd
         for (index = 0; index < totalsize; index++)
         {
            UpdateVar[index] += beta_slow[scratchstep] * ScratchVar[index];
         }
      }
    }

  }

  return;
}
