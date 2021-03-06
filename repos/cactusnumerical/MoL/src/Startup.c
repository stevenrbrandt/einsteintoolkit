 /*@@
   @file      Startup.c
   @date      Wed May 22 02:09:19 2002
   @author    
   @desc 
   Register the startup banner.
   The external variables are also declared here. These are 
   the arrays containing the variable indexes and right hand sides, 
   and the number of each type of variable currently in use (the 
   parameters only give the maximum possible for each type).
   @enddesc 
   @version   $Header$
 @@*/

#include <stddef.h>

#include "cctk.h"
#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_Startup_c);

/********************************************************************
 ********************     External Variables   **********************
 ********************************************************************/

CCTK_INT *restrict EvolvedVariableIndex = NULL;
CCTK_INT *restrict EvolvedVariableIndexSlow = NULL;
CCTK_INT *restrict RHSVariableIndex = NULL;
CCTK_INT *restrict RHSVariableIndexSlow = NULL;
CCTK_INT *restrict ConstrainedVariableIndex = NULL;
CCTK_INT *restrict SandRVariableIndex = NULL;

CCTK_INT MoLNumEvolvedVariables = 0;
CCTK_INT MoLNumEvolvedVariablesSlow = 0;
CCTK_INT MoLNumConstrainedVariables = 0;
CCTK_INT MoLNumSandRVariables = 0;


CCTK_INT *restrict EvolvedArrayVariableIndex = NULL;
CCTK_INT *restrict RHSArrayVariableIndex = NULL;
CCTK_INT *restrict ConstrainedArrayVariableIndex = NULL;
CCTK_INT *restrict SandRArrayVariableIndex = NULL;

CCTK_INT MoLNumEvolvedArrayVariables = 0;
CCTK_INT MoLNumConstrainedArrayVariables = 0;
CCTK_INT MoLNumSandRArrayVariables = 0;


CCTK_INT ScheduleStatus = 0;

CCTK_REAL *restrict ArrayScratchSpace = NULL;
CCTK_INT *restrict ArrayScratchSizes = NULL;
CCTK_INT CurrentArrayScratchSize = 0;

CCTK_REAL *restrict ArrayErrorSpace = NULL;
CCTK_INT *restrict ArrayErrorSizes = NULL;
CCTK_INT CurrentArrayErrorSize = 0;

CCTK_INT MoLMaxNumRegisteredVariables = 0;

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

int MoL_Startup(void);

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
   @routine    MoL_Startup
   @date       Wed May 22 02:17:17 2002
   @author     Ian Hawke
   @desc 
   Register the startup banner with the flesh.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int MoL_Startup(void) 
{
  
  const char *banner = "MoL: Generalized time integration.";
  
  CCTK_RegisterBanner(banner);
  
  return 0;
  
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
