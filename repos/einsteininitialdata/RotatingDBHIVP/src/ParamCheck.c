 /*@@
   @file      ParamCheck.c
   @date      Fri Apr 26 18:03:09 2002
   @author    Tom Goodale
   @desc 
              Check the parameters for RotatingDBHIVP
   @enddesc 
   @version $Header: /arrangements/PerturbedBH2/RotatingDBHIVP/src/ParamCheck.c,v 1.1 2002/04/30 guzman
 @@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(RotatingDBHIVP_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void RotatingDBHIVP_ParamCheck(CCTK_ARGUMENTS);
void RotatingDBHIVP_ConformalCheck(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

void RotatingDBHIVP_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* Do we know how to deal with this type of metric ? */
  if(! CCTK_EQUALS(metric_type, "physical") &&
     ! CCTK_EQUALS(metric_type, "static conformal"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
  }

  if(CCTK_EQUALS(metric_type, "static conformal") &&
     ! CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs"))
  {
    CCTK_PARAMWARN("ADMConstraints can currently only work with a physical metric or a static conf");
  }
   
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    RotatingDBHIVP_ConformalCheck
   @date       Tu Apr 30 2002
   @author     Goodale-Guzman
   @desc
   Check that the initial data has setup enough derivatives.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

 @@*/
void RotatingDBHIVP_ConformalCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_EQUALS(metric_type, "physical"))
  {
    /* Use this internally so we don't need to check the parameter again. */
    *conformal_state = 0;
  }
  else if(CCTK_EQUALS(metric_type, "static conformal"))
  {
    if(*conformal_state < 3)
    {
      CCTK_WARN(0, "RotatingDBHIVP needs second derivatives of the conformal factor when running with a static conformal metric");
    }
  }

}
/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

