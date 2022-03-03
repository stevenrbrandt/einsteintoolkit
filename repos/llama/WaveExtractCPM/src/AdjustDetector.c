 /*@@
   @file      AdjustDetector.c
   @date      Mon Nov 25 13:01:43 2002
   @author    Frank Herrmann
   @desc 
              Adjust next detector, increase the variable current_detector
   @enddesc 
   @version  $Id: AdjustDetector.c 6 2004-06-30 18:08:24Z herrmann $
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void WavExtrCPM_AdjustDetector(CCTK_ARGUMENTS);

 /*@@
   @routine    WavExtrCPM_AdjustDetector
   @date       Mon Nov 25 13:02:44 2002
   @author     Frank Herrmann
   @desc 
               Adjust next detector. This function decreases the variable 
               current_detector by one.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/
void WavExtrCPM_AdjustDetector(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose>4)
    CCTK_INFO("Adjust next detector");

  if (*do_nothing == 1)
    *do_nothing = 0;

  if (*current_detector != 0) 
    *current_detector=*current_detector-1;

  if (*current_detector < 0)
    *current_detector=0;

  if (verbose >2)
    CCTK_VInfo(CCTK_THORNSTRING,"Done with detector %d, next will go to Detector No. %d",*current_detector +1,*current_detector);
}
