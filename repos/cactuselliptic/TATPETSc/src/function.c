/* (C) 2001-04-18 Erik Schnetter <schnetter@uni-tuebingen.de> */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <petscdmda.h>
#include <petscsnes.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "TATPETSc.h"



int TATPETSc_function (SNES snes, Vec x, Vec f, void *userptr)
{
  DECLARE_CCTK_PARAMETERS;
  
#if 1
  userdata *user = (userdata*)userptr;
#else
  DMMG dmmg = (DMMG)userptr;
  userdata *user = (userdata*)dmmg->user;
#endif
  const cGH *cctkGH = user->cctkGH;
  
  assert (user->magic==MAGIC);
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_function");
  
  ++user->funcall_count;
  
  TATPETSc_copy (x, userptr, TATcopyout, TATcopyvars);
  int ierr = (user->bnd) (cctkGH, -1, user->data);
  int ierr2 = (user->fun) (cctkGH, -1, user->data);
  if (!ierr) ierr = ierr2;
  TATPETSc_copy (f, userptr, TATcopyin, TATcopyvals);
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_function done.");
  
  return ierr;
}
