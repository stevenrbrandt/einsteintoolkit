/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


extern "C" void EM_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void EM_SelectBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(El_group_bound, "none"  ) ||
      CCTK_EQUALS(El_group_bound, "static") ||
      CCTK_EQUALS(El_group_bound, "flat"  ) ||
      CCTK_EQUALS(El_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::El_group", El_group_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register El_group_bound BC for EM::El_group!");
  }
  
  if (CCTK_EQUALS(B_group_bound, "none"  ) ||
      CCTK_EQUALS(B_group_bound, "static") ||
      CCTK_EQUALS(B_group_bound, "flat"  ) ||
      CCTK_EQUALS(B_group_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::B_group", B_group_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register B_group_bound BC for EM::B_group!");
  }
  
  if (CCTK_EQUALS(El1_bound, "none"  ) ||
      CCTK_EQUALS(El1_bound, "static") ||
      CCTK_EQUALS(El1_bound, "flat"  ) ||
      CCTK_EQUALS(El1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::El1", El1_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register El1_bound BC for EM::El1!");
  }
  
  if (CCTK_EQUALS(El2_bound, "none"  ) ||
      CCTK_EQUALS(El2_bound, "static") ||
      CCTK_EQUALS(El2_bound, "flat"  ) ||
      CCTK_EQUALS(El2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::El2", El2_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register El2_bound BC for EM::El2!");
  }
  
  if (CCTK_EQUALS(El3_bound, "none"  ) ||
      CCTK_EQUALS(El3_bound, "static") ||
      CCTK_EQUALS(El3_bound, "flat"  ) ||
      CCTK_EQUALS(El3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::El3", El3_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register El3_bound BC for EM::El3!");
  }
  
  if (CCTK_EQUALS(B1_bound, "none"  ) ||
      CCTK_EQUALS(B1_bound, "static") ||
      CCTK_EQUALS(B1_bound, "flat"  ) ||
      CCTK_EQUALS(B1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::B1", B1_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register B1_bound BC for EM::B1!");
  }
  
  if (CCTK_EQUALS(B2_bound, "none"  ) ||
      CCTK_EQUALS(B2_bound, "static") ||
      CCTK_EQUALS(B2_bound, "flat"  ) ||
      CCTK_EQUALS(B2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::B2", B2_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register B2_bound BC for EM::B2!");
  }
  
  if (CCTK_EQUALS(B3_bound, "none"  ) ||
      CCTK_EQUALS(B3_bound, "static") ||
      CCTK_EQUALS(B3_bound, "flat"  ) ||
      CCTK_EQUALS(B3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "EM::B3", B3_bound);
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register B3_bound BC for EM::B3!");
  }
  
  if (CCTK_EQUALS(El_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_El_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El_group_bound < 0) handle_El_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El_group_bound , El_group_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_El_group_bound ,El_group_bound_speed, "SPEED") < 0)
       CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El_group_bound, 
                      "EM::El_group", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::El_group!");
  
  }
  
  if (CCTK_EQUALS(B_group_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B_group_bound < 0) handle_B_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B_group_bound , B_group_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B_group_bound ,B_group_bound_speed, "SPEED") < 0)
       CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B_group_bound, 
                      "EM::B_group", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::B_group!");
  
  }
  
  if (CCTK_EQUALS(El1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_El1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El1_bound < 0) handle_El1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El1_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El1_bound , El1_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_El1_bound ,El1_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El1_bound, 
                      "EM::El1", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::El1!");
  
  }
  
  if (CCTK_EQUALS(El2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_El2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El2_bound < 0) handle_El2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El2_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El2_bound , El2_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_El2_bound ,El2_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El2_bound, 
                      "EM::El2", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::El2!");
  
  }
  
  if (CCTK_EQUALS(El3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_El3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El3_bound < 0) handle_El3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El3_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El3_bound , El3_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_El3_bound ,El3_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El3_bound, 
                      "EM::El3", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::El3!");
  
  }
  
  if (CCTK_EQUALS(B1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B1_bound < 0) handle_B1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B1_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B1_bound , B1_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B1_bound ,B1_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B1_bound, 
                      "EM::B1", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::B1!");
  
  }
  
  if (CCTK_EQUALS(B2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B2_bound < 0) handle_B2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B2_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B2_bound , B2_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B2_bound ,B2_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B2_bound, 
                      "EM::B2", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::B2!");
  
  }
  
  if (CCTK_EQUALS(B3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B3_bound < 0) handle_B3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B3_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B3_bound , B3_bound_limit, "LIMIT") < 0)
       CCTK_WARN(0, "could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B3_bound ,B3_bound_speed, "SPEED") < 0)
        CCTK_WARN(0, "could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B3_bound, 
                      "EM::B3", "Radiation");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Radiation BC for EM::B3!");
  
  }
  
  if (CCTK_EQUALS(El_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_El_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El_group_bound < 0) handle_El_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El_group_bound ,El_group_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El_group_bound, 
                      "EM::El_group", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Scalar BC for EM::El_group!");
  
  }
  
  if (CCTK_EQUALS(B_group_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B_group_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B_group_bound < 0) handle_B_group_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B_group_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B_group_bound ,B_group_bound_scalar, "SCALAR") < 0)
        CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B_group_bound, 
                      "EM::B_group", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Failed to register Scalar BC for EM::B_group!");
  
  }
  
  if (CCTK_EQUALS(El1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_El1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El1_bound < 0) handle_El1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El1_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El1_bound ,El1_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El1_bound, 
                      "EM::El1", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for EM::El1!");
  
  }
  
  if (CCTK_EQUALS(El2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_El2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El2_bound < 0) handle_El2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El2_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El2_bound ,El2_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El2_bound, 
                      "EM::El2", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for EM::El2!");
  
  }
  
  if (CCTK_EQUALS(El3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_El3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_El3_bound < 0) handle_El3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_El3_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_El3_bound ,El3_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_El3_bound, 
                      "EM::El3", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for EM::El3!");
  
  }
  
  if (CCTK_EQUALS(B1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B1_bound < 0) handle_B1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B1_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B1_bound ,B1_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B1_bound, 
                      "EM::B1", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for EM::B1!");
  
  }
  
  if (CCTK_EQUALS(B2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B2_bound < 0) handle_B2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B2_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B2_bound ,B2_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B2_bound, 
                      "EM::B2", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for EM::B2!");
  
  }
  
  if (CCTK_EQUALS(B3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B3_bound < 0) handle_B3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B3_bound < 0) CCTK_WARN(0, "could not create table!");
    if (Util_TableSetReal(handle_B3_bound ,B3_bound_scalar, "SCALAR") < 0)
      CCTK_WARN(0, "could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B3_bound, 
                      "EM::B3", "scalar");
  
    if (ierr < 0)
       CCTK_WARN(0, "Error in registering Scalar BC for EM::B3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#EM::El_group_bound       = "skip"
#$bound$#EM::El_group_bound_speed = 1.0
#$bound$#EM::El_group_bound_limit = 0.0
#$bound$#EM::El_group_bound_scalar = 0.0

#$bound$#EM::B_group_bound       = "skip"
#$bound$#EM::B_group_bound_speed = 1.0
#$bound$#EM::B_group_bound_limit = 0.0
#$bound$#EM::B_group_bound_scalar = 0.0

#$bound$#EM::El1_bound       = "skip"
#$bound$#EM::El1_bound_speed = 1.0
#$bound$#EM::El1_bound_limit = 0.0
#$bound$#EM::El1_bound_scalar = 0.0

#$bound$#EM::El2_bound       = "skip"
#$bound$#EM::El2_bound_speed = 1.0
#$bound$#EM::El2_bound_limit = 0.0
#$bound$#EM::El2_bound_scalar = 0.0

#$bound$#EM::El3_bound       = "skip"
#$bound$#EM::El3_bound_speed = 1.0
#$bound$#EM::El3_bound_limit = 0.0
#$bound$#EM::El3_bound_scalar = 0.0

#$bound$#EM::B1_bound       = "skip"
#$bound$#EM::B1_bound_speed = 1.0
#$bound$#EM::B1_bound_limit = 0.0
#$bound$#EM::B1_bound_scalar = 0.0

#$bound$#EM::B2_bound       = "skip"
#$bound$#EM::B2_bound_speed = 1.0
#$bound$#EM::B2_bound_limit = 0.0
#$bound$#EM::B2_bound_scalar = 0.0

#$bound$#EM::B3_bound       = "skip"
#$bound$#EM::B3_bound_speed = 1.0
#$bound$#EM::B3_bound_limit = 0.0
#$bound$#EM::B3_bound_scalar = 0.0

*/

