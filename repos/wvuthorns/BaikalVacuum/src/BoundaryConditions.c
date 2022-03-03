
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"

// Set `none` boundary conditions on BSSN RHSs, as these are set via NewRad.
void BaikalVacuum_BoundaryConditions_evolved_gfs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::aDD00GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::aDD01GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::aDD02GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::aDD11GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::aDD12GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::aDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::aDD22GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::alphaGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::alphaGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::betU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::betU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::betU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::betU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::betU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::betU2GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::cfGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::cfGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD00GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::hDD00GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD01GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::hDD01GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD02GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::hDD02GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD11GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::hDD11GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD12GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::hDD12GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::hDD22GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::hDD22GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::lambdaU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::lambdaU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::lambdaU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::lambdaU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::lambdaU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::lambdaU2GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::trKGF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::trKGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::vetU0GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::vetU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::vetU1GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::vetU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalVacuum::vetU2GF", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::vetU2GF!");

}

// Set `flat` boundary conditions on BSSN constraints, similar to what Lean does.
void BaikalVacuum_BoundaryConditions_aux_gfs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  CCTK_INT bndsize;  // <- bndsize = number of ghostzones
  if(FD_order == 4) {
    bndsize = 3;
  } else if(FD_order == 6) {
    bndsize = 4;
  } else if(FD_order == 8) {
    bndsize = 5;
  } else {
    CCTK_ERROR("Error: Inside BaikalVacuum_BoundaryConditions_aux_gfs(): Chose unsupported FD_order.");
  }

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::HGF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::HGF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::MU0GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::MU0GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::MU1GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::MU1GF!");

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize, -1, "BaikalVacuum::MU2GF", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalVacuum::MU2GF!");
}
