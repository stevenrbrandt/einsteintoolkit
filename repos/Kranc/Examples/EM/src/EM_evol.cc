/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"

namespace EM {

extern "C" void EM_evol_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % EM_evol_calc_every != EM_evol_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EM::B_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for EM::B_grouprhs.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "EM::El_grouprhs","flat");
  if (ierr < 0)
    CCTK_WARN(1, "Failed to register flat BC for EM::El_grouprhs.");
  return;
}

static void EM_evol_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3(EM_evol,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED = B1[index];
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED = B2[index];
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED = B3[index];
    CCTK_REAL El1L CCTK_ATTRIBUTE_UNUSED = El1[index];
    CCTK_REAL El2L CCTK_ATTRIBUTE_UNUSED = El2[index];
    CCTK_REAL El3L CCTK_ATTRIBUTE_UNUSED = El3[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    const CCTK_REAL PDstandard2nd2B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B1[index]);
    const CCTK_REAL PDstandard2nd3B1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B1[index]);
    const CCTK_REAL PDstandard2nd1B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B2[index]);
    const CCTK_REAL PDstandard2nd3B2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&B2[index]);
    const CCTK_REAL PDstandard2nd1B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&B3[index]);
    const CCTK_REAL PDstandard2nd2B3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&B3[index]);
    const CCTK_REAL PDstandard2nd2El1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&El1[index]);
    const CCTK_REAL PDstandard2nd3El1 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&El1[index]);
    const CCTK_REAL PDstandard2nd1El2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&El2[index]);
    const CCTK_REAL PDstandard2nd3El2 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd3(&El2[index]);
    const CCTK_REAL PDstandard2nd1El3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1(&El3[index]);
    const CCTK_REAL PDstandard2nd2El3 CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2(&El3[index]);
    /* Calculate temporaries and grid functions */
    CCTK_REAL El1rhsL CCTK_ATTRIBUTE_UNUSED = PDstandard2nd2B3 - 
      PDstandard2nd3B2;
    
    CCTK_REAL El2rhsL CCTK_ATTRIBUTE_UNUSED = -PDstandard2nd1B3 + 
      PDstandard2nd3B1;
    
    CCTK_REAL El3rhsL CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1B2 - 
      PDstandard2nd2B1;
    
    CCTK_REAL B1rhsL CCTK_ATTRIBUTE_UNUSED = -PDstandard2nd2El3 + 
      PDstandard2nd3El2;
    
    CCTK_REAL B2rhsL CCTK_ATTRIBUTE_UNUSED = PDstandard2nd1El3 - 
      PDstandard2nd3El1;
    
    CCTK_REAL B3rhsL CCTK_ATTRIBUTE_UNUSED = -PDstandard2nd1El2 + 
      PDstandard2nd2El1;
    /* Copy local copies back to grid functions */
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    El1rhs[index] = El1rhsL;
    El2rhs[index] = El2rhsL;
    El3rhs[index] = El3rhsL;
  }
  CCTK_ENDLOOP3(EM_evol);
}
extern "C" void EM_evol(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering EM_evol_Body");
  }
  if (cctk_iteration % EM_evol_calc_every != EM_evol_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "EM::B_group",
    "EM::B_grouprhs",
    "EM::El_group",
    "EM::El_grouprhs"};
  AssertGroupStorage(cctkGH, "EM_evol", 4, groups);
  
  EnsureStencilFits(cctkGH, "EM_evol", 1, 1, 1);
  
  LoopOverInterior(cctkGH, EM_evol_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving EM_evol_Body");
  }
}

} // namespace EM
