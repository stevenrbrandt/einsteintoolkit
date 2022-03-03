/* (C) 2001-04-18 Erik Schnetter <schnetter@uni-tuebingen.de> */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>

#include <mpi.h>

#include <petscdmda.h>
#include <petscsnes.h>
#include <petscversion.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <util_ErrorCodes.h>
#include <util_Table.h>

#include "TATPETSc.h"



int TATPETSc_solve (const cGH *cctkGH,
		    const int *var, const int *val, int nvars,
		    int options_table,
		    int (*fun) (const cGH *cctkGH, int options_table,
                                void *data),
		    int (*bnd) (const cGH *cctkGH, int options_table,
                                void *data),
		    void *data)
{
  DECLARE_CCTK_PARAMETERS;
  
  PetscInt ierr;
  
  
  
  if (! CCTK_IsThornActive(CCTK_THORNSTRING)) {
    CCTK_WARN (0, "Thorn " CCTK_THORNSTRING " has not been activated.  It is therefore not possible to call TATPETSc_solve.");
  }
  
  
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_solve");
  
  
  
  /* application data */
  userdata user;
  user.magic = MAGIC;
  user.cctkGH = cctkGH;
  
  /* world communicator */
  MPI_Comm comm = PETSC_COMM_WORLD;
  user.comm = comm;
  int rank, size;
  MPI_Comm_size (comm, &size);
  MPI_Comm_rank (comm, &rank);
  
  
  
  /* Check and copy arguments */
  
  assert (nvars>0);
  assert (var);
  assert (val);
  user.nvars = nvars;
  user.var = malloc(sizeof(*user.var) * user.nvars);
  assert (user.var);
  user.val = malloc(sizeof(*user.val) * user.nvars);
  assert (user.val);
  for (int n=0; n<nvars; ++n) {
    assert (var[n]>=0 && var[n]<CCTK_NumVars());
    user.var[n] = var[n];
    assert (val[n]>=0 && val[n]<CCTK_NumVars());
    user.val[n] = val[n];
  }
  
  cGroupDynamicData dyndata;
  {
    int group = CCTK_GroupIndexFromVarI (user.var[0]);
    assert (group>=0);
    CCTK_GroupDynamicData (cctkGH, group, &user.dyndata);
  }
  
  assert (user.dyndata.dim <= DIM);
  
  /* Check grid variable properties */
  for (int n=0; n<nvars; ++n) {
    int group = CCTK_GroupIndexFromVarI (user.var[n]);
    assert (group>=0);
    CCTK_GroupDynamicData (cctkGH, group, &dyndata);
    assert (dyndata.dim == user.dyndata.dim);
    for (int d=0; d<user.dyndata.dim; ++d) {
      assert (dyndata.dim == user.dyndata.dim);
      assert (dyndata.gsh[d] == user.dyndata.gsh[d]);
      assert (dyndata.lsh[d] == user.dyndata.lsh[d]);
      assert (dyndata.lbnd[d] == user.dyndata.lbnd[d]);
      assert (dyndata.ubnd[d] == user.dyndata.ubnd[d]);
      assert (dyndata.bbox[2*d] == user.dyndata.bbox[2*d]);
      assert (dyndata.bbox[2*d+1] == user.dyndata.bbox[2*d+1]);
      assert (dyndata.nghostzones[d] == user.dyndata.nghostzones[d]);
    }
    group = CCTK_GroupIndexFromVarI (user.val[n]);
    assert (group>=0);
    CCTK_GroupDynamicData (cctkGH, group, &dyndata);
    assert (dyndata.dim == user.dyndata.dim);
    for (int d=0; d<user.dyndata.dim; ++d) {
      assert (dyndata.dim == user.dyndata.dim);
      assert (dyndata.gsh[d] == user.dyndata.gsh[d]);
      assert (dyndata.lsh[d] == user.dyndata.lsh[d]);
      assert (dyndata.lbnd[d] == user.dyndata.lbnd[d]);
      assert (dyndata.ubnd[d] == user.dyndata.ubnd[d]);
      assert (dyndata.bbox[2*d] == user.dyndata.bbox[2*d]);
      assert (dyndata.bbox[2*d+1] == user.dyndata.bbox[2*d+1]);
      assert (dyndata.nghostzones[d] == user.dyndata.nghostzones[d]);
    }
  }
  
  const int dim = user.dyndata.dim;
  assert (dim>=0 && dim<=DIM);
  
  /* optional arguments */
  CCTK_INT periodic[DIM];
  ierr = Util_TableGetIntArray (options_table, DIM, periodic, "periodic");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (int d=0; d<dim; ++d) periodic[d] = 0;
    ierr = dim;
  }
  assert (ierr == dim);
  
  CCTK_INT solvebnds[2*DIM];
  ierr = Util_TableGetIntArray (options_table, 2*DIM, solvebnds, "solvebnds");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (int d=0; d<2*dim; ++d) solvebnds[d] = 0;
    ierr = 2*dim;
  }
  assert (ierr == 2*dim);
  
  CCTK_INT nboundaryzones[2*DIM];
  ierr = Util_TableGetIntArray (options_table, 2*DIM, nboundaryzones,
                                "nboundaryzones");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    for (int d=0; d<2*dim; ++d) {
      nboundaryzones[d] = user.dyndata.nghostzones[d/2];
    }
    ierr = 2*dim;
  }
  assert (ierr == 2*dim);
  
  CCTK_INT sw;
  ierr = Util_TableGetInt (options_table, &sw, "stencil_width");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    sw = 1;
    ierr = 1;
  }
  assert (ierr == 1);
  
  CCTK_FPOINTER ptmp;
  ierr = Util_TableGetFnPointer (options_table, &ptmp, "jacobian");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    ptmp = NULL;
    ierr = 1;
  }
  assert (ierr == 1);
  int (*jac) (const cGH *cctkGH, int options_table, void *data);
  jac = ptmp;
  
  CCTK_INT * restrict jac0;
  if (jac) {
    jac0 = malloc (nvars * sizeof *jac0);
    assert (jac0);
    ierr = Util_TableGetIntArray (options_table, nvars, jac0, "jac");
    assert (ierr == nvars);
  } else {
    jac0 = NULL;
  }
  
  ierr = Util_TableGetFnPointer (options_table, &ptmp, "get_coloring");
  if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    ptmp = 0;
    ierr = 1;
  }
  assert (ierr == 1);
  int (*get_coloring) (DM da, ISColoring *iscoloring, Mat *J, void *data);
  get_coloring = ptmp;
  
  for (int d=0; d<2*user.dyndata.dim; ++d) {
    user.nboundaryzones[d] = solvebnds[d] ? 0 : nboundaryzones[d];
  }
  
  assert (fun);
  user.fun = fun;
  assert (bnd);
  user.bnd = bnd;
  user.jac = jac;
  user.jac0 = jac0;
  user.data = data;
  user.funcall_count = 0;
  user.jaccall_count = 0;
  
  
  
  assert (sw>=0);
#if 0
  assert (mglevels>=2);
#endif
  
  
  
  /* Deinstall error handler temporarily */
  ierr = PetscPopErrorHandler ();
  CHKERRQ(ierr);
  
  
  
  /* Determine Cactus' processor layout of grid variables */
  
  const int *lbnd = user.dyndata.lbnd;
  const int *lsh  = user.dyndata.lsh;
  
  const int nprocs = CCTK_nProcs(cctkGH);
  assert (nprocs>0);
  
  /* layout: [nprocs][dim] as [nprocs*dim + dim] */
  int *all_lbnd = malloc(nprocs * dim * sizeof(*all_lbnd));
  assert (all_lbnd);
  MPI_Allgather ((void*)lbnd, dim, MPI_INT, all_lbnd, dim, MPI_INT, comm);
  int *all_lsh = malloc(nprocs * dim * sizeof(*all_lsh));
  assert (all_lsh);
  MPI_Allgather ((void*)lsh, dim, MPI_INT, all_lsh, dim, MPI_INT, comm);
  
  /* number of processors per direction */
  int nprocs_dim[DIM];
  for (int d=0; d<dim; ++d) {
    nprocs_dim[d] = 1;
    for (int n=1; n<nprocs; ++n) {
      if (all_lbnd[n*dim+d] > all_lbnd[(n-1)*dim+d]) ++nprocs_dim[d];
      if (all_lbnd[n*dim+d] < all_lbnd[(n-1)*dim+d]) break;
    }
  }
  
  /* product_prefix of nprocs_dim */
  PetscInt proc_offset_dim[DIM]; /* [dim] */
  if (dim>0) {
    proc_offset_dim[0] = 1;
    /* assume x loops fastest */
    for (int d=1; d<dim; ++d) {
      proc_offset_dim[d] = proc_offset_dim[d-1] * nprocs_dim[d-1];
    }
  }
  
  /* number of grid points per direction and per processor */
  PetscInt *npoints_proc[DIM];	/* [dim][nprocs_dim[dim]] */
  for (int d=0; d<dim; ++d) {
    npoints_proc[d] = malloc(nprocs_dim[d] * sizeof(*npoints_proc[d]));
    assert (npoints_proc[d]);
    PetscInt sum = 0;
    for (int n=0; n<nprocs_dim[d]; ++n) {
      npoints_proc[d][n] = all_lsh[(proc_offset_dim[d]*n)*dim+d];
      sum += npoints_proc[d][n] - 2*user.dyndata.nghostzones[d];
    }
    assert (sum == user.dyndata.gsh[d] - 2*user.dyndata.nghostzones[d]);
  }
  
  /* check all assumptions (for the current processor) */
  for (int d=0; d<dim; ++d) {
    PetscInt n = (CCTK_MyProc(cctkGH) / proc_offset_dim[d]) % nprocs_dim[d];
    assert (npoints_proc[d][n] == lsh[d]);
  }
  
  
  
  /* Create distributed array object */
  
  /* matrix */
  PetscInt NI[DIM];             /* number of grid points */
  PetscInt *lx[DIM];            /* local number of grid points */
  Mat J;			/* Jacobian */
  
  for (int d=0; d<user.dyndata.dim; ++d) {
    /* global number of grid points */
    NI[d] = user.dyndata.gsh[d] - (user.nboundaryzones[2*d]
                                   + user.nboundaryzones[2*d+1]);
    /* local number of grid points */
    lx[d] = malloc(nprocs_dim[d] * sizeof(*lx[d]));
    PetscInt sum = 0;
    for (int n=0; n<nprocs_dim[d]; ++n) {
      lx[d][n] = (npoints_proc[d][n]
                  - (n==0
                     ? user.nboundaryzones[2*d  ]
                     : user.dyndata.nghostzones[d])
                  - (n==nprocs_dim[d]-1
                     ? user.nboundaryzones[2*d+1]
                     : user.dyndata.nghostzones[d]));
      sum += lx[d][n];
    }
    assert (sum == NI[d]);
  }
  
  /* distributed array */
  DM da;
  switch (user.dyndata.dim) {
  case 1:
    if (veryverbose) CCTK_INFO ("DMDACreate1d");
    ierr = DMDACreate1d (comm,
                         periodic[0] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
                         NI[0],
                         nvars, /* degrees of freedom */
                         sw,    /* stencil width */
                         lx[0],
                         &da);
    CHKERRQ(ierr);
    ierr = DMSetUp(da);
    CHKERRQ(ierr);
    break;
  case 2:
    if (veryverbose) CCTK_INFO ("DMDACreate2d");
    ierr = DMDACreate2d (comm,
                         periodic[0] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
                         periodic[1] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
                         DMDA_STENCIL_BOX,
                         NI[0],NI[1],
                         nprocs_dim[0],nprocs_dim[1],
                         nvars, /* degrees of freedom */
                         sw,    /* stencil width */
                         lx[0],lx[1],
                         &da);
    CHKERRQ(ierr);
    ierr = DMSetUp(da);
    CHKERRQ(ierr);
    break;
  case 3:
    if (veryverbose) CCTK_INFO ("DMDACreate3d");
    ierr = DMDACreate3d (comm,
                         periodic[0] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
                         periodic[1] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
                         periodic[2] ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE,
                         DMDA_STENCIL_BOX,
                         NI[0],NI[1],NI[2],
                         nprocs_dim[0],nprocs_dim[1],nprocs_dim[2],
                         nvars, /* degrees of freedom */
                         sw,    /* stencil width */
                         lx[0],lx[1],lx[2],
                         &da);
    CHKERRQ(ierr);
    ierr = DMSetUp(da);
    CHKERRQ(ierr);
    break;
  default:
    CCTK_VERROR ("Unsupported dimensionality: %d", (int)user.dyndata.dim);
  }
  user.da = da;
  
  for (int n=0; n<nvars; ++n) {
    ierr = DMDASetFieldName (da, n, CCTK_VarName(var[n]));
    CHKERRQ(ierr);
  }
  
  
  
#if 1
  
  
  
  /* Extract global and local vectors from DMDA */
  
  if (veryverbose) CCTK_INFO ("DMCreateGlobalVector");
  Vec x;			/* solution */
  ierr = DMCreateGlobalVector (da, &x);
  CHKERRQ(ierr);
  Vec f;			/* function value */
  ierr = VecDuplicate (x, &f);
  CHKERRQ(ierr);
  
  
  
  /* Create nonlinear solver context */
  
  if (veryverbose) CCTK_INFO ("SNESCreate");
  SNES snes;
  ierr = SNESCreate (comm, &snes);
  CHKERRQ(ierr);
  
  
  
  /* Set function evaluation routine and vector */
  
  if (veryverbose) CCTK_INFO ("SNESSetFunction");
  ierr = SNESSetFunction (snes, f, TATPETSc_function, &user);
  CHKERRQ(ierr);
  
  
  
  /* matrix coloring */
  ISColoring iscoloring;
  MatFDColoring matfdcoloring;
  
  if (jac) {
    /* Calculate Jacobian directly through a user given function */
    
    if (veryverbose) CCTK_INFO ("DMCreateMatrix");
    ierr = DMCreateMatrix (da, &J);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("DMCreateColoring");
    ierr = DMCreateColoring (da, IS_COLORING_GLOBAL, &iscoloring);
    CHKERRQ(ierr);
    ierr = TATPETSc_jacobian (snes, x, J, J, &user);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("SNESSetJacobian");
    ierr = SNESSetJacobian (snes, J, J, TATPETSc_jacobian, &user);
    CHKERRQ(ierr);
    
  } else {
    /* Approximate Jacobian numerically (automatically) */
    
    if (!get_coloring) {
      if (veryverbose) CCTK_INFO ("DMCreateMatrix");
      ierr = DMCreateMatrix (da, &J);
      CHKERRQ(ierr);
      if (veryverbose) CCTK_INFO ("DMCreateColoring");
      ierr = DMCreateColoring (da, IS_COLORING_GLOBAL, &iscoloring);
      CHKERRQ(ierr);
    } else {
      if (veryverbose) CCTK_INFO ("get_coloring");
      ierr = get_coloring (da, &iscoloring, &J, data);
      CHKERRQ(ierr);
    }
    if (veryverbose) CCTK_INFO ("MatFDColoringCreate");
    ierr = MatFDColoringCreate (J, iscoloring, &matfdcoloring);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("MatFDColoringSetFunction");
    ierr = MatFDColoringSetFunction (matfdcoloring,
				     (int(*)(void))TATPETSc_function, &user);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("MatFDColoringSetFromOptions");
    ierr = MatFDColoringSetFromOptions (matfdcoloring);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("MatFDColoringSetUp");
    ierr = MatFDColoringSetUp (J, iscoloring, matfdcoloring);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("ISColoringDestroy");
    ierr = ISColoringDestroy (&iscoloring);
    CHKERRQ(ierr);
    if (veryverbose) CCTK_INFO ("SNESSetJacobian");
    ierr = SNESSetJacobian (snes, J, J, SNESComputeJacobianDefaultColor,
			    matfdcoloring);
    CHKERRQ(ierr);
    
  }
  
  
  
  /* Customise solver */
  
  if (veryverbose) CCTK_INFO ("SNESSetFromOptions");
  ierr = SNESSetFromOptions (snes);
  CHKERRQ(ierr);
  
  
  
  /* Make initial guess */
  
  ierr = TATPETSc_copy (x, &user, TATcopyin, TATcopyvars);
  CHKERRQ(ierr);
  
  
  
#else
  
  
  
  /* Create multigrid solver */
  
  if (veryverbose) CCTK_INFO ("DMMGCreate");
  DMMG *dmmg;
  ierr =  DMMGCreate (PETSC_COMM_WORLD, mglevels, &user, &dmmg);
  CHKERRQ(ierr);
  
  if (veryverbose) CCTK_INFO ("DMMGSetDM");
  ierr = DMMGSetDM (dmmg, (DM)da);
  CHKERRQ(ierr);
  
  if (veryverbose) CCTK_INFO ("DMMGSetSNES");
  ierr = DMMGSetSNES (dmmg, TATPETSc_function, PETSC_NULL);
  CHKERRQ(ierr);
  
  if (veryverbose) CCTK_INFO ("DMMGGetSNES");
  snes = DMMGGetSNES (dmmg);
  assert (snes);
  
  Vec x;			/* solution */
  if (veryverbose) CCTK_INFO ("SNESGetSolution");
  ierr = SNESGetSolution (snes, &x);
  CHKERRQ(ierr);
  assert (x);
  
  Vec f;			/* function value */
  if (veryverbose) CCTK_INFO ("SNESGetFunction");
  ierr = SNESGetFunction (snes, &f, 0, 0);
  CHKERRQ(ierr);
  assert (f);
  
#if 0
  if (veryverbose) CCTK_INFO ("DMMGSetInitialGuess");
  ierr = DMMGSetInitialGuess (dmmg, ?);
  CHKERRQ(ierr);
#endif
  
  
  
#endif
  
  
  
  /* Calculate initial residual */
  
  if (verbose) {
    ierr = TATPETSc_function (snes, x, f, &user);
    CHKERRQ(ierr);
    double residual;
    ierr = VecNorm (f, NORM_2, &residual);
    CHKERRQ(ierr);
    CCTK_VInfo (CCTK_THORNSTRING, "Initial residual L2-norm: %g", residual);
  }
  
  
  
  /* Solve the nonlinear system */
  
  long clk_tck = sysconf(_SC_CLK_TCK);
  assert (clk_tck>0);
  struct tms tms_buffer;
  clock_t ticks0 = times(&tms_buffer);
  assert (ticks0!=(clock_t)-1);
  clock_t time0 = tms_buffer.tms_utime + tms_buffer.tms_stime;
  
#if 1
  
  if (verbose) CCTK_INFO ("SNESSolve");
  
  ierr = SNESSolve (snes, PETSC_NULL, x);
  CHKERRQ(ierr);

#else

  if (verbose) CCTK_INFO ("DMMGSolve");
  ierr = DMMGSolve (dmmg);
  CHKERRQ(ierr);

#endif
  
  PetscInt iters;               /* iterations needed */
  ierr = SNESGetIterationNumber (snes, &iters);
  CHKERRQ(ierr);
  PetscInt liniters;            /* linear iterations needed */
  ierr = SNESGetLinearSolveIterations (snes, &liniters);
  CHKERRQ(ierr);
  
  clock_t ticks1 = times(&tms_buffer);
  assert (ticks1!=(clock_t)-1);
  clock_t time1 = tms_buffer.tms_utime + tms_buffer.tms_stime;
  
  if (verbose) {
    CCTK_VInfo (CCTK_THORNSTRING, "Number of Newton iterations: %d",
                (int)iters);
    CCTK_VInfo (CCTK_THORNSTRING, "Number of linear iterations: %d",
                (int)liniters);
    CCTK_VInfo (CCTK_THORNSTRING, "Number of function calls:    %d",
		user.funcall_count);
    if (jac) {
      CCTK_VInfo (CCTK_THORNSTRING, "Number of Jacobian calls:    %d",
		  user.jaccall_count);
    }
    CCTK_VInfo (CCTK_THORNSTRING, "Wall time: %g sec, CPU time: %g sec",
		(double)(ticks1-ticks0)/(double)clk_tck,
		(double)(time1-time0)/(double)clk_tck);
  }
  
  
  
  /* Analyse result */
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason (snes, &reason);
  CHKERRQ(ierr);
  
  const char *msg;
  switch (reason) {
  case SNES_CONVERGED_FNORM_ABS:      msg = "SNES_CONVERGED_FNORM_ABS: ||F|| < atol"; break;
  case SNES_CONVERGED_FNORM_RELATIVE: msg = "SNES_CONVERGED_FNORM_RELATIVE: ||F|| < rtol*||F_initial||"; break;
  case SNES_CONVERGED_SNORM_RELATIVE: msg = "SNES_CONVERGED_SNORM_RELATIVE: Newton computed step size small; || delta x || < stol || x ||"; break;
  case SNES_CONVERGED_ITS:            msg = "SNES_CONVERGED_ITS: maximum iterations reached"; break;
  case SNES_CONVERGED_TR_DELTA:       msg = "SNES_CONVERGED_TR_DELTA: diverged"; break;
  case SNES_DIVERGED_FUNCTION_DOMAIN: msg = "SNES_DIVERGED_FUNCTION_DOMAIN: the new x location passed the function is not in the domain of F"; break;
  case SNES_DIVERGED_FUNCTION_COUNT:  msg = "SNES_DIVERGED_FUNCTION_COUNT"; break;
  case SNES_DIVERGED_LINEAR_SOLVE:    msg = "SNES_DIVERGED_LINEAR_SOLVE: the linear solve failed"; break;
  case SNES_DIVERGED_FNORM_NAN:       msg = "SNES_DIVERGED_FNORM_NAN"; break;
  case SNES_DIVERGED_MAX_IT:          msg = "SNES_DIVERGED_MAX_IT"; break;
  case SNES_DIVERGED_LINE_SEARCH:     msg = "SNES_DIVERGED_LINE_SEARCH: the line search failed"; break;
  case SNES_DIVERGED_INNER:           msg = "SNES_DIVERGED_INNER: inner solve failed"; break;
  case SNES_DIVERGED_LOCAL_MIN:       msg = "SNES_DIVERGED_LOCAL_MIN: || J^T b || is small, implies converged to local minimum of F()"; break;
  case SNES_CONVERGED_ITERATING:      msg = "ITERATING"; break;
  default:                            msg = "(unknown reason)";
  }
  
  if ((int)reason<=0) {
    CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING, "SNESsolve diverged for reason %d: %s", (int)reason, msg);
  } else {
    if (verbose) {
      CCTK_VInfo (CCTK_THORNSTRING, "SNESsolve converged for reason %d: %s", (int)reason, msg);
    }
  }
  
  
  
  /* Calculate final residual */
  
  ierr = TATPETSc_function (snes, x, f, &user);
  CHKERRQ(ierr);
  if (verbose) {
    double residual;		/* L2-norm of residual */
    ierr = VecNorm (f, NORM_2, &residual);
    CHKERRQ(ierr);
    CCTK_VInfo (CCTK_THORNSTRING, "Final residual L2-norm: %g", residual);
  }
  
  
  
  /* Get solution and residual */
  
  ierr = TATPETSc_copy (x, &user, TATcopyout, TATcopyvars);
  CHKERRQ(ierr);
  ierr = TATPETSc_copy (f, &user, TATcopyout, TATcopyvals);
  CHKERRQ(ierr);
  
  
  
  /* Free work space */
  
#if 1
  
  if (veryverbose) CCTK_INFO ("VecDestroy");
  ierr = VecDestroy (&x);
  CHKERRQ(ierr);
  ierr = VecDestroy (&f);
  CHKERRQ(ierr);
  
  if (veryverbose) CCTK_INFO ("MatDestroy");
  ierr = MatDestroy (&J);
  CHKERRQ(ierr);
  
  if (!jac) {
    if (veryverbose) CCTK_INFO ("MatFDColoringDestroy");
    ierr = MatFDColoringDestroy (&matfdcoloring);
    CHKERRQ(ierr);
  }
  
  if (veryverbose) CCTK_INFO ("SNESDestroy");
  ierr = SNESDestroy (&snes);
  CHKERRQ(ierr);
  
#else
  
  if (veryverbose) CCTK_INFO ("DMMGDestroy");
  ierr = DMMGDestroy (&dmmg);
  CHKERRQ(ierr);
  
#endif
  
  if (veryverbose) CCTK_INFO ("DMDestroy");
  ierr = DMDestroy (&da);
  CHKERRQ(ierr);
  
  user.nvars = 0;
  free (user.var);
  user.var = 0;
  free (user.val);
  user.val = 0;
  free (user.jac0);
  user.jac0 = 0;
  
  free (all_lbnd);
  free (all_lsh);
  for (int d=0; d<user.dyndata.dim; ++d) free (npoints_proc[d]);
  for (int d=0; d<user.dyndata.dim; ++d) free (lx[d]);
  
  /* Reinstall error handler */
  ierr = PetscPushErrorHandler (TATPETSc_error_handler, 0);
  CHKERRQ(ierr);
  
  if (veryverbose) CCTK_INFO ("*** TATPETSc_solve done.");
  
  return reason<0;
}
