/* (C) 2001-04-18 Erik Schnetter <schnetter@uni-tuebingen.de> */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <petsc.h>
#include <petscversion.h>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "TATPETSc.h"
#include <TATelliptic.h>



/* A new error handler that does nothing */
PetscErrorCode TATPETSc_error_handler (MPI_Comm comm, int line,
                                       const char *func, const char *file,
                                       PetscErrorCode n, PetscErrorType p,
                                       const char *mess, void *ctx)
{
  return n;
}



int TATPETSc_initialize (void)
{
  DECLARE_CCTK_PARAMETERS;
  
  /* world communicator */
  MPI_Comm comm;
  
  static int argc;
  static char *args;
  static char **argv;
  
  PetscInt ierr;
  
  
  
  if (veryverbose) CCTK_INFO ("PetscInitialize");
  
  
  
  /* Create PETSc command line options */
  
  args = strdup(options);
  assert (args);
  
  int len = strlen(args);
  argc = 2;
  for (int i=0; i<len; ++i) {
    if (args[i]==' ' || args[i]=='\n') {
      args[i] = '\0';
      ++argc;
    }
  }
  
  argv = malloc(sizeof(*argv) * (argc+1));
  assert (argv);
  argv[0] = "Cactus";
  argv[1] = args;
  argc = 2;
  for (int i=1; i<len; ++i) {
    if (args[i-1]=='\0') {
      argv[argc] = &args[i];
      ++argc;
    }
  }
  argv[argc] = 0;
  
  if (veryverbose) {
    CCTK_INFO ("PETSc command line arguments:");
    for (int i=0; i<argc; ++i) {
      CCTK_VInfo (CCTK_THORNSTRING, "   %d: %s", i, argv[i]);
    }
    CCTK_INFO ("End of PETSc command line arguments.");
  }
  
  
  
  /* Initialise PETSc */
  if (CCTK_IsFunctionAliased ("GetMPICommWorld")) {
    comm = * (MPI_Comm const *) GetMPICommWorld (NULL);
  } else {
    comm = MPI_COMM_WORLD;
  }
  PETSC_COMM_WORLD = comm;
  
  ierr = PetscInitialize (&argc, &argv, PETSC_NULL, PETSC_NULL);
  CHKERRQ(ierr);
  
  
  
  /* Install a new error handler */
  ierr = PetscPushErrorHandler (TATPETSc_error_handler, 0);
  CHKERRQ(ierr);
  
  
  
  /* Register the solver */
  ierr = TATelliptic_RegisterSolver (TATPETSc_solve, "TATPETSc");
  assert (!ierr);



  return 0;
}
