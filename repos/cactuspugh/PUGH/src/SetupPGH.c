 /*@@
   @file      SetupPGH.c
   @date      Fri Feb 21 10:18:17 1997
   @author    Paul Walker
   @desc
              Initializes the Pugh Grid Hierachy.
   @enddesc
   @version   $Id$
 @@*/

/* #define DEBUG_PUGH 1 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "cctk.h"
#include "cctk_Parameters.h"
#include "pugh.h"
#include "pughi.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusPUGH_PUGH_SetupPGH_c);


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
int PUGH_SetTopology (int dim, const int topology[]);


/********************************************************************
 *********************     Internal Routines   **********************
 ********************************************************************/
static int Setup_nProcs (pGH *pughGH, int dim);


/********************************************************************
 *********************     Static Variables    **********************
 ********************************************************************/
/* static integers which contain partition information if set */
static int topology_1d_x = 0;
static int topology_2d_x = 0;
static int topology_2d_y = 0;
static int topology_3d_x = 0;
static int topology_3d_y = 0;
static int topology_3d_z = 0;


 /*@@
   @routine    PUGH_SetupPGH
   @date       Fri Feb 21 10:21:36 1997
   @author     Gabrielle Allen, Tom Goodale, Paul Walker
   @desc
               Sets up the PGH distribution based on processor number and
               problem size.
               This includes setting an optimal domain decomposition for each
               dimension, and setting ghosts and overlaps.
   @enddesc

   @var        callerid
   @vdesc      back reference to the calling cGH structure
   @vtype      void *
   @vio        in
   @endvar
   @var        dim
   @vdesc      number of dimensions to set up
   @vtype      int
   @vio        in
   @endvar
   @var        nsize
   @vdesc      number of grid points in every dimension
   @vtype      int *
   @vio        in
   @endvar
   @var        nghostzones
   @vdesc      number of ghostzones in every dimension
   @vtype      int *
   @vio        in
   @endvar
   @var        perme
   @vdesc      periodic boundary flags for every dimension
   @vtype      int *
   @vio        in
   @endvar

   @returntype pGH *
   @returndesc
               the pointer to the new pGH structure
   @endreturndesc
@@*/
pGH *PUGH_SetupPGH (void *callerid,
                    int dim,
                    int *nsize,
                    int *nghostzones,
                    int *perme)
{
  DECLARE_CCTK_PARAMETERS
  pGH *pughGH;
  int i, idim;
  int *nprocs;
  int type;
  const void *allow_mixeddim_gfs_ptr;
  CCTK_INT allow_mixeddim_gfs;


  /* Allocate for myself */
  pughGH = malloc (sizeof (pGH));

  /* Set an identifier for my parent */
  pughGH->callerid = callerid;

  pughGH->GFExtras     = malloc (dim * sizeof (pGExtras *));
  pughGH->Connectivity = malloc (dim*sizeof(pConnectivity*));

  /* Set the total number of processors */
  Setup_nProcs (pughGH, dim);

  /* Determine whether mixed-dimensional grid functions need to be
     supported */
  allow_mixeddim_gfs_ptr =
    CCTK_ParameterGet ("allow_mixeddim_gfs", "Cactus", &type);
  if (! allow_mixeddim_gfs_ptr)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "internal error -- flesh parameter Cactus::allow_mixeddim_gfs does not exist");
  }
  if (type != PARAMETER_BOOLEAN)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "internal error -- flesh parameter Cactus::allow_mixeddim_gfs has wrong type");
  }
  allow_mixeddim_gfs = *(const CCTK_INT*)allow_mixeddim_gfs_ptr;

  /* Set up connectivity and extras for each dimension */
  for (idim = 1; idim <= dim; idim++)
  {
    /* Do not set up topologies for lower-dimensional grid functions
       if these are not allowed.  This prevents problems if such
       topologies do not exist, e.g. on small grids.  */
    if (allow_mixeddim_gfs || idim == dim)
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Setting up a topology for %d dimensions", idim);
      nprocs = malloc (idim * sizeof (int));

      PUGH_SetupDefaultTopology (idim, nprocs);

      /* Check that there are enough grid points in this dimension
       * to make parallelising it worthwhile
       */
      for (i = 0; i < idim; i++)
      {
        if ((! nprocs[i]) && (abs (nsize[i]) <= 2 * nghostzones[i] + 1))
        {
          nprocs[i] = 1;
        }
      }

      pughGH->Connectivity[idim-1] =
        PUGH_SetupConnectivity (idim, pughGH->nprocs,
                                nsize, nghostzones,
                                nprocs, perme);
      free(nprocs);

      pughGH->GFExtras[idim-1] =
        PUGH_SetupPGExtras (1, idim, perme,
                            nsize, nghostzones,
                            pughGH->nprocs,
                            pughGH->Connectivity[idim-1]->nprocs,
                            pughGH->Connectivity[idim-1]->neighbours,
                            pughGH->myproc);
    }
    else
    {
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Not setting up a topology for %d dimensions", idim);
      pughGH->Connectivity[idim-1] = NULL;
      pughGH->GFExtras[idim-1] = NULL;
    }
  }

  /* create the timer for communication time */
  pughGH->comm_time = timer_output ? CCTK_TimerCreateI () : -1;

  pughGH->active = 1;
  pughGH->commmodel = PUGH_ALLOCATEDBUFFERS;
  pughGH->identity_string = NULL;
  pughGH->level = 0;
  pughGH->mglevel = 0;
  pughGH->convlevel = 0;
  pughGH->forceSync = 0;
  pughGH->nvariables = 0;
  pughGH->narrays = 0;
  pughGH->activetimelevels = NULL;
  pughGH->maxtimelevels = NULL;
  pughGH->variables = NULL;

  return (pughGH);
}


 /*@@
   @routine    PUGH_Terminate
   @date       Tue Apr 18 15:21:42 2000
   @author     Tom Goodale
   @desc
               PUGH's termination routine to destroy a pGH structure.
   @enddesc

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
void PUGH_Terminate (cGH *GH)
{
  pGH *pughGH = PUGH_pGH (GH);


  PUGH_DestroyPGH (&pughGH);
}


 /*@@
   @routine    PUGH_DestroyPGH
   @date       Thu Aug 21 11:38:10 1997
   @author     Paul Walker
   @desc
               Destroys a GH object.
   @enddesc

   @var        GHin
   @vdesc      address of PUGH GH extensions object to be destroyed
   @vtype      pGH **
   @vio        in
   @endvar
@@*/
void PUGH_DestroyPGH (pGH **GHin)
{
  pGH    *GH;
  pGA    *GA;
  cGroup pgroup;
  int i, group, numgroups, this_var, variable;


  GH = *GHin;

#ifdef CCTK_MPI
  CACTUS_MPI_ERROR (MPI_Comm_free (&GH->PUGH_COMM_WORLD));

  CACTUS_MPI_ERROR (MPI_Type_free (&GH->PUGH_mpi_complex));
#ifdef HAVE_CCTK_REAL4
  CACTUS_MPI_ERROR (MPI_Type_free (&GH->PUGH_mpi_complex8));
#endif
#ifdef HAVE_CCTK_REAL8
  CACTUS_MPI_ERROR (MPI_Type_free (&GH->PUGH_mpi_complex16));
#endif
#ifdef HAVE_CCTK_REAL16
  CACTUS_MPI_ERROR (MPI_Type_free (&GH->PUGH_mpi_complex32));
#endif
#endif

  /* Great. Now go about the work of destroying me. */
  variable = 0;
  numgroups = CCTK_NumGroups ();
  for (group = 0; group < numgroups; group++)
  {
#ifdef DEBUG_PUGH
    printf ("Calling Destroying Group %s\n", CCTK_GroupName (group));
    fflush (stdout);
#endif

    CCTK_GroupData (group, &pgroup);
    if (! pgroup.numvars)
    {
      continue;
    }

    GA = GH->variables[variable][0];

    /* Destroy group comm buffers */
    if (GA->groupcomm)
    {
      if (GA->groupcomm->commflag != PUGH_NOCOMM)
      {
        PUGH_DisableGArrayGroupComm (GH, variable, GA->groupcomm);
      }
      PUGH_DestroyComm (&GA->groupcomm);
    }

    /* Destroy the group's connectivity and extras structure
       for CCTK_ARRAY groups.
       Remember that the connectivity and extras for CCTK_GF types
       are shared between all such groups and are destroyed later. */
    if (pgroup.dim == 0 || GA->connectivity != GH->Connectivity[pgroup.dim-1])
    {
      PUGH_DestroyConnectivity (&GA->connectivity);
    }
    if (pgroup.dim == 0 || GA->extras != GH->GFExtras[pgroup.dim-1])
    {
      PUGH_DestroyPGExtras (&GA->extras);
    }

    for (this_var = 0; this_var < pgroup.numvars; this_var++, variable++)
    {
      for (i = 0 ; i < pgroup.numtimelevels; i++)
      {
        PUGH_DestroyGArray (&(((pGA ***)GH->variables)[variable][i]));
      }
      free (GH->variables[variable]);
    }
  }

  for (i = 0; i < GH->dim; i++)
  {
    PUGH_DestroyConnectivity (GH->Connectivity + i);
    PUGH_DestroyPGExtras (GH->GFExtras + i);
  }

  free (GH->identity_string);
  free (GH->Connectivity);
  free (GH->GFExtras);
  free (GH->activetimelevels);
  free (GH->maxtimelevels);
  free (GH->variables);
  free (GH);
  *GHin = NULL;
}


 /*@@
   @routine    PUGH_SetTopology
   @date       Wed 02 Apr 2003
   @author     Fokke Dijkstra
   @desc
               Set/reset the static topology integer variables.

               Note that using such static variables here instead of the
               processor_topology_[123]d_[xyz] parameters is just a hack because
               we don't want to make those steerable (this could make
               people think they can change the processor topology at
               runtime).
   @enddesc

   @var        dim
   @vdesc      number of dimensions for the processor topology (1,2 or 3)
   @vtype      int
   @vio        in
   @endvar
   @var        topology
   @vdesc      list with number of processors in the x,y and z direction
   @vtype      int[]
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0: OK
               -1: invalid dimension
   @endreturndesc
@@*/
int PUGH_SetTopology (int dim, const int topology[])
{
  int retval;


  retval = 0;
  if (dim == 1)
  {
    topology_1d_x = topology[0];
  }
  else if (dim == 2)
  {
    topology_2d_x = topology[0];
    topology_2d_y = topology[1];
  }
  else if (dim == 3)
  {
    topology_3d_x = topology[0];
    topology_3d_y = topology[1];
    topology_3d_z = topology[2];
  }
  else
  {
    CCTK_WARN (1, "Only 1D, 2D, and 3D supported");
    retval = -1;
  }

  return (retval);
}


 /*@@
   @routine    PUGH_SetupDefaultTopology
   @date       Sun 23 Jan 1999
   @author     Gabrielle Allen
   @desc
               Set the processor decomposition from PUGH parameters.
   @enddesc

   @var        dim
   @vdesc      topology dimension
   @vtype      int
   @vio        in
   @endvar
   @var        nprocs
   @vdesc      number of processors in each dimension
   @vtype      int []
   @vio        out
   @endvar

   @returntype int
   @returndesc
               PUGH_SUCCESS for success, or PUGH_ERROR for invalid dimension
   @endreturndesc
@@*/
int PUGH_SetupDefaultTopology (int dim, int nprocs[])
{
  DECLARE_CCTK_PARAMETERS
  int retval;


  retval = PUGH_SUCCESS;

  switch (dim)
  {
    case 1:
      nprocs[0] = topology_1d_x ? topology_1d_x : processor_topology_1d_x;
      break;
    case 2:
      nprocs[0] = topology_2d_x ? topology_2d_x : processor_topology_2d_x;
      nprocs[1] = topology_2d_y ? topology_2d_y : processor_topology_2d_y;
      break;
    case 3:
      nprocs[0] = topology_3d_x ? topology_3d_x : processor_topology_3d_x;
      nprocs[1] = topology_3d_y ? topology_3d_y : processor_topology_3d_y;
      nprocs[2] = topology_3d_z ? topology_3d_z : processor_topology_3d_z;
      break;
    default:
      memset (nprocs, 0, dim * sizeof (int));
      retval = PUGH_ERROR;
  }

  return (retval);
}


 /*@@
   @routine    PUGH_ParallelInit
   @date       Tue Apr 18 15:21:42 2000
   @author     Tom Goodale
   @desc
               PUGH overloadable routine for CCTK_ParallelInit().
   @enddesc
   @calls

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        unused
   @endvar

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
int PUGH_ParallelInit (cGH *GH)
{
  /* avoid compiler warning about unused parameter */
  GH = GH;

  return (0);
}


 /*@@
   @routine    PUGH_Exit
   @date       Saturday July 15 2000
   @author     Gabrielle Allen
   @desc
               PUGH overloadable routine for CCTK_Exit().
   @enddesc
   @calls      MPI_Finalize
               exit

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar
   @var        retval
   @vdesc      return code to exit with
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               This function should never return.
               But if it does it will return the given return code.
   @endreturndesc
@@*/
int PUGH_Exit (cGH *GH, int retval)
{
  /* avoid compiler warning about unused parameter */
  GH = GH;

#ifdef CCTK_MPI
  CACTUS_MPI_ERROR (MPI_Finalize ());
#endif
  exit (retval);
  return (retval);
}


 /*@@
   @routine    PUGH_Abort
   @date       Saturday July 15 2000
   @author     Gabrielle Allen
   @desc
               PUGH overloadable routine for CCTK_Abort().
   @enddesc
   @calls      MPI_Abort
               exit

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      cGH *
   @vio        in
   @endvar
   @var        retval
   @vdesc      return code to abort with
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               This function should never return.
               But if it does it will return the return code to abort with.
   @endreturndesc
@@*/
int PUGH_Abort (cGH *GH, int retval)
{
#ifdef CCTK_MPI
  /* flush stdout and stderr before calling MPI_Abort() because
     some MPI implementations simply kill other MPI processes */
  fflush (stdout);
  fflush (stderr);
  CACTUS_MPI_ERROR (MPI_Abort (GH ? PUGH_pGH (GH)->PUGH_COMM_WORLD :
                                    MPI_COMM_WORLD, retval));
#else
  /* avoid compiler warning about unused parameter */
  GH = GH;

  /* FIXME */
  /*abort();*/
#endif
  exit (retval);
  return (retval);
}


 /*@@
   @routine    PUGH_MyProc
   @date       Tue Jan 23 1999
   @author     Gabrielle Allen
   @desc
               PUGH overloadable routine for CCTK_MyProc().
   @enddesc
   @calls      MPI_Comm_rank

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the processor number of the caller
   @endreturndesc
@@*/
int PUGH_MyProc (const cGH *GH)
{
  int myproc;


#ifdef CCTK_MPI
  CACTUS_MPI_ERROR (MPI_Comm_rank (GH ? PUGH_pGH (GH)->PUGH_COMM_WORLD :
                                        MPI_COMM_WORLD, &myproc));
#else
  GH = GH;
  myproc = 0;
#endif

  return (myproc);
}


 /*@@
   @routine    PUGH_nProcs
   @date       Tue Jan 23 1999
   @author     Gabrielle Allen
   @desc
               PUGH overloadable routine for CCTK_nProcs().
   @enddesc
   @calls      MPI_Comm_size

   @var        GH
   @vdesc      Pointer to CCTK grid hierarchy
   @vtype      const cGH *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the total number of processors
   @endreturndesc
@@*/
int PUGH_nProcs (const cGH *GH)
{
  int nprocs;


#ifdef CCTK_MPI
  CACTUS_MPI_ERROR (MPI_Comm_size (GH ? PUGH_pGH (GH)->PUGH_COMM_WORLD :
                                        MPI_COMM_WORLD, &nprocs));
#else
  GH = GH;
  nprocs = 1;
#endif

  return (nprocs);
}


#if 0
void pGH_DumpInfo(pGH *GH)
{
  int j,k;
  printf("INFO: \n");
  printf("   myproc: %d/%d \n",GH->myproc,GH->nprocs);
  printf("   ownership: \n");
  for (j=0;j<GH->dim;j++)
  {
    printf("     GH->ownership[0/1][%d]: %d %d \n",
           j,GH->GFExtras[GH->dim-1]->ownership[0][j],GH->GFExtras[GH->dim-1]->ownership[1][j]);
  }

  for (j=0;j<GH->dim;j++)
  {
    for (k=0;k<2*GH->dim;k++)
    {
      printf("     GH->GFExtras[GH->dim-1]->ghosts[0/1][%d][%d]: %d %d \n",
             k,j,GH->GFExtras[GH->dim-1]->ghosts[0][k][j],GH->GFExtras[GH->dim-1]->ghosts[1][k][j]);
    }
  }

  for (j=0;j<GH->dim;j++)
  {
    for (k=0;k<2*GH->dim;k++)
    {
      printf("     GH->GFExtras[GH->dim-1]->overlap[0/1][%d][%d]: %d %d \n",
             k,j,GH->GFExtras[GH->dim-1]->overlap[0][j][k],GH->GFExtras[GH->dim-1]->overlap[1][k][j]);
    }
  }
}
#endif


/********************************************************************
 *********************     Internal Routines   **********************
 ********************************************************************/
 /*@@
   @routine    Setup_nProcs
   @date       Tue Apr 18 15:21:42 2000
   @author     Tom Goodale
   @desc
               Setup PUGH_COMM_WORLD communicator and make sure that CCTK_INT
               and CCTK_REAL sizes are the same on all processors.
               Define MPI datatypes for CCTK_REAL and CCTK_COMPLEX variables.
   @enddesc

   @var        pughGH
   @vdesc      Pointer to PUGH grid hierarchy
   @vtype      pGH *
   @vio        in
   @endvar
   @var        dim
   @vdesc      dimension of processor topology to set up
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success
   @endreturndesc
@@*/
static int Setup_nProcs (pGH *pughGH, int dim)
{
#ifdef CCTK_MPI
  CCTK_REAL4 sum_sizes [2], compiled_sizes [2];


  /* Set up my communicator. This would allow us to run pugh on
     a subset of processors at a later date if, for instance, we
     were using panda or what not.
   */
  CACTUS_MPI_ERROR (MPI_Comm_dup (MPI_COMM_WORLD, &pughGH->PUGH_COMM_WORLD));

  CACTUS_MPI_ERROR (MPI_Comm_size (pughGH->PUGH_COMM_WORLD, &pughGH->nprocs));
  CACTUS_MPI_ERROR (MPI_Comm_rank (pughGH->PUGH_COMM_WORLD, &pughGH->myproc));

  /* check that all executables uses the same integer and fp precision
     within a metacomputing environment
     NOTE: We cannot use CCTK_INT4s in MPI_Allreduce() since
           (although they have the same size) they might be correspond to
           different MPI datatypes.
           CCTK_REAL4 should always refer to MPI_FLOAT. */
  compiled_sizes [0] = sizeof (CCTK_INT);
  compiled_sizes [1] = sizeof (CCTK_REAL);
  CACTUS_MPI_ERROR (MPI_Allreduce (compiled_sizes, sum_sizes, 2, PUGH_MPI_REAL4,
                                   MPI_SUM, pughGH->PUGH_COMM_WORLD));
  if (compiled_sizes [0] * pughGH->nprocs != sum_sizes [0])
  {
    CCTK_WARN (0, "Cannot run executables with different precision for "
                  "CCTK_INTs within a metacomputing environment !\n"
                  "Please configure with unique CCTK_INTEGER_PRECISION !");
  }

  if (compiled_sizes [1] * pughGH->nprocs != sum_sizes [1])
  {
    CCTK_WARN (0, "Cannot run executables with different precision for "
                  "CCTK_REALs within a metacomputing environment !\n"
                  "Please configure with unique CCTK_REAL_PRECISION !");
  }

  /* define the complex datatype as a concatanation of 2 PUGH_MPI_REALs */
  CACTUS_MPI_ERROR (MPI_Type_contiguous (2, PUGH_MPI_REAL,
                                         &pughGH->PUGH_mpi_complex));
  CACTUS_MPI_ERROR (MPI_Type_commit (&pughGH->PUGH_mpi_complex));

  /* dito for fixed-precision reals */
#ifdef HAVE_CCTK_REAL4
  CACTUS_MPI_ERROR (MPI_Type_contiguous (2, PUGH_MPI_REAL4,
                                         &pughGH->PUGH_mpi_complex8));
  CACTUS_MPI_ERROR (MPI_Type_commit (&pughGH->PUGH_mpi_complex8));
#endif
#ifdef HAVE_CCTK_REAL8
  CACTUS_MPI_ERROR (MPI_Type_contiguous (2, PUGH_MPI_REAL8,
                                         &pughGH->PUGH_mpi_complex16));
  CACTUS_MPI_ERROR (MPI_Type_commit (&pughGH->PUGH_mpi_complex16));
#endif
#ifdef HAVE_CCTK_REAL16
  CACTUS_MPI_ERROR (MPI_Type_contiguous (2, PUGH_MPI_REAL16,
                                         &pughGH->PUGH_mpi_complex32));
  CACTUS_MPI_ERROR (MPI_Type_commit (&pughGH->PUGH_mpi_complex32));
#endif

#else
  pughGH->nprocs = 1;
  pughGH->myproc = 0;
#endif

  pughGH->dim = dim;

  return (0);
}
