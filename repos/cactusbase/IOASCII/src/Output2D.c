/*@@
  @file      Output2D.c
  @date      Thu May 11 2000
  @author    Thomas Radke
  @desc
             Functions to deal with 2D ASCII output of variables
  @enddesc
  @version   $Id$
@@*/

#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "util_String.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"
#include "ioASCIIGH.h"

/* the rcs ID and its dummy function to use it */
static const char *rcsid = "$Id$";
CCTK_FILEVERSION(CactusBase_IOASCII_Output2D_c)

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static int CheckOutputVar(int vindex);
static void SetOutputFlag(int vindex, const char *optstring, void *arg);

/*@@
   @routine    IOASCII_Output2DGH
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Loops over all variables and outputs them if necessary
   @enddesc
   @calls      IOASCII_TimeFor2D
               IOASCII_Write2D

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               the number of variables which were output at this iteration
               (or 0 if it wasn't time to output yet)
   @endreturndesc
@@*/
int IOASCII_Output2DGH(const cGH *GH) {
  int vindex, retval;
  const asciiioGH *myGH;

  retval = 0;
  myGH = (const asciiioGH *)CCTK_GHExtension(GH, "IOASCII");

  /* loop over all variables */
  for (vindex = CCTK_NumVars() - 1; vindex >= 0; vindex--) {
    if (IOASCII_TimeFor2D(GH, vindex) &&
        IOASCII_Write2D(GH, vindex, CCTK_VarName(vindex)) == 0) {
      /* register variable as having 2D output this iteration */
      myGH->out2D_last[vindex] = GH->cctk_iteration;
      retval++;
    }
  }

  return (retval);
}

/*@@
   @routine    IOASCII_Output2DVarAs
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Unconditional output of a variable
               using the IOASCII 2D output method
   @enddesc
   @calls      IOASCII_Write2D

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        fullname
   @vdesc      complete name of variable to output
   @vtype      const char *
   @vio        in
   @endvar
   @var        alias
   @vdesc      alias name of variable to output
               (used to generate output filename)
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOASCII_Write2D, or<BR>
               -1 if variable cannot be output by IOASCII_2D
   @endreturndesc
@@*/
int IOASCII_Output2DVarAs(const cGH *GH, const char *fullname,
                          const char *alias) {
  int vindex, retval;

  retval = -1;
  vindex = CCTK_VarIndex(fullname);
  if (vindex < 0) {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "There is no such variable %s.  Ignoring IOASCII 2D output.",
               fullname);
  } else if (CheckOutputVar(vindex) == 0) {
    retval = IOASCII_Write2D(GH, vindex, alias);
  }

  return (retval);
}

/*@@
  @routine    IOASCII_TimeFor2D
  @date       Sat March 6 1999
  @author     Gabrielle Allen
  @desc
              Decides if it is time to output a variable
              using the IOASCII 2D output method
  @enddesc

  @var        GH
  @vdesc      pointer to CCTK GH
  @vtype      const cGH *
  @vio        in
  @endvar
  @var        vindex
  @vdesc      index of variable
  @vtype      int
  @vio        in
  @endvar

  @returntype int
  @returndesc
              1 if output should take place at this iteration, or<BR>
              0 if not
  @endreturndesc
@@*/
int IOASCII_TimeFor2D(const cGH *GH, int vindex) {
  int retval;
  char *fullname;
  asciiioGH *myGH;

  myGH = (asciiioGH *)CCTK_GHExtension(GH, "IOASCII");
  IOASCII_CheckSteerableParameters2D(myGH);

  /* check if this variable should be output */
  retval = myGH->out2D_every[vindex] > 0 &&
           GH->cctk_iteration % myGH->out2D_every[vindex] == 0;
  if (retval) {
    /* check if variable wasn't already output this iteration */
    if (myGH->out2D_last[vindex] == GH->cctk_iteration) {
      fullname = CCTK_FullName(vindex);
      CCTK_VWarn(5, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Already done IOASCII 2D output for '%s' in current "
                 "iteration (probably via triggers)",
                 fullname);
      free(fullname);
      retval = 0;
    }
  }

  return (retval);
}

/*@@
   @routine    IOASCII_TriggerOutput2D
   @date       Sat March 6 1999
   @author     Gabrielle Allen
   @desc
               Triggers the output of a variable
               using the IOASCII 2D output method
   @enddesc
   @calls      IOASCII_Write2D

   @var        GH
   @vdesc      pointer to CCTK GH
   @vtype      const cGH *
   @vio        in
   @endvar
   @var        vindex
   @vdesc      index of variable to output
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               return code of @seeroutine IOASCII_Write2D
   @endreturndesc
@@*/
int IOASCII_TriggerOutput2D(const cGH *GH, int vindex) {
  int retval;
  asciiioGH *myGH;

  /* do the 2D output */
  retval = IOASCII_Write2D(GH, vindex, CCTK_VarName(vindex));
  if (retval == 0) {
    /* register variable as having 2D output this iteration */
    myGH = (asciiioGH *)CCTK_GHExtension(GH, "IOASCII");
    myGH->out2D_last[vindex] = GH->cctk_iteration;
  }

  return (retval);
}

/*@@
   @routine    IOASCII_CheckSteerableParameters2D
   @date       Tue 8 June 2004
   @author     Thomas Radke
   @desc
               Check if steerable IOASCII 2D parameters have changed
   @enddesc
   @calls      CCTK_TraverseString

   @var        myGH
   @vdesc      Pointer to IOASCII GH
   @vtype      asciiioGH *
   @vio        inout
   @endvar
@@*/
void IOASCII_CheckSteerableParameters2D(asciiioGH *myGH) {
  int i, num_vars;
  char *fullname, *msg;
  DECLARE_CCTK_PARAMETERS

  /* how often to output */
  i = myGH->out2D_every_default;
  myGH->out2D_every_default = out2D_every >= 0 ? out2D_every : out_every;

  /* report if frequency changed */
  if (myGH->out2D_every_default != i && !CCTK_Equals(verbose, "none")) {
    if (myGH->out2D_every_default > 0) {
      CCTK_VInfo(CCTK_THORNSTRING, "Periodic 2D output every %d "
                                   "iterations",
                 myGH->out2D_every_default);
    } else {
      CCTK_INFO("Periodic 2D output turned off");
    }
  }

  /* re-parse the 'out2D_vars' parameter if it was changed */
  if (strcmp(out2D_vars, myGH->out2D_vars) || myGH->out2D_every_default != i) {
    num_vars = CCTK_NumVars();
    memset(myGH->out2D_every, 0, num_vars * sizeof(CCTK_INT));
    if (CCTK_TraverseString(out2D_vars, SetOutputFlag, myGH,
                            CCTK_GROUP_OR_VAR) < 0) {
      CCTK_WARN(myGH->stop_on_parse_errors ? 0 : 1,
                "error while parsing parameter 'IOASCII::out2D_vars'");
    }

    if (myGH->out2D_every_default == i || !CCTK_Equals(verbose, "none")) {
      msg = NULL;
      for (i = 0; i < num_vars; i++) {
        if (myGH->out2D_every[i]) {
          fullname = CCTK_FullName(i);
          if (!msg) {
            Util_asprintf(&msg, "Periodic 2D output requested for "
                                "'%s'",
                          fullname);
          } else {
            char *tmp = msg;
            Util_asprintf(&msg, "%s, '%s'", msg, fullname);
            free(tmp);
          }
          free(fullname);
        }
      }
      if (msg) {
        CCTK_INFO(msg);
        free(msg);
      }
    }

    /* save the last setting of 'out2D_vars' parameter */
    free(myGH->out2D_vars);
    myGH->out2D_vars = strdup(out2D_vars);
  }
}

/**************************** local functions ******************************/
/* check if this variable can be output (static conditions) */
static int CheckOutputVar(int vindex) {
  int groupindex;
  cGroup groupinfo;
  char *fullname;
  const char *errormsg;

  /* get the variable group information */
  groupindex = CCTK_GroupIndexFromVarI(vindex);
  CCTK_GroupData(groupindex, &groupinfo);

  errormsg = NULL;
  if (groupinfo.dim < 2 || groupinfo.dim > 3) {
    errormsg = "dim != [2,3]";
  } else if (groupinfo.grouptype != CCTK_GF &&
             groupinfo.grouptype != CCTK_ARRAY) {
    errormsg = "not a grid function or array";
  }

  if (errormsg) {
    fullname = CCTK_FullName(vindex);
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "No IOASCII 2D output for '%s': %s", fullname, errormsg);
    free(fullname);
  }

  return (errormsg != NULL);
}

/* callback for CCTK_TraverseString() to set the output flag
   for the given variable */
static void SetOutputFlag(int vindex, const char *optstring, void *arg) {
  const asciiioGH *myGH = (const asciiioGH *)arg;

  if (CheckOutputVar(vindex) == 0) {
    myGH->out2D_every[vindex] = myGH->out2D_every_default;

    if (optstring) {
      IOUtil_ParseOutputFrequency("2D IOASCII", "IOASCII::out2D_vars",
                                  myGH->stop_on_parse_errors, vindex, optstring,
                                  &myGH->out2D_every[vindex], NULL);
    }
  }
}
