
# Step 1: Import needed core NRPy+ and Python modules
from outputC import lhrh         # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import os, sys                   # Standard Python modules for multiplatform OS-level functions
# We need the function keep_param__return_type() from this module:
import BaikalETK.BaikalETK_ETK_ccl_files_codegen as ccl

make_code_defn_list = []
def append_to_make_code_defn_list(filename):
    if filename not in make_code_defn_list:
        make_code_defn_list.append(filename)
    return filename

def driver_C_codes(Csrcdict, ThornName,
                   rhs_list,evol_gfs_list,aux_gfs_list,auxevol_gfs_list,
                   LapseCondition = "OnePlusLog", enable_stress_energy_source_terms=False):
    # First the ETK banner code, proudly showing the NRPy+ banner
    import NRPy_logo as logo
    outstr = """
#include <stdio.h>

void BaikalETK_Banner()
{
    """
    logostr = logo.print_logo(print_to_stdout=False)
    outstr += "printf(\"BaikalETK: another Einstein Toolkit thorn generated by\\n\");\n"
    for line in logostr.splitlines():
        outstr += "    printf(\""+line+"\\n\");\n"
    outstr += "}\n"

    # Finally add C code string to dictionaries (Python dictionaries are immutable)
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("Banner.c")] = outstr.replace("BaikalETK",ThornName)

    # Then the RegisterSlicing() function, needed for other ETK thorns
    outstr = """
#include "cctk.h"

#include "Slicing.h"

int BaikalETK_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("BaikalETK");
  return 0;
}"""

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("RegisterSlicing.c")] = outstr.replace("BaikalETK",ThornName)

    # Next BaikalETK_Symmetry_registration(): Register symmetries

    full_gfs_list = []
    full_gfs_list.extend(evol_gfs_list)
    full_gfs_list.extend(auxevol_gfs_list)
    full_gfs_list.extend(aux_gfs_list)

    outstr = """
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void BaikalETK_Symmetry_registration(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Stores gridfunction parity across x=0, y=0, and z=0 planes, respectively
  int sym[3];

  // Next register parities for each gridfunction based on its name
  //    (to ensure this algorithm is robust, gridfunctions with integers
  //     in their base names are forbidden in NRPy+).
"""
    outstr += ""
    for gfname in full_gfs_list:
        gfname_without_GFsuffix = gfname[:-2]
        # Do not add T4UU gridfunctions if enable_stress_energy_source_terms==False:
        if not (enable_stress_energy_source_terms == False and "T4UU" in gfname_without_GFsuffix):
            outstr += """
      // Default to scalar symmetry:
      sym[0] = 1; sym[1] = 1; sym[2] = 1;
      // Now modify sym[0], sym[1], and/or sym[2] as needed
      //    to account for gridfunction parity across
      //    x=0, y=0, and/or z=0 planes, respectively
"""
            # If gridfunction name does not end in a digit, by NRPy+ syntax, it must be a scalar
            if gfname_without_GFsuffix[len(gfname_without_GFsuffix) - 1].isdigit() == False:
                outstr += "      // (this gridfunction is a scalar -- no need to change default sym[]'s!)\n"
            elif len(gfname_without_GFsuffix) > 2:
                # Rank-1 indexed expression (e.g., vector)
                if gfname_without_GFsuffix[len(gfname_without_GFsuffix) - 2].isdigit() == False:
                    if int(gfname_without_GFsuffix[-1]) > 2:
                        print("Error: Found invalid gridfunction name: " + gfname)
                        sys.exit(1)
                    symidx = gfname_without_GFsuffix[-1]
                    outstr += "      sym[" + symidx + "] = -1;\n"
                # Rank-2 indexed expression
                elif gfname_without_GFsuffix[len(gfname_without_GFsuffix) - 2].isdigit() == True:
                    if len(gfname_without_GFsuffix) > 3 and gfname_without_GFsuffix[
                        len(gfname_without_GFsuffix) - 3].isdigit() == True:
                        print(
                            "Error: Found a Rank-3 or above gridfunction: " + gfname + ", which is at the moment unsupported.")
                        print("It should be easy to support this if desired.")
                        sys.exit(1)
                    symidx0 = gfname_without_GFsuffix[-2]
                    if "T4UU" in gfname: symidx0 = str(int(symidx0)-1) # E.g., T4UU23 is T4UUyz, corresponding to directions 1,2
                    if int(symidx0) >= 0: outstr += "      sym[" + symidx0 + "] *= -1;\n"
                    symidx1 = gfname_without_GFsuffix[-1]
                    if "T4UU" in gfname: symidx1 = str(int(symidx1)-1) # E.g., T4UU23 is T4UUyz, corresponding to directions 1,2
                    if int(symidx1) >= 0: outstr += "      sym[" + symidx1 + "] *= -1;\n"
            else:
                print(
                    "Don't know how you got this far with a gridfunction named " + gfname + ", but I'll take no more of this nonsense.")
                print("   Please follow best-practices and rename your gridfunction to be more descriptive")
                sys.exit(1)
            outstr += "      SetCartSymVN(cctkGH, sym, \"BaikalETK::" + gfname + "\");\n"
    outstr += "}\n"

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("Symmetry_registration_oldCartGrid3D.c")] = \
        outstr.replace("BaikalETK",ThornName)

    # Next set RHSs to zero
    outstr = """
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void BaikalETK_zero_rhss(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
"""
    set_rhss_to_zero = ""
    for gf in rhs_list:
        set_rhss_to_zero += gf+"[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;\n"

    outstr += lp.loop(["i2","i1","i0"],["0", "0", "0"],
                      ["cctk_lsh[2]","cctk_lsh[1]","cctk_lsh[0]"],
                      ["1","1","1"],
                      ["#pragma omp parallel for","","",],"",set_rhss_to_zero)
    outstr += "}\n"
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("zero_rhss.c")] = outstr.replace("BaikalETK",ThornName)

    # Next floor the lapse
    outstr = """
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )

void BaikalETK_floor_the_lapse(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
"""
    outstr += lp.loop(["i2", "i1", "i0"], ["0", "0", "0"],
                      ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                      ["1", "1", "1"],
                      ["#pragma omp parallel for", "", "", ], "", """
alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)] = MAX(alphaGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)], lapse_floor);
""")
    outstr += """
}
#undef MAX
#endif\n"""
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("floor_the_lapse.c")] = outstr.replace("BaikalETK", ThornName)

    # Next registration with the Method of Lines thorn
    outstr = """
//--------------------------------------------------------------------------
// Register with the Method of Lines time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// MoL documentation located in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

void BaikalETK_MoL_registration(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  // Register evolution & RHS gridfunction groups with MoL, so it knows

  group = CCTK_GroupIndex("BaikalETK::evol_variables");
  rhs = CCTK_GroupIndex("BaikalETK::evol_variables_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
}
"""
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("MoL_registration.c")] = outstr.replace("BaikalETK",ThornName)

    # Next register with the boundary conditions thorns.
    # PART 1: Set BC type to "none" for all variables
    # Since we choose NewRad boundary conditions, we must register all
    #   gridfunctions to have boundary type "none". This is because
    #   NewRad is seen by the rest of the Toolkit as a modification to the
    #   RHSs.

    # This code is based on Kranc's McLachlan/ML_BSSN/src/Boundaries.cc code.
    outstr = """
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"

// Set `none` boundary conditions on BSSN RHSs, as these are set via NewRad.
void BaikalETK_BoundaryConditions_evolved_gfs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
"""
    for gf in evol_gfs_list:
        outstr += """
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalETK::"""+gf+"""", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalETK::"""+gf+"""!");
"""
    outstr += """
}

// Set `flat` boundary conditions on BSSN constraints, similar to what Lean does.
void BaikalETK_BoundaryConditions_aux_gfs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

  CCTK_INT bndsize[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];

  GetBoundarySpecification(6, bndsize, is_internal, is_staggered, shiftout);

"""
    for gf in aux_gfs_list:
        outstr += """
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, bndsize[0], -1, "BaikalETK::"""+gf+"""", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalETK::"""+gf+"""!");
"""
    outstr += "}\n"

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("BoundaryConditions.c")] = outstr.replace("BaikalETK",ThornName)

    # PART 2: Set C code for calling NewRad BCs
    #   As explained in lean_public/LeanBSSNMoL/src/calc_bssn_rhs.F90,
    #   the function NewRad_Apply takes the following arguments:
    #   NewRad_Apply(cctkGH, var, rhs, var0, v0, radpower),
    #     which implement the boundary condition:
    #       var  =  var_at_infinite_r + u(r-var_char_speed*t)/r^var_radpower
    #  Obviously for var_radpower>0, var_at_infinite_r is the value of
    #    the variable at r->infinity. var_char_speed is the propagation
    #    speed at the outer boundary, and var_radpower is the radial
    #    falloff rate.

    outstr = """
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_NewRad(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

"""
    for gf in evol_gfs_list:
        var_at_infinite_r = "0.0"
        var_char_speed    = "1.0"
        var_radpower      = "1.0"

        if "alpha" in gf:
            var_at_infinite_r = "1.0"
            if LapseCondition == "OnePlusLog":
                var_char_speed = "sqrt(2.0)"
            else:
                pass # 1.0 (default) is fine
        if "aDD" in gf or "trK" in gf: # consistent with Lean code.
            var_radpower = "2.0"

        outstr += "  NewRad_Apply(cctkGH, "+gf+", "+gf.replace("GF","")+"_rhsGF, "+var_at_infinite_r+", "+var_char_speed+", "+var_radpower+");\n"
    outstr += "}\n"

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("BoundaryCondition_NewRad.c")] = outstr.replace("BaikalETK",ThornName)

    # First we convert from ADM to BSSN, as is required to convert initial data
    #    (given using) ADM quantities, to the BSSN evolved variables
    import BSSN.ADM_Numerical_Spherical_or_Cartesian_to_BSSNCurvilinear as atob
    IDhDD,IDaDD,IDtrK,IDvetU,IDbetU,IDalpha,IDcf,IDlambdaU = \
        atob.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear("Cartesian","DoNotOutputADMInputFunction",os.path.join(ThornName,"src"))

    # Store the original list of registered gridfunctions; we'll want to unregister
    #   all the *SphorCart* gridfunctions after we're finished with them below.
    orig_glb_gridfcs_list = []
    for gf in gri.glb_gridfcs_list:
        orig_glb_gridfcs_list.append(gf)

    # We ignore the return values for the following register_gridfunctions...() calls,
    #   as they are unused.
    gri.register_gridfunctions(                 "AUXEVOL", "alphaSphorCart")
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL", "betaSphorCartU")
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL", "BSphorCartU")
    ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "gammaSphorCartDD", "sym01")
    ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "KSphorCartDD", "sym01")

    # ADM to BSSN conversion, used for converting ADM initial data into a form readable by this thorn.
    # ADM to BSSN, Part 1: Set up function call and pointers to ADM gridfunctions
    outstr = """
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_ADM_to_BSSN(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL *alphaSphorCartGF = alp;
"""
    # It's ugly if we output code in the following ordering, so we'll first
    #   output to a string and then sort the string to beautify the code a bit.
    outstrtmp = []
    for i in range(3):
        outstrtmp.append("    CCTK_REAL *betaSphorCartU"+str(i)+"GF = beta"+chr(ord('x')+i)+";\n")
        outstrtmp.append("    CCTK_REAL *BSphorCartU"+str(i)+"GF = dtbeta"+chr(ord('x')+i)+";\n")
        for j in range(i,3):
            outstrtmp.append("    CCTK_REAL *gammaSphorCartDD"+str(i)+str(j)+"GF = g"+chr(ord('x')+i)+chr(ord('x')+j)+";\n")
            outstrtmp.append("    CCTK_REAL *KSphorCartDD"+str(i)+str(j)+"GF = k"+chr(ord('x')+i)+chr(ord('x')+j)+";\n")
    outstrtmp.sort()
    for line in outstrtmp:
        outstr += line

    # ADM to BSSN, Part 2: Set up ADM to BSSN conversions for BSSN gridfunctions that do not require
    #                      finite-difference derivatives (i.e., all gridfunctions except lambda^i (=Gamma^i
    #                      in non-covariant BSSN)):
    #                      h_{ij}, a_{ij}, trK, vet^i=beta^i,bet^i=B^i, cf (conformal factor), and alpha

    # Output finite difference stencils as inlined expressions.
    #   We do this instead of outputting as FD functions, as this function
    #   does not take long to compile, and we have already output all the
    #   FD functions to file, so if this one contains new FD functions,
    #   the compile will fail.
    par.set_parval_from_str("finite_difference::enable_FD_functions", False)

    all_but_lambdaU_expressions = [
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD00"), rhs=IDhDD[0][0]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD01"), rhs=IDhDD[0][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD02"), rhs=IDhDD[0][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD11"), rhs=IDhDD[1][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD12"), rhs=IDhDD[1][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD22"), rhs=IDhDD[2][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "aDD00"), rhs=IDaDD[0][0]),
        lhrh(lhs=gri.gfaccess("in_gfs", "aDD01"), rhs=IDaDD[0][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "aDD02"), rhs=IDaDD[0][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "aDD11"), rhs=IDaDD[1][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "aDD12"), rhs=IDaDD[1][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "aDD22"), rhs=IDaDD[2][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "trK"), rhs=IDtrK),
        lhrh(lhs=gri.gfaccess("in_gfs", "vetU0"), rhs=IDvetU[0]),
        lhrh(lhs=gri.gfaccess("in_gfs", "vetU1"), rhs=IDvetU[1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "vetU2"), rhs=IDvetU[2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "betU0"), rhs=IDbetU[0]),
        lhrh(lhs=gri.gfaccess("in_gfs", "betU1"), rhs=IDbetU[1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "betU2"), rhs=IDbetU[2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "alpha"), rhs=IDalpha),
        lhrh(lhs=gri.gfaccess("in_gfs", "cf"), rhs=IDcf)]

    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
    all_but_lambdaU_outC = fin.FD_outputC("returnstring", all_but_lambdaU_expressions, outCparams)
    outstr += lp.loop(["i2", "i1", "i0"], ["0", "0", "0"], ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                      ["1", "1", "1"], ["#pragma omp parallel for", "", ""], "    ", all_but_lambdaU_outC)

    # ADM to BSSN, Part 3: Set up ADM to BSSN conversions for BSSN gridfunctions defined from
    #                      finite-difference derivatives: lambda^i, which is Gamma^i in non-covariant BSSN):
    outstr += """
    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
"""

    path = os.path.join(ThornName,"src")
    BaikalETK_src_filelist = []
    for _root, _dirs, files in os.walk(path):  # _root, _dirs unused.
        for filename in files:
            BaikalETK_src_filelist.append(filename)
    BaikalETK_src_filelist.sort() # Sort the list in place.

    BSSN_FD_orders_output = []
    for filename in BaikalETK_src_filelist:
        if "BSSN_RHSs_" in filename:
            array = filename.replace(".","_").split("_")
            FDorder =  int(array[-2])
            if FDorder not in BSSN_FD_orders_output:
                BSSN_FD_orders_output.append(FDorder)
    BSSN_FD_orders_output.sort()

    for current_FD_order in BSSN_FD_orders_output:
        # Store original finite-differencing order:
        orig_FD_order = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
        # Set new finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", current_FD_order)

        outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
        lambdaU_expressions = [lhrh(lhs=gri.gfaccess("in_gfs","lambdaU0"),rhs=IDlambdaU[0]),
                               lhrh(lhs=gri.gfaccess("in_gfs","lambdaU1"),rhs=IDlambdaU[1]),
                               lhrh(lhs=gri.gfaccess("in_gfs","lambdaU2"),rhs=IDlambdaU[2])]
        lambdaU_expressions_FDout = fin.FD_outputC("returnstring",lambdaU_expressions, outCparams)

        lambdafile = "ADM_to_BSSN__compute_lambdaU_FD_order_"+str(current_FD_order)+".h"
        with open(os.path.join(ThornName,"src",lambdafile),"w") as file:
            file.write(lp.loop(["i2","i1","i0"],["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],
                       ["cctk_lsh[2]-cctk_nghostzones[2]","cctk_lsh[1]-cctk_nghostzones[1]","cctk_lsh[0]-cctk_nghostzones[0]"],
                       ["1","1","1"],["#pragma omp parallel for","",""],"",lambdaU_expressions_FDout))

        outstr += "    if(FD_order == "+str(current_FD_order)+") {\n"
        outstr += "        #include \""+lambdafile+"\"\n"
        outstr += "    }\n"
        # Restore original FD order
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", orig_FD_order)

    outstr += """
    ExtrapolateGammas(cctkGH,lambdaU0GF);
    ExtrapolateGammas(cctkGH,lambdaU1GF);
    ExtrapolateGammas(cctkGH,lambdaU2GF);
}
"""

    # Unregister the *SphorCartGF's.
    gri.glb_gridfcs_list = orig_glb_gridfcs_list

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("ADM_to_BSSN.c")] = outstr.replace("BaikalETK",ThornName)

    import BSSN.ADM_in_terms_of_BSSN as btoa
    import BSSN.BSSN_quantities as Bq

    btoa.ADM_in_terms_of_BSSN()
    Bq.BSSN_basic_tensors() # Gives us betaU & BU

    outstr = """
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_BSSN_to_ADM(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

"""
    btoa_lhrh = []
    for i in range(3):
        for j in range(i,3):
            btoa_lhrh.append(lhrh(lhs="g"+chr(ord('x')+i)+chr(ord('x')+j)+"[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                                  rhs=btoa.gammaDD[i][j]))
    for i in range(3):
        for j in range(i,3):
            btoa_lhrh.append(lhrh(lhs="k"+chr(ord('x')+i)+chr(ord('x')+j)+"[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                                  rhs=btoa.KDD[i][j]))
    btoa_lhrh.append(lhrh(lhs="alp[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",rhs=Bq.alpha))

    for i in range(3):
        btoa_lhrh.append(lhrh(lhs="beta"+chr(ord('x')+i)+"[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                                  rhs=Bq.betaU[i]))

    for i in range(3):
        btoa_lhrh.append(lhrh(lhs="dtbeta"+chr(ord('x')+i)+"[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                                  rhs=Bq.BU[i]))

    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
    bssn_to_adm_Ccode = fin.FD_outputC("returnstring",btoa_lhrh, outCparams)
    outstr += lp.loop(["i2","i1","i0"],["0","0","0"],["cctk_lsh[2]","cctk_lsh[1]","cctk_lsh[0]"],
                       ["1","1","1"],["#pragma omp parallel for","",""],"",bssn_to_adm_Ccode)

    outstr += "}\n"

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("BSSN_to_ADM.c")] = outstr.replace("BaikalETK",ThornName)

    ###########################
    ###########################
    # BSSN_RHSs and Ricci driver functions
    ###########################
    common_includes = """
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "SIMD/SIMD_intrinsics.h"
#include "finite_difference_functions.h"
"""
    common_preloop = """
    DECLARE_CCTK_ARGUMENTS;
    const CCTK_REAL NOSIMDinvdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const REAL_SIMD_ARRAY invdx0 = ConstSIMD(NOSIMDinvdx0);
    const CCTK_REAL NOSIMDinvdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const REAL_SIMD_ARRAY invdx1 = ConstSIMD(NOSIMDinvdx1);
    const CCTK_REAL NOSIMDinvdx2 = 1.0/CCTK_DELTA_SPACE(2);
    const REAL_SIMD_ARRAY invdx2 = ConstSIMD(NOSIMDinvdx2);
"""

    # Output the driver code for computing the Ricci tensor:
    outstr = common_includes
    for order in BSSN_FD_orders_output:
        outstr += """extern void """+ThornName+"_"+"BSSN_Ricci_FD_order_"+str(order)+"(CCTK_ARGUMENTS);\n"
    outstr += """
void BaikalETK_driver_pt1_BSSN_Ricci(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    const CCTK_INT *FD_order = CCTK_ParameterGet("FD_order","BaikalETK",NULL);
"""
    for order in BSSN_FD_orders_output:
        outstr += "    if(*FD_order == "+str(order)+") {\n"
        outstr += "        "+ThornName+"_"+"BSSN_Ricci_FD_order_"+str(order)+"(CCTK_PASS_CTOC);\n"
        outstr += "    }\n"
    outstr += "} // END FUNCTION\n"
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("driver_pt1_BSSN_Ricci.c")] = outstr.replace("BaikalETK",ThornName)

    # Create functions for the largest C kernels (BSSN RHSs and Ricci) and output
    #    the .h files to .c files with function wrappers; delete original .h files
    for filename in BaikalETK_src_filelist:
        if ("BSSN_Ricci_FD_order_") in filename and (".h" in filename):
            outstr = common_includes + "void BaikalETK_"+filename.replace(".h","")+"(CCTK_ARGUMENTS) {\n"
            outstr += common_preloop
            with open(os.path.join(path,filename), "r") as currfile:
                outstr += currfile.read()
            # Now that we've inserted the contents of the kernel into this file,
            #     we delete the file containing the kernel
            os.remove(os.path.join(path,filename))
            outstr += "} // END FUNCTION\n"
            # Add C code string to dictionary (Python dictionaries are immutable)
            Csrcdict[append_to_make_code_defn_list(filename.replace(".h",".c"))] = outstr.replace("BaikalETK",ThornName)
    ###########################
    # Output BSSN RHSs driver function
    outstr = common_includes
    for filename in BaikalETK_src_filelist:
        if ("BSSN_RHSs_" in filename) and (".h" in filename):
            outstr += """extern void """ + ThornName+"_"+filename.replace(".h", "(CCTK_ARGUMENTS);") + "\n"

    outstr += """
void BaikalETK_driver_pt2_BSSN_RHSs(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;

    const CCTK_INT *FD_order = CCTK_ParameterGet("FD_order","BaikalETK",NULL);

"""
    for filename in BaikalETK_src_filelist:
        if ("BSSN_RHSs_" in filename) and (".h" in filename):
            array = filename.replace(".", "_").split("_")
            outstr += "    if(*FD_order == " + str(array[-2]) + ") {\n"
            outstr += "        " + ThornName+"_"+filename.replace(".h", "(CCTK_PASS_CTOC);") + "\n"
            outstr += "    }\n"
    outstr += "} // END FUNCTION\n"
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("driver_pt2_BSSN_RHSs.c")] = outstr.replace("BaikalETK", ThornName)

    def SIMD_declare_C_params():
        SIMD_declare_C_params_str = ""
        for i in range(len(par.glb_Cparams_list)):
            # keep_param is a boolean indicating whether we should accept or reject
            #    the parameter. singleparstring will contain the string indicating
            #    the variable type.
            keep_param, singleparstring = ccl.keep_param__return_type(par.glb_Cparams_list[i])

            if (keep_param) and ("CCTK_REAL" in singleparstring):
                parname = par.glb_Cparams_list[i].parname
                SIMD_declare_C_params_str += "    const "+singleparstring + "*NOSIMD"+parname+\
                                             " = CCTK_ParameterGet(\""+parname+"\",\"BaikalETK\",NULL);\n"
                SIMD_declare_C_params_str += "    const REAL_SIMD_ARRAY "+parname+" = ConstSIMD(*NOSIMD"+parname+");\n"
        return SIMD_declare_C_params_str

    # Create functions for the largest C kernels (BSSN RHSs and Ricci) and output
    #    the .h files to .c files with function wrappers; delete original .h files
    path = os.path.join(ThornName, "src")
    for filename in BaikalETK_src_filelist:
        if ("BSSN_RHSs_" in filename) and (".h" in filename):
            outstr = common_includes + "void BaikalETK_"+filename.replace(".h","")+"(CCTK_ARGUMENTS) {\n"
            outstr += common_preloop+SIMD_declare_C_params()
            with open(os.path.join(path,filename), "r") as currfile:
                outstr += currfile.read()
            # Now that we've inserted the contents of the kernel into this file,
            #     we delete the file containing the kernel
            os.remove(os.path.join(path,filename))
            outstr += "} // END FUNCTION\n"
            # Add C code string to dictionary (Python dictionaries are immutable)
            Csrcdict[append_to_make_code_defn_list(filename.replace(".h",".c"))] = outstr.replace("BaikalETK",ThornName)

    # Next, the driver for enforcing detgammabar = detgammahat constraint:
    outstr = common_includes + """
void BaikalETK_enforce_detgammahat_constraint(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    #include "enforcedetgammabar_constraint.h"
}
"""
    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("driver_enforcedetgammabar_constraint.c")] = \
        outstr.replace("BaikalETK", ThornName)

    # Next, the driver for computing the BSSN Hamiltonian & momentum constraints
    outstr = """
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "SIMD/SIMD_intrinsics.h" // Contains needed definition of REAL_SIMD_ARRAY
#include "finite_difference_functions.h"

void BaikalETK_BSSN_constraints(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
"""
    for filename in BaikalETK_src_filelist:
        if "BSSN_constraints_" in filename:
            array = filename.replace(".","_").split("_")
            outstr += "    if(FD_order == "+str(array[-2])+") {\n"
            outstr += "        #include \""+filename+"\"\n"
            outstr += "    }\n"
    outstr += "}\n"

    # Add C code string to dictionary (Python dictionaries are immutable)
    Csrcdict[append_to_make_code_defn_list("driver_BSSN_constraints.c")] = outstr.replace("BaikalETK",ThornName)

    if enable_stress_energy_source_terms == True:
        # Declare T4DD as a set of gridfunctions. These won't
        #    actually appear in interface.ccl, as interface.ccl
        #    was set above. Thus before calling the code output
        #    by FD_outputC(), we'll have to set pointers
        #    to the actual gridfunctions they reference.
        #    (In this case the eTab's.)
        T4DD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","T4DD","sym01",DIM=4)
        import BSSN.ADMBSSN_tofrom_4metric as AB4m
        AB4m.g4UU_ito_BSSN_or_ADM("BSSN")

        T4UUraised = ixp.zerorank2(DIM=4)
        for mu in range(4):
            for nu in range(4):
                for delta in range(4):
                    for gamma in range(4):
                        T4UUraised[mu][nu] += AB4m.g4UU[mu][delta]*AB4m.g4UU[nu][gamma]*T4DD[delta][gamma]

        T4UU_expressions = [
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU00"),rhs=T4UUraised[0][0]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU01"),rhs=T4UUraised[0][1]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU02"),rhs=T4UUraised[0][2]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU03"),rhs=T4UUraised[0][3]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU11"),rhs=T4UUraised[1][1]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU12"),rhs=T4UUraised[1][2]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU13"),rhs=T4UUraised[1][3]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU22"),rhs=T4UUraised[2][2]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU23"),rhs=T4UUraised[2][3]),
            lhrh(lhs=gri.gfaccess("in_gfs","T4UU33"),rhs=T4UUraised[3][3])]

        outCparams = "outCverbose=False,includebraces=False,preindent=2,enable_SIMD=True"
        T4UUstr = fin.FD_outputC("returnstring",T4UU_expressions, outCparams)
        T4UUstr_loop = lp.loop(["i2","i1","i0"],["0","0","0"],["cctk_lsh[2]","cctk_lsh[1]","cctk_lsh[0]"],
                               ["1","1","SIMD_width"],["#pragma omp parallel for","",""],"",T4UUstr)

        outstr = """
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SIMD/SIMD_intrinsics.h"

void BaikalETK_driver_BSSN_T4UU(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL *restrict T4DD00GF = eTtt;
    const CCTK_REAL *restrict T4DD01GF = eTtx;
    const CCTK_REAL *restrict T4DD02GF = eTty;
    const CCTK_REAL *restrict T4DD03GF = eTtz;
    const CCTK_REAL *restrict T4DD11GF = eTxx;
    const CCTK_REAL *restrict T4DD12GF = eTxy;
    const CCTK_REAL *restrict T4DD13GF = eTxz;
    const CCTK_REAL *restrict T4DD22GF = eTyy;
    const CCTK_REAL *restrict T4DD23GF = eTyz;
    const CCTK_REAL *restrict T4DD33GF = eTzz;
"""+T4UUstr_loop+"""
}\n"""

        # Add C code string to dictionary (Python dictionaries are immutable)
        Csrcdict[append_to_make_code_defn_list("driver_BSSN_T4UU.c")] = outstr.replace("BaikalETK",ThornName)