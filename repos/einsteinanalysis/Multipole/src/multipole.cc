#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

#include <cstdio>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "multipole.hh"
#include "interpolate.hh"
#include "utils.hh"
#include "sphericalharmonic.hh"

using namespace std;

static void fill_variable(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg != 0);

  vector<Multipole::variable_desc> &vs = *(vector<Multipole::variable_desc> * ) callback_arg;

  Multipole::variable_desc v;

  v.index = idx;

  // Default values if there is no option string or if the options are
  // not present
  v.imag_index = -1;
  v.spin_weight = 0;
  v.name = string(CCTK_VarName(v.index));

  if (optstring != 0)
  {
    int table = Util_TableCreateFromString(optstring);

    if (table >= 0)
    {
      const int buffer_length = 256;
      char buffer[buffer_length];
      
      Util_TableGetInt(table, &v.spin_weight , "sw");
      if (Util_TableGetString(table, buffer_length, buffer , "cmplx") >= 0)
      {
        v.imag_index = CCTK_VarIndex(buffer);
      }
      if (Util_TableGetString(table, buffer_length, buffer , "name") >= 0)
      {
        v.name = string(buffer);
      }

      const int ierr = Util_TableDestroy(table);
      if (ierr)
      {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not destroy table: %d", ierr);
      }
    }
  }
  vs.push_back(v);
}

static void parse_variables_string(const string &var_string,
                     vector<Multipole::variable_desc> &vars)
{
  int ierr = CCTK_TraverseString(var_string.c_str(), fill_variable, &vars, CCTK_GROUP_OR_VAR);
  assert(ierr >= 0);
}

static void output_modes(CCTK_ARGUMENTS, const vector<Multipole::variable_desc> &vars,
                     const CCTK_REAL radii[], const Multipole::mode_array& modes)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_ascii)
  {
    if (CCTK_MyProc(cctkGH) == 0)
    {
      for (int v = 0; v < modes.get_nvars(); v++)
      {
        for(int i=0; i<modes.get_nradii(); i++)
        {
          const CCTK_REAL rad = radii[i];
          for (int l=0; l <= modes.get_lmax(); l++)
          {
            for (int m=-l; m <= l; m++)
            {
              ostringstream name;
              name << "mp_" << vars[v].name << "_l" << l << "_m" << m <<
                "_r" << setiosflags(ios::fixed) << setprecision(2) <<  rad << ".asc";
              Multipole_OutputComplexToFile(CCTK_PASS_CTOC, name.str(),
                                            modes(v, i, l, m, 0), modes(v, i, l, m, 1));
            }
          }
        }
      }
    }
  }
  if (output_hdf5)
  {
    if (CCTK_MyProc(cctkGH) == 0)
    {
      Multipole_OutputComplexToH5File(CCTK_PASS_CTOC, vars, radii, modes);
    }
  }
}

static void output_1D(CCTK_ARGUMENTS, const Multipole::variable_desc &v, CCTK_REAL rad,
                      const vector<CCTK_REAL> &th, const vector<CCTK_REAL> &ph,
                      const vector<CCTK_REAL> &real, const vector<CCTK_REAL> &imag)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (CCTK_MyProc(cctkGH) == 0 && output_ascii)
  {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0)
    {
      ostringstream real_base;
      real_base << "mp_" << string(CCTK_VarName(v.index)) << "_r" << setiosflags(ios::fixed) << setprecision(2) << rad;
      Multipole_Output1D(CCTK_PASS_CTOC, real_base.str()+string(".th.asc"), th, ph, mp_theta, real);
      Multipole_Output1D(CCTK_PASS_CTOC, real_base.str()+string(".ph.asc"), th, ph, mp_phi, real);

      if (v.imag_index != -1)
      {
        ostringstream imag_base;
        imag_base << "mp_" << string(CCTK_VarName(v.imag_index)) << "_r" << setiosflags(ios::fixed) << setprecision(2) << rad;
        Multipole_Output1D(CCTK_PASS_CTOC, imag_base.str()+string(".th.asc"), th, ph, mp_theta, imag);
        Multipole_Output1D(CCTK_PASS_CTOC, imag_base.str()+string(".ph.asc"), th, ph, mp_phi, imag);
      }
    }
  }
}

bool int_in_array(int a, const vector<int> &array)
{
  for (size_t i = 0; i < array.size(); i++)
  {
    if (array[i] == a)
      return true;
  }
  return false;
}

int find_int_in_array(int a, const vector<int> &array)
{
  for (size_t i = 0; i < array.size(); i++)
  {
    if (array[i] == a)
      return i;
  }
  return -1;
}


static
void get_spin_weights(const vector<Multipole::variable_desc> &vars, vector<int> &spin_weights)
{
  for (size_t i = 0; i < vars.size(); i++)
  {
    if (!int_in_array(vars[i].spin_weight, spin_weights))
    {
      spin_weights.push_back(vars[i].spin_weight);
    }
  }
}

// For backward compatibility we allow the user to set l_mode instead
// of l_max, but if it is left at the default of -1, l_max is used.
static int get_l_max()
{
  DECLARE_CCTK_PARAMETERS;
  return l_mode == -1 ? l_max : l_mode;
}

static
void setup_harmonics(const vector<int> &spin_weights,
                     int lmax,
                     vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph, int array_size,
                     vector<vector<vector<vector<CCTK_REAL>>>> &reY,
                     vector<vector<vector<vector<CCTK_REAL>>>> &imY)
{
  th.resize(array_size);
  ph.resize(array_size);

  reY.resize(spin_weights.size());
  imY.resize(spin_weights.size());
  for (size_t si = 0; si < spin_weights.size(); si++)
  {
    int sw = spin_weights[si];

    reY[si].resize(lmax+1);
    imY[si].resize(lmax+1);
    for (int l = 0; l <= lmax; l++)
    {
      reY[si][l].resize(2*lmax+1);
      imY[si][l].resize(2*lmax+1);
      for (int m = -l; m <= lmax; m++)
      {
        reY[si][l][m+l].resize(array_size);
        imY[si][l][m+l].resize(array_size);
        Multipole_HarmonicSetup(sw, l, m, th, ph,
                                reY[si][l][m+l], imY[si][l][m+l]);
      }
    }
  }
}

extern "C" void Multipole_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if (l_mode != -1)
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter l_mode is deprecated. Use l_max instead.  For compatibility, l_max = l_mode is being used.");
  }

  if (!CCTK_Equals(mode_type, "deprecated"))
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter mode_type is deprecated and is no longer used.  All modes will be computed.");
  }

  if (l_min != -1)
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter l_min is deprecated and is no longer used.  Modes from l = 0 will be computed.");
  }

  if (m_mode != -100)
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter m_mode is deprecated. All m modes will be computed.");
  }
}

extern "C" void Multipole_Calc(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  static vector<CCTK_REAL> xs, ys, zs;
  static vector<CCTK_REAL> xhat, yhat, zhat;
  static vector<CCTK_REAL> th, ph;
  static vector<CCTK_REAL> real, imag;
  static vector<vector<vector<vector<CCTK_REAL>>>> reY;
  static vector<vector<vector<vector<CCTK_REAL>>>> imY;
  static vector<Multipole::variable_desc> vars;
  static vector<int> spin_weights;

  static bool initialized = false;

  const int array_size=(ntheta+1)*(nphi+1);

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  int lmax = get_l_max();

  if (!initialized)
  {
    real.resize(array_size);
    imag.resize(array_size);
    th.resize(array_size);
    ph.resize(array_size);
    xs.resize(array_size);
    ys.resize(array_size);
    zs.resize(array_size);
    xhat.resize(array_size);
    yhat.resize(array_size);
    zhat.resize(array_size);
  
    parse_variables_string(string(variables), vars);
    get_spin_weights(vars, spin_weights);
    Multipole_CoordSetup(xhat, yhat, zhat, th, ph);
    setup_harmonics(spin_weights, lmax, th, ph,
                    array_size, reY, imY);
    initialized = true;
  }

  Multipole::mode_array modes(vars.size(), nradii, lmax);
  for (size_t v = 0; v < vars.size(); v++)
  {
    //assert(vars[v].spin_weight == -2);

    int si = find_int_in_array(vars[v].spin_weight, spin_weights);
    assert(si != -1);

    for(int i=0; i<nradii; i++)
    {
      // Compute x^i = r * \hat x^i
      Multipole_ScaleCartesian(radius[i], xhat, yhat, zhat, xs, ys, zs);
      
      // Interpolate Psi4r and Psi4i
      Multipole_Interp(CCTK_PASS_CTOC, xs, ys, zs, vars[v].index, vars[v].imag_index, 
                       real, imag);
      for (int l=0; l <= lmax; l++)
      {
        for (int m=-l; m <= l; m++)
        {
          // Integrate sYlm (real + i imag) over the sphere at radius r
          Multipole_Integrate(reY[si][l][m+l], imY[si][l][m+l],
                              real, imag, th, ph,
                              &modes(v, i, l, m, 0), &modes(v, i, l, m, 1));
    
        }//loop over m
      }//loop over l
      output_1D(CCTK_PASS_CTOC, vars[v], radius[i], th, ph, real, imag);
    }//loop over radii
  }//loop over variables
  output_modes(CCTK_PASS_CTOC, vars, radius, modes);
}
