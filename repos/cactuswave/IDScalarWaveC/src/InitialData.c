 /*@@
   @file    InitialData.c
   @date
   @author  Werner Benger
   @desc
            Initial data for the 3D Wave Equation
            Derived from Tom Goodale
   @enddesc
   @version $Id$
 @@*/

#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#ifndef DECLARE_CCTK_ARGUMENTS_CHECKED
#define DECLARE_CCTK_ARGUMENTS_CHECKED(func) DECLARE_CCTK_ARGUMENTS
#endif

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_IDScalarWaveC_InitialData_c)

#define SQR(val) ((val) * (val))

void IDScalarWaveC_InitialData(CCTK_ARGUMENTS);


 /*@@
   @routine    IDScalarWaveC_InitialData
   @date
   @author     Tom Goodale
   @desc
               Set up initial data for the wave equation
   @enddesc
   @history
   @hdate Mon Oct 11 11:48:03 1999 @hauthor Werner Benger
   @hdesc  Converted to C++
   @hdate Mon Oct 11 11:48:20 1999 @hauthor Tom Goodale
   @hdesc Added the rest of the initial data.
   @hdate Thu Feb 17 09:22:01 2000 @hauthor Tom Goodale
   @hdesc Converted to C
   @endhistory

@@*/

void IDScalarWaveC_InitialData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_CHECKED(IDScalarWaveC_InitialData)
  DECLARE_CCTK_PARAMETERS

  int i,j,k;

  CCTK_REAL dt;
  CCTK_REAL omega;
  int vindex;
  CCTK_REAL X, Y, Z, R;
  CCTK_REAL pi;

  dt = CCTK_DELTA_TIME;

  if(CCTK_Equals(initial_data, "plane"))
  {
    omega = sqrt(SQR(kx)+SQR(ky)+SQR(kz));

    for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
      {
        for(i=0; i<cctk_lsh[0]; i++)
        {
          vindex =  CCTK_GFINDEX3D(cctkGH,i,j,k);

          phi[vindex]     = amplitude*cos(kx*x[vindex]+ky*y[vindex]+kz*z[vindex]+omega*cctk_time);
          phi_p[vindex] = amplitude*cos(kx*x[vindex]+ky*y[vindex]+kz*z[vindex]+omega*(cctk_time-dt));
        }
      }
    }
  }
  else if(CCTK_Equals(initial_data, "gaussian"))
  {
    for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
      {
        for(i=0; i<cctk_lsh[0]; i++)
        {
          vindex =  CCTK_GFINDEX3D(cctkGH,i,j,k);

          X = x[vindex];
          Y = y[vindex];
          Z = z[vindex];

          R = sqrt(X*X + Y*Y + Z*Z);

          phi[vindex] = amplitude*exp( - SQR( (R - radius) / sigma ) );

          if (R == 0.0)
          {
            phi_p[vindex] = amplitude*(1.0 - 2.0*dt*dt/sigma/sigma)*exp(-dt*dt/sigma/sigma);
          }
          else
          {
            phi_p[vindex] = amplitude/2.0*(R-dt)/R*
              exp( - SQR( (R - radius - dt)/ sigma ) )
              + amplitude/2.0*(R+dt)/R*
              exp( - SQR( (R - radius + dt)/ sigma ) );
          }
        }
      }
    }
  }
  else if(CCTK_Equals(initial_data, "box"))
  {
    pi = 4.0*atan(1.0);
    omega = sqrt(SQR(kx)+SQR(ky)+SQR(kz));

    for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
      {
        for(i=0; i<cctk_lsh[0]; i++)
        {
          vindex =  CCTK_GFINDEX3D(cctkGH,i,j,k);

          phi[vindex] = amplitude*sin(kx*(x[vindex]-0.5)*pi)*
                                 sin(ky*(y[vindex]-0.5)*pi)*
                                 sin(kz*(z[vindex]-0.5)*pi)*
                                 cos(omega*cctk_time*pi);

          phi_p[vindex] = amplitude*sin(kx*(x[vindex]-0.5)*pi)*
                                       sin(ky*(y[vindex]-0.5)*pi)*
                                       sin(kz*(z[vindex]-0.5)*pi)*
                                       cos(omega*(cctk_time-dt)*pi);
        }
      }
    }
  }
  else if (CCTK_Equals(initial_data, "none"))
  {
    for(k=0; k<cctk_lsh[2]; k++)
    {
      for(j=0; j<cctk_lsh[1]; j++)
      {
        for(i=0; i<cctk_lsh[0]; i++)
        {
          vindex =  CCTK_GFINDEX3D(cctkGH,i,j,k);

          phi[vindex] = phi_p[vindex] = 0.0;
        }
      }
    }
  }
}
