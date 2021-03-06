#include <stdio.h>
#include <stdlib.h>
#include "NRPyEOS.h"

int main(int argc, char **argv) {

  // Step 0: Check for correct usage
  if( argc != 2 ) {
    fprintf(stderr,"(NRPyEOS - minimal) Correct usage is ./minimal eos_file_path\n");
    exit(1);
  }

  // Step 1: Initialize the EOS struct
  NRPyEOS_params eos_params;
  NRPyEOS_readtable_set_EOS_params(argv[1],&eos_params);

  // Step 2: Perform one interpolation
  const double rho = 1e-6;
  const double Y_e = 0.4;
  const double T   = 1.1;
  double P, eps;
  NRPyEOS_P_and_eps_from_rho_Ye_T(&eos_params,rho,Y_e,T,&P,&eps);

  // Step 3: Print information
  printf("(NRPyEOS) Density    : %.15e\n",rho);
  printf("(NRPyEOS) e- fraction: %.15e\n",Y_e);
  printf("(NRPyEOS) Temperature: %.15e\n",T);
  printf("(NRPyEOS) Pressure   : %.15e\n",P);
  printf("(NRPyEOS) Energy     : %.15e\n",eps);

  // Step 4: Free memory
  NRPyEOS_free_memory(&eos_params);

  // All done!
  return 0;
}
