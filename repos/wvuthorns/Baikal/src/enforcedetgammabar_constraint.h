#pragma omp parallel for
for (int i2 = 0; i2 < cctk_lsh[2]; i2++) {
    #include "rfm_files/rfm_struct__read2.h"
    for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
        #include "rfm_files/rfm_struct__read1.h"
        for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
            #include "rfm_files/rfm_struct__read0.h"
            /*
             * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
             */
            const double hDD00 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD01 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD02 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD11 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD12 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            const double hDD22 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
            /*
             * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
             */
            const double FDPart3_0 = hDD00 + 1;
            const double FDPart3_1 = hDD22 + 1;
            const double FDPart3_2 = hDD11 + 1;
            const double FDPart3_3 = cbrt(fabs(1)/(FDPart3_0*FDPart3_1*FDPart3_2 - FDPart3_0*((hDD12)*(hDD12)) - FDPart3_1*((hDD01)*(hDD01)) - FDPart3_2*((hDD02)*(hDD02)) + 2*hDD01*hDD02*hDD12));
            hDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_0*FDPart3_3 - 1;
            hDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_3*hDD01;
            hDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_3*hDD02;
            hDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_2*FDPart3_3 - 1;
            hDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_3*hDD12;
            hDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_3 - 1;
            
        } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
    } // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
} // END LOOP: for (int i2 = 0; i2 < cctk_lsh[2]; i2++)
