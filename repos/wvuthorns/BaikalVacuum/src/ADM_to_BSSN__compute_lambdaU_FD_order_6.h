#pragma omp parallel for
for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++) {
    for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++) {
        for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++) {
               /*
                * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                */
               const double hDD00_i0_i1_i2m3 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-3)];
               const double hDD00_i0_i1_i2m2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
               const double hDD00_i0_i1_i2m1 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
               const double hDD00_i0_i1m3_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1-3,i2)];
               const double hDD00_i0_i1m2_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
               const double hDD00_i0_i1m1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
               const double hDD00_i0m3_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0-3,i1,i2)];
               const double hDD00_i0m2_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
               const double hDD00_i0m1_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
               const double hDD00 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
               const double hDD00_i0p1_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
               const double hDD00_i0p2_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
               const double hDD00_i0p3_i1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0+3,i1,i2)];
               const double hDD00_i0_i1p1_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
               const double hDD00_i0_i1p2_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
               const double hDD00_i0_i1p3_i2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1+3,i2)];
               const double hDD00_i0_i1_i2p1 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
               const double hDD00_i0_i1_i2p2 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
               const double hDD00_i0_i1_i2p3 = hDD00GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+3)];
               const double hDD01_i0_i1_i2m3 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-3)];
               const double hDD01_i0_i1_i2m2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
               const double hDD01_i0_i1_i2m1 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
               const double hDD01_i0_i1m3_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1-3,i2)];
               const double hDD01_i0_i1m2_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
               const double hDD01_i0_i1m1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
               const double hDD01_i0m3_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0-3,i1,i2)];
               const double hDD01_i0m2_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
               const double hDD01_i0m1_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
               const double hDD01 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
               const double hDD01_i0p1_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
               const double hDD01_i0p2_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
               const double hDD01_i0p3_i1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0+3,i1,i2)];
               const double hDD01_i0_i1p1_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
               const double hDD01_i0_i1p2_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
               const double hDD01_i0_i1p3_i2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1+3,i2)];
               const double hDD01_i0_i1_i2p1 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
               const double hDD01_i0_i1_i2p2 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
               const double hDD01_i0_i1_i2p3 = hDD01GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+3)];
               const double hDD02_i0_i1_i2m3 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-3)];
               const double hDD02_i0_i1_i2m2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
               const double hDD02_i0_i1_i2m1 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
               const double hDD02_i0_i1m3_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1-3,i2)];
               const double hDD02_i0_i1m2_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
               const double hDD02_i0_i1m1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
               const double hDD02_i0m3_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0-3,i1,i2)];
               const double hDD02_i0m2_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
               const double hDD02_i0m1_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
               const double hDD02 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
               const double hDD02_i0p1_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
               const double hDD02_i0p2_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
               const double hDD02_i0p3_i1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0+3,i1,i2)];
               const double hDD02_i0_i1p1_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
               const double hDD02_i0_i1p2_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
               const double hDD02_i0_i1p3_i2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1+3,i2)];
               const double hDD02_i0_i1_i2p1 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
               const double hDD02_i0_i1_i2p2 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
               const double hDD02_i0_i1_i2p3 = hDD02GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+3)];
               const double hDD11_i0_i1_i2m3 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-3)];
               const double hDD11_i0_i1_i2m2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
               const double hDD11_i0_i1_i2m1 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
               const double hDD11_i0_i1m3_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1-3,i2)];
               const double hDD11_i0_i1m2_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
               const double hDD11_i0_i1m1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
               const double hDD11_i0m3_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0-3,i1,i2)];
               const double hDD11_i0m2_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
               const double hDD11_i0m1_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
               const double hDD11 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
               const double hDD11_i0p1_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
               const double hDD11_i0p2_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
               const double hDD11_i0p3_i1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0+3,i1,i2)];
               const double hDD11_i0_i1p1_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
               const double hDD11_i0_i1p2_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
               const double hDD11_i0_i1p3_i2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1+3,i2)];
               const double hDD11_i0_i1_i2p1 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
               const double hDD11_i0_i1_i2p2 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
               const double hDD11_i0_i1_i2p3 = hDD11GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+3)];
               const double hDD12_i0_i1_i2m3 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-3)];
               const double hDD12_i0_i1_i2m2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
               const double hDD12_i0_i1_i2m1 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
               const double hDD12_i0_i1m3_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1-3,i2)];
               const double hDD12_i0_i1m2_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
               const double hDD12_i0_i1m1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
               const double hDD12_i0m3_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0-3,i1,i2)];
               const double hDD12_i0m2_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
               const double hDD12_i0m1_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
               const double hDD12 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
               const double hDD12_i0p1_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
               const double hDD12_i0p2_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
               const double hDD12_i0p3_i1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0+3,i1,i2)];
               const double hDD12_i0_i1p1_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
               const double hDD12_i0_i1p2_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
               const double hDD12_i0_i1p3_i2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1+3,i2)];
               const double hDD12_i0_i1_i2p1 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
               const double hDD12_i0_i1_i2p2 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
               const double hDD12_i0_i1_i2p3 = hDD12GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+3)];
               const double hDD22_i0_i1_i2m3 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-3)];
               const double hDD22_i0_i1_i2m2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-2)];
               const double hDD22_i0_i1_i2m1 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2-1)];
               const double hDD22_i0_i1m3_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1-3,i2)];
               const double hDD22_i0_i1m2_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1-2,i2)];
               const double hDD22_i0_i1m1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1-1,i2)];
               const double hDD22_i0m3_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0-3,i1,i2)];
               const double hDD22_i0m2_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0-2,i1,i2)];
               const double hDD22_i0m1_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0-1,i1,i2)];
               const double hDD22 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
               const double hDD22_i0p1_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0+1,i1,i2)];
               const double hDD22_i0p2_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0+2,i1,i2)];
               const double hDD22_i0p3_i1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0+3,i1,i2)];
               const double hDD22_i0_i1p1_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1+1,i2)];
               const double hDD22_i0_i1p2_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1+2,i2)];
               const double hDD22_i0_i1p3_i2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1+3,i2)];
               const double hDD22_i0_i1_i2p1 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+1)];
               const double hDD22_i0_i1_i2p2 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+2)];
               const double hDD22_i0_i1_i2p3 = hDD22GF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2+3)];
               const double FDPart1_Rational_3_4 = 3.0/4.0;
               const double FDPart1_Rational_3_20 = 3.0/20.0;
               const double FDPart1_Rational_1_60 = 1.0/60.0;
               const double hDD_dD000 = invdx0*(FDPart1_Rational_1_60*(-hDD00_i0m3_i1_i2 + hDD00_i0p3_i1_i2) + FDPart1_Rational_3_20*(hDD00_i0m2_i1_i2 - hDD00_i0p2_i1_i2) + FDPart1_Rational_3_4*(-hDD00_i0m1_i1_i2 + hDD00_i0p1_i1_i2));
               const double hDD_dD001 = invdx1*(FDPart1_Rational_1_60*(-hDD00_i0_i1m3_i2 + hDD00_i0_i1p3_i2) + FDPart1_Rational_3_20*(hDD00_i0_i1m2_i2 - hDD00_i0_i1p2_i2) + FDPart1_Rational_3_4*(-hDD00_i0_i1m1_i2 + hDD00_i0_i1p1_i2));
               const double hDD_dD002 = invdx2*(FDPart1_Rational_1_60*(-hDD00_i0_i1_i2m3 + hDD00_i0_i1_i2p3) + FDPart1_Rational_3_20*(hDD00_i0_i1_i2m2 - hDD00_i0_i1_i2p2) + FDPart1_Rational_3_4*(-hDD00_i0_i1_i2m1 + hDD00_i0_i1_i2p1));
               const double hDD_dD010 = invdx0*(FDPart1_Rational_1_60*(-hDD01_i0m3_i1_i2 + hDD01_i0p3_i1_i2) + FDPart1_Rational_3_20*(hDD01_i0m2_i1_i2 - hDD01_i0p2_i1_i2) + FDPart1_Rational_3_4*(-hDD01_i0m1_i1_i2 + hDD01_i0p1_i1_i2));
               const double hDD_dD011 = invdx1*(FDPart1_Rational_1_60*(-hDD01_i0_i1m3_i2 + hDD01_i0_i1p3_i2) + FDPart1_Rational_3_20*(hDD01_i0_i1m2_i2 - hDD01_i0_i1p2_i2) + FDPart1_Rational_3_4*(-hDD01_i0_i1m1_i2 + hDD01_i0_i1p1_i2));
               const double hDD_dD012 = invdx2*(FDPart1_Rational_1_60*(-hDD01_i0_i1_i2m3 + hDD01_i0_i1_i2p3) + FDPart1_Rational_3_20*(hDD01_i0_i1_i2m2 - hDD01_i0_i1_i2p2) + FDPart1_Rational_3_4*(-hDD01_i0_i1_i2m1 + hDD01_i0_i1_i2p1));
               const double hDD_dD020 = invdx0*(FDPart1_Rational_1_60*(-hDD02_i0m3_i1_i2 + hDD02_i0p3_i1_i2) + FDPart1_Rational_3_20*(hDD02_i0m2_i1_i2 - hDD02_i0p2_i1_i2) + FDPart1_Rational_3_4*(-hDD02_i0m1_i1_i2 + hDD02_i0p1_i1_i2));
               const double hDD_dD021 = invdx1*(FDPart1_Rational_1_60*(-hDD02_i0_i1m3_i2 + hDD02_i0_i1p3_i2) + FDPart1_Rational_3_20*(hDD02_i0_i1m2_i2 - hDD02_i0_i1p2_i2) + FDPart1_Rational_3_4*(-hDD02_i0_i1m1_i2 + hDD02_i0_i1p1_i2));
               const double hDD_dD022 = invdx2*(FDPart1_Rational_1_60*(-hDD02_i0_i1_i2m3 + hDD02_i0_i1_i2p3) + FDPart1_Rational_3_20*(hDD02_i0_i1_i2m2 - hDD02_i0_i1_i2p2) + FDPart1_Rational_3_4*(-hDD02_i0_i1_i2m1 + hDD02_i0_i1_i2p1));
               const double hDD_dD110 = invdx0*(FDPart1_Rational_1_60*(-hDD11_i0m3_i1_i2 + hDD11_i0p3_i1_i2) + FDPart1_Rational_3_20*(hDD11_i0m2_i1_i2 - hDD11_i0p2_i1_i2) + FDPart1_Rational_3_4*(-hDD11_i0m1_i1_i2 + hDD11_i0p1_i1_i2));
               const double hDD_dD111 = invdx1*(FDPart1_Rational_1_60*(-hDD11_i0_i1m3_i2 + hDD11_i0_i1p3_i2) + FDPart1_Rational_3_20*(hDD11_i0_i1m2_i2 - hDD11_i0_i1p2_i2) + FDPart1_Rational_3_4*(-hDD11_i0_i1m1_i2 + hDD11_i0_i1p1_i2));
               const double hDD_dD112 = invdx2*(FDPart1_Rational_1_60*(-hDD11_i0_i1_i2m3 + hDD11_i0_i1_i2p3) + FDPart1_Rational_3_20*(hDD11_i0_i1_i2m2 - hDD11_i0_i1_i2p2) + FDPart1_Rational_3_4*(-hDD11_i0_i1_i2m1 + hDD11_i0_i1_i2p1));
               const double hDD_dD120 = invdx0*(FDPart1_Rational_1_60*(-hDD12_i0m3_i1_i2 + hDD12_i0p3_i1_i2) + FDPart1_Rational_3_20*(hDD12_i0m2_i1_i2 - hDD12_i0p2_i1_i2) + FDPart1_Rational_3_4*(-hDD12_i0m1_i1_i2 + hDD12_i0p1_i1_i2));
               const double hDD_dD121 = invdx1*(FDPart1_Rational_1_60*(-hDD12_i0_i1m3_i2 + hDD12_i0_i1p3_i2) + FDPart1_Rational_3_20*(hDD12_i0_i1m2_i2 - hDD12_i0_i1p2_i2) + FDPart1_Rational_3_4*(-hDD12_i0_i1m1_i2 + hDD12_i0_i1p1_i2));
               const double hDD_dD122 = invdx2*(FDPart1_Rational_1_60*(-hDD12_i0_i1_i2m3 + hDD12_i0_i1_i2p3) + FDPart1_Rational_3_20*(hDD12_i0_i1_i2m2 - hDD12_i0_i1_i2p2) + FDPart1_Rational_3_4*(-hDD12_i0_i1_i2m1 + hDD12_i0_i1_i2p1));
               const double hDD_dD220 = invdx0*(FDPart1_Rational_1_60*(-hDD22_i0m3_i1_i2 + hDD22_i0p3_i1_i2) + FDPart1_Rational_3_20*(hDD22_i0m2_i1_i2 - hDD22_i0p2_i1_i2) + FDPart1_Rational_3_4*(-hDD22_i0m1_i1_i2 + hDD22_i0p1_i1_i2));
               const double hDD_dD221 = invdx1*(FDPart1_Rational_1_60*(-hDD22_i0_i1m3_i2 + hDD22_i0_i1p3_i2) + FDPart1_Rational_3_20*(hDD22_i0_i1m2_i2 - hDD22_i0_i1p2_i2) + FDPart1_Rational_3_4*(-hDD22_i0_i1m1_i2 + hDD22_i0_i1p1_i2));
               const double hDD_dD222 = invdx2*(FDPart1_Rational_1_60*(-hDD22_i0_i1_i2m3 + hDD22_i0_i1_i2p3) + FDPart1_Rational_3_20*(hDD22_i0_i1_i2m2 - hDD22_i0_i1_i2p2) + FDPart1_Rational_3_4*(-hDD22_i0_i1_i2m1 + hDD22_i0_i1_i2p1));
               /*
                * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                */
               const double FDPart3_0 = hDD22 + 1;
               const double FDPart3_1 = -FDPart3_0*hDD01 + hDD02*hDD12;
               const double FDPart3_5 = hDD11 + 1;
               const double FDPart3_7 = hDD00 + 1;
               const double FDPart3_9 = (1.0/(FDPart3_0*FDPart3_5*FDPart3_7 - FDPart3_0*((hDD01)*(hDD01)) - FDPart3_5*((hDD02)*(hDD02)) - FDPart3_7*((hDD12)*(hDD12)) + 2*hDD01*hDD02*hDD12));
               const double FDPart3_10 = (1.0/2.0)*FDPart3_9;
               const double FDPart3_11 = FDPart3_1*FDPart3_10;
               const double FDPart3_12 = -FDPart3_5*hDD02 + hDD01*hDD12;
               const double FDPart3_13 = FDPart3_10*FDPart3_12;
               const double FDPart3_14 = hDD_dD012 + hDD_dD021 - hDD_dD120;
               const double FDPart3_15 = FDPart3_9*(FDPart3_0*FDPart3_5 - ((hDD12)*(hDD12)));
               const double FDPart3_16 = (1.0/2.0)*FDPart3_15;
               const double FDPart3_17 = -FDPart3_7*hDD12 + hDD01*hDD02;
               const double FDPart3_18 = 2*FDPart3_9;
               const double FDPart3_19 = FDPart3_17*FDPart3_18;
               const double FDPart3_20 = hDD_dD012 - hDD_dD021 + hDD_dD120;
               const double FDPart3_21 = FDPart3_12*FDPart3_18;
               const double FDPart3_22 = -hDD_dD012 + hDD_dD021 + hDD_dD120;
               const double FDPart3_23 = FDPart3_1*FDPart3_18;
               const double FDPart3_24 = 2*hDD_dD122 - hDD_dD221;
               const double FDPart3_25 = 2*hDD_dD022 - hDD_dD220;
               const double FDPart3_26 = FDPart3_9*(FDPart3_5*FDPart3_7 - ((hDD01)*(hDD01)));
               const double FDPart3_27 = -hDD_dD112 + 2*hDD_dD121;
               const double FDPart3_28 = 2*hDD_dD011 - hDD_dD110;
               const double FDPart3_29 = FDPart3_9*(FDPart3_0*FDPart3_7 - ((hDD02)*(hDD02)));
               const double FDPart3_30 = -hDD_dD001 + 2*hDD_dD010;
               const double FDPart3_31 = -hDD_dD002 + 2*hDD_dD020;
               const double FDPart3_32 = FDPart3_10*FDPart3_17;
               const double FDPart3_33 = (1.0/2.0)*FDPart3_29;
               const double FDPart3_34 = (1.0/2.0)*FDPart3_26;
               lambdaU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_15*(FDPart3_11*FDPart3_30 + FDPart3_13*FDPart3_31 + FDPart3_16*hDD_dD000) + FDPart3_19*(FDPart3_11*hDD_dD112 + FDPart3_13*hDD_dD221 + FDPart3_14*FDPart3_16) + FDPart3_21*(FDPart3_11*FDPart3_20 + FDPart3_13*hDD_dD220 + FDPart3_16*hDD_dD002) + FDPart3_23*(FDPart3_11*hDD_dD110 + FDPart3_13*FDPart3_22 + FDPart3_16*hDD_dD001) + FDPart3_26*(FDPart3_11*FDPart3_24 + FDPart3_13*hDD_dD222 + FDPart3_16*FDPart3_25) + FDPart3_29*(FDPart3_11*hDD_dD111 + FDPart3_13*FDPart3_27 + FDPart3_16*FDPart3_28);
               lambdaU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_15*(FDPart3_11*hDD_dD000 + FDPart3_30*FDPart3_33 + FDPart3_31*FDPart3_32) + FDPart3_19*(FDPart3_11*FDPart3_14 + FDPart3_32*hDD_dD221 + FDPart3_33*hDD_dD112) + FDPart3_21*(FDPart3_11*hDD_dD002 + FDPart3_20*FDPart3_33 + FDPart3_32*hDD_dD220) + FDPart3_23*(FDPart3_11*hDD_dD001 + FDPart3_22*FDPart3_32 + FDPart3_33*hDD_dD110) + FDPart3_26*(FDPart3_11*FDPart3_25 + FDPart3_24*FDPart3_33 + FDPart3_32*hDD_dD222) + FDPart3_29*(FDPart3_11*FDPart3_28 + FDPart3_27*FDPart3_32 + FDPart3_33*hDD_dD111);
               lambdaU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_15*(FDPart3_13*hDD_dD000 + FDPart3_30*FDPart3_32 + FDPart3_31*FDPart3_34) + FDPart3_19*(FDPart3_13*FDPart3_14 + FDPart3_32*hDD_dD112 + FDPart3_34*hDD_dD221) + FDPart3_21*(FDPart3_13*hDD_dD002 + FDPart3_20*FDPart3_32 + FDPart3_34*hDD_dD220) + FDPart3_23*(FDPart3_13*hDD_dD001 + FDPart3_22*FDPart3_34 + FDPart3_32*hDD_dD110) + FDPart3_26*(FDPart3_13*FDPart3_25 + FDPart3_24*FDPart3_32 + FDPart3_34*hDD_dD222) + FDPart3_29*(FDPart3_13*FDPart3_28 + FDPart3_27*FDPart3_34 + FDPart3_32*hDD_dD111);
            
        } // END LOOP: for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++)
    } // END LOOP: for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++)
} // END LOOP: for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++)
