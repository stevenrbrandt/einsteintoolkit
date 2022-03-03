static const char* const differencing =
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth1(u) (kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))\n"
"#else\n"
"#  define PDstandardNth1(u) (PDstandardNth1_impl(u,p1o12dx,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth1_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dx, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return kmul(p1o12dx,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,-1,0,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,1,0,0),ksub(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth2(u) (kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))\n"
"#else\n"
"#  define PDstandardNth2(u) (PDstandardNth2_impl(u,p1o12dy,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth2_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dy, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return kmul(p1o12dy,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,-1,0),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,1,0),ksub(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth3(u) (kmul(p1o12dz,kmadd(ToReal(-8),KRANC_GFOFFSET3D(u,0,0,-1),kmadd(ToReal(8),KRANC_GFOFFSET3D(u,0,0,1),ksub(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))\n"
"#else\n"
"#  define PDstandardNth3(u) (PDstandardNth3_impl(u,p1o12dz,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth3_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o12dz, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return PDstandardNth2_impl(u, p1o12dz, cdk, cdj);\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth11(u) (kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0))))))\n"
"#else\n"
"#  define PDstandardNth11(u) (PDstandardNth11_impl(u,pm1o12dx2,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth11_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dx2, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return kmul(pm1o12dx2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,-1,0,0),KRANC_GFOFFSET3D(u,1,0,0)),kadd(KRANC_GFOFFSET3D(u,-2,0,0),KRANC_GFOFFSET3D(u,2,0,0)))));\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth22(u) (kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0))))))\n"
"#else\n"
"#  define PDstandardNth22(u) (PDstandardNth22_impl(u,pm1o12dy2,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth22_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dy2, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return kmul(pm1o12dy2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,-1,0),KRANC_GFOFFSET3D(u,0,1,0)),kadd(KRANC_GFOFFSET3D(u,0,-2,0),KRANC_GFOFFSET3D(u,0,2,0)))));\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth33(u) (kmul(pm1o12dz2,kmadd(ToReal(30),KRANC_GFOFFSET3D(u,0,0,0),kmadd(ToReal(-16),kadd(KRANC_GFOFFSET3D(u,0,0,-1),KRANC_GFOFFSET3D(u,0,0,1)),kadd(KRANC_GFOFFSET3D(u,0,0,-2),KRANC_GFOFFSET3D(u,0,0,2))))))\n"
"#else\n"
"#  define PDstandardNth33(u) (PDstandardNth33_impl(u,pm1o12dz2,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth33_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC pm1o12dz2, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return PDstandardNth22_impl(u, pm1o12dz2, cdk, cdj);\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth12(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))\n"
"#else\n"
"#  define PDstandardNth12(u) (PDstandardNth12_impl(u,p1o144dxdy,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth12_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0))))))));\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth13(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))\n"
"#else\n"
"#  define PDstandardNth13(u) (PDstandardNth13_impl(u,p1o144dxdz,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth13_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return PDstandardNth12_impl(u, p1o144dxdz, cdk, cdj);\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth21(u) (kmul(p1o144dxdy,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,1,0),KRANC_GFOFFSET3D(u,1,-1,0)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,-1,0),KRANC_GFOFFSET3D(u,1,1,0)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,2,0),kadd(KRANC_GFOFFSET3D(u,1,-2,0),kadd(KRANC_GFOFFSET3D(u,-2,1,0),KRANC_GFOFFSET3D(u,2,-1,0)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,-2,0),kadd(KRANC_GFOFFSET3D(u,1,2,0),kadd(KRANC_GFOFFSET3D(u,-2,-1,0),KRANC_GFOFFSET3D(u,2,1,0)))),kadd(KRANC_GFOFFSET3D(u,-2,-2,0),ksub(ksub(KRANC_GFOFFSET3D(u,2,2,0),KRANC_GFOFFSET3D(u,2,-2,0)),KRANC_GFOFFSET3D(u,-2,2,0)))))))))\n"
"#else\n"
"#  define PDstandardNth21(u) (PDstandardNth21_impl(u,p1o144dxdy,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth21_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdy, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return PDstandardNth12_impl(u, p1o144dxdy, cdj, cdk);\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth23(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))\n"
"#else\n"
"#  define PDstandardNth23(u) (PDstandardNth23_impl(u,p1o144dydz,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth23_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2))))))));\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth31(u) (kmul(p1o144dxdz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,-1,0,1),KRANC_GFOFFSET3D(u,1,0,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,-1,0,-1),KRANC_GFOFFSET3D(u,1,0,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,-1,0,2),kadd(KRANC_GFOFFSET3D(u,1,0,-2),kadd(KRANC_GFOFFSET3D(u,-2,0,1),KRANC_GFOFFSET3D(u,2,0,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,-1,0,-2),kadd(KRANC_GFOFFSET3D(u,1,0,2),kadd(KRANC_GFOFFSET3D(u,-2,0,-1),KRANC_GFOFFSET3D(u,2,0,1)))),kadd(KRANC_GFOFFSET3D(u,-2,0,-2),ksub(ksub(KRANC_GFOFFSET3D(u,2,0,2),KRANC_GFOFFSET3D(u,2,0,-2)),KRANC_GFOFFSET3D(u,-2,0,2)))))))))\n"
"#else\n"
"#  define PDstandardNth31(u) (PDstandardNth31_impl(u,p1o144dxdz,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth31_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dxdz, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return PDstandardNth12_impl(u, p1o144dxdz, cdk, cdj);\n"
"}\n"
"#endif\n"
"\n"
"#ifndef KRANC_DIFF_FUNCTIONS\n"
"#  define PDstandardNth32(u) (kmul(p1o144dydz,kmadd(ToReal(-64),kadd(KRANC_GFOFFSET3D(u,0,-1,1),KRANC_GFOFFSET3D(u,0,1,-1)),kmadd(ToReal(64),kadd(KRANC_GFOFFSET3D(u,0,-1,-1),KRANC_GFOFFSET3D(u,0,1,1)),kmadd(ToReal(8),kadd(KRANC_GFOFFSET3D(u,0,-1,2),kadd(KRANC_GFOFFSET3D(u,0,1,-2),kadd(KRANC_GFOFFSET3D(u,0,-2,1),KRANC_GFOFFSET3D(u,0,2,-1)))),kmadd(ToReal(-8),kadd(KRANC_GFOFFSET3D(u,0,-1,-2),kadd(KRANC_GFOFFSET3D(u,0,1,2),kadd(KRANC_GFOFFSET3D(u,0,-2,-1),KRANC_GFOFFSET3D(u,0,2,1)))),kadd(KRANC_GFOFFSET3D(u,0,-2,-2),ksub(ksub(KRANC_GFOFFSET3D(u,0,2,2),KRANC_GFOFFSET3D(u,0,2,-2)),KRANC_GFOFFSET3D(u,0,-2,2)))))))))\n"
"#else\n"
"#  define PDstandardNth32(u) (PDstandardNth32_impl(u,p1o144dydz,cdj,cdk))\n"
"static CCTK_REAL_VEC PDstandardNth32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk) CCTK_ATTRIBUTE_NOINLINE CCTK_ATTRIBUTE_UNUSED;\n"
"static CCTK_REAL_VEC PDstandardNth32_impl(const CCTK_REAL* restrict const u, const CCTK_REAL_VEC p1o144dydz, const ptrdiff_t cdj, const ptrdiff_t cdk)\n"
"{\n"
"  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL);\n"
"  return PDstandardNth23_impl(u, p1o144dydz, cdj, cdk);\n"
"}\n"
"#endif\n"
"\n"
""
;