//
// Tuple types that enable fast evaluation on GPU
//
#ifndef TUPLETYPESCUDA_H
#define TUPLETYPESCUDA_H

#ifdef NAMD_CUDA
#include <cuda_runtime.h>  // float3

struct CudaBond {
  int i, j, itype;
  // int ivir;
  float scale;
  float3 ioffsetXYZ;
};

struct CudaAngle {
  int i, j, k, itype;
  // int ivir, kvir;
  float scale;
  float3 ioffsetXYZ;
  float3 koffsetXYZ;
};

struct CudaDihedral {
  int i, j, k, l, itype;
  // int ivir, jvir, lvir;
  float scale;
  float3 ioffsetXYZ;
  float3 joffsetXYZ;
  float3 loffsetXYZ;
};

struct CudaExclusion {
  int i, j, vdwtypei, vdwtypej;
  // int ivir;
  float3 ioffsetXYZ;
};

struct CudaCrossterm {
  int i1, i2, i3, i4, i5, i6, i7, i8, itype;
  float scale;
  float3 offset12XYZ;
  float3 offset23XYZ;
  float3 offset34XYZ;
  float3 offset56XYZ;
  float3 offset67XYZ;
  float3 offset78XYZ;
};

struct CudaBondValue {
  float k;   //  Force constant for the bond
  float x0;  //  Rest distance for the bond
  float x1;  //  Upper wall for harmonic wall potential (with x0 lower wall)
};

struct CudaAngleValue {
  float k;   //  Force constant for angle
  float theta0;  //  Rest angle for angle
  float k_ub;  //  Urey-Bradley force constant
  float r_ub;  //  Urey-Bradley distance
  int normal; // Whether we use harmonic (0) or cos-based (1) angle terms
};

struct CudaDihedralValue {
  float k;     //  Force constant
  float delta; //  Phase shift
  int n;       //  Periodicity*2, if n low bit is set to 0, this is the last in multiplicity
};

// struct CudaCrosstermData { float d00,d01,d10,d11; };

struct CudaCrosstermValue {
  enum {dim=24};
  float4 c[dim][dim][4]; // bicubic interpolation coefficients
};

#endif

#endif // TUPLETYPESCUDA_H
