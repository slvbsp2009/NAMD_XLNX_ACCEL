#ifdef WIN32
#define _USE_MATH_DEFINES
#define __thread __declspec(thread)
#endif
#include <math.h>
#include <cuda.h>
#include <cub/cub.cuh>
#include "ComputeBondedCUDAKernel.h"

#ifdef FUTURE_CUDADEVICE
#include "CudaDevice.h"
#else
#include "DeviceCUDA.h"
extern __thread DeviceCUDA *deviceCUDA;
#endif

__device__ inline long long int lliroundf(float f)
{
    long long int l;
    asm("cvt.rni.s64.f32  %0, %1;" : "=l"(l) : "f"(f));
    return l;
}

__device__ inline unsigned long long int llitoulli(long long int l)
{
    unsigned long long int u;
    asm("mov.b64    %0, %1;" : "=l"(u) : "l"(l));
    return u;
}

template <typename T>
__forceinline__ __device__
void convertForces(const float x, const float y, const float z,
  T &Tx, T &Ty, T &Tz) {

  Tx = (T)(x);
  Ty = (T)(y);
  Tz = (T)(z);
}

template <>
__forceinline__ __device__
void convertForces<long long int>(const float x, const float y, const float z,
  long long int &Tx, long long int &Ty, long long int &Tz) {

  Tx = lliroundf(x*float_to_force);
  Ty = lliroundf(y*float_to_force);
  Tz = lliroundf(z*float_to_force);
}

template <typename T>
__forceinline__ __device__
void calcComponentForces(
  float fij,
  const float dx, const float dy, const float dz,
  T &fxij, T &fyij, T &fzij) {

  fxij = (T)(fij*dx);
  fyij = (T)(fij*dy);
  fzij = (T)(fij*dz);

}

template <>
__forceinline__ __device__
void calcComponentForces<long long int>(
  float fij,
  const float dx, const float dy, const float dz,
  long long int &fxij,
  long long int &fyij,
  long long int &fzij) {

  fij *= float_to_force;
  fxij = lliroundf(fij*dx);
  fyij = lliroundf(fij*dy);
  fzij = lliroundf(fij*dz);

}

template <typename T>
__forceinline__ __device__
void calcComponentForces(
  float fij1,
  const float dx1, const float dy1, const float dz1,
  float fij2,
  const float dx2, const float dy2, const float dz2,
  T &fxij, T &fyij, T &fzij) {

  fxij = (T)(fij1*dx1 + fij2*dx2);
  fyij = (T)(fij1*dy1 + fij2*dy2);
  fzij = (T)(fij1*dz1 + fij2*dz2);
}

template <>
__forceinline__ __device__
void calcComponentForces<long long int>(
  float fij1,
  const float dx1,
  const float dy1,
  const float dz1,
  float fij2,
  const float dx2,
  const float dy2,
  const float dz2,
  long long int &fxij,
  long long int &fyij,
  long long int &fzij) {

  fij1 *= float_to_force;
  fij2 *= float_to_force;
  fxij = lliroundf(fij1*dx1 + fij2*dx2);
  fyij = lliroundf(fij1*dy1 + fij2*dy2);
  fzij = lliroundf(fij1*dz1 + fij2*dz2);
}

template <typename T>
__forceinline__ __device__
void storeForces(const T fx, const T fy, const T fz,
     const int ind, const int stride,
     T* force) {
  // The generic version can not be used
}

// Template specialization for 64bit integer = "long long int"
template <>
__forceinline__ __device__ 
void storeForces <long long int> (const long long int fx,
          const long long int fy,
          const long long int fz,
          const int ind, const int stride,
          long long int* force) {
  atomicAdd((unsigned long long int *)&force[ind           ], llitoulli(fx));
  atomicAdd((unsigned long long int *)&force[ind + stride  ], llitoulli(fy));
  atomicAdd((unsigned long long int *)&force[ind + stride*2], llitoulli(fz));
}

//
// Calculates bonds
//
template <typename T, bool doEnergy, bool doVirial>
__device__ void bondForce(
  const int index,
  const CudaBond* __restrict__ bondList,
  const CudaBondValue* __restrict__ bondValues,
  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  T* __restrict__ force, double &energy,
#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double>& virial
#else
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virial
#endif
  ) {

  CudaBond bl = bondList[index];
  CudaBondValue bondValue = bondValues[bl.itype];
  if (bondValue.x0 == 0.0f) return;

  float shx = bl.ioffsetXYZ.x*lata.x + bl.ioffsetXYZ.y*latb.x + bl.ioffsetXYZ.z*latc.x;
  float shy = bl.ioffsetXYZ.x*lata.y + bl.ioffsetXYZ.y*latb.y + bl.ioffsetXYZ.z*latc.y;
  float shz = bl.ioffsetXYZ.x*lata.z + bl.ioffsetXYZ.y*latb.z + bl.ioffsetXYZ.z*latc.z;

  float4 xyzqi = xyzq[bl.i];
  float4 xyzqj = xyzq[bl.j];

  float xij = xyzqi.x + shx - xyzqj.x;
  float yij = xyzqi.y + shy - xyzqj.y;
  float zij = xyzqi.z + shz - xyzqj.z;

  float r = sqrtf(xij*xij + yij*yij + zij*zij);

  float db = r - bondValue.x0;
  if (bondValue.x1) {
    // in this case, the bond represents a harmonic wall potential
    // where x0 is the lower wall and x1 is the upper
    db = (r > bondValue.x1 ? r - bondValue.x1 : (r > bondValue.x0 ? 0 : db));
  }
  float fij = db * bondValue.k * bl.scale;
 
  if (doEnergy) {
    energy += (double)(fij*db);
  }
  fij *= -2.0f/r;
  
  T T_fxij, T_fyij, T_fzij;
  calcComponentForces<T>(fij, xij, yij, zij, T_fxij, T_fyij, T_fzij);
  
  // Store forces
  storeForces<T>(T_fxij, T_fyij, T_fzij, bl.i, stride, force);
  storeForces<T>(-T_fxij, -T_fyij, -T_fzij, bl.j, stride, force);

  // Store virial
  if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
    float fxij = fij*xij;
    float fyij = fij*yij;
    float fzij = fij*zij;
    virial.xx = (fxij*xij);
    virial.xy = (fxij*yij);
    virial.xz = (fxij*zij);
    virial.yx = (fyij*xij);
    virial.yy = (fyij*yij);
    virial.yz = (fyij*zij);
    virial.zx = (fzij*xij);
    virial.zy = (fzij*yij);
    virial.zz = (fzij*zij);
#endif
  }
}

//
// Calculates modified exclusions
//
template <typename T, bool doEnergy, bool doVirial, bool doElect>
__device__ void modifiedExclusionForce(
  const int index,
  const CudaExclusion* __restrict__ exclusionList,
  const bool doSlow,
  const float one_scale14,                // 1 - scale14
  const int vdwCoefTableWidth,
#if __CUDA_ARCH__ >= 350
  const float2* __restrict__ vdwCoefTable,
#else
  cudaTextureObject_t vdwCoefTableTex, 
#endif
  cudaTextureObject_t forceTableTex, cudaTextureObject_t energyTableTex,
  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  const float cutoff2,
  double &energyVdw,
  T* __restrict__ forceNbond, double &energyNbond,
  T* __restrict__ forceSlow, double &energySlow,
#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double>& virialNbond,
  ComputeBondedCUDAKernel::BondedVirial<double>& virialSlow
#else
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virialNbond,
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virialSlow
#endif
  ) {

  CudaExclusion bl = exclusionList[index];

  float shx = bl.ioffsetXYZ.x*lata.x + bl.ioffsetXYZ.y*latb.x + bl.ioffsetXYZ.z*latc.x;
  float shy = bl.ioffsetXYZ.x*lata.y + bl.ioffsetXYZ.y*latb.y + bl.ioffsetXYZ.z*latc.y;
  float shz = bl.ioffsetXYZ.x*lata.z + bl.ioffsetXYZ.y*latb.z + bl.ioffsetXYZ.z*latc.z;

  float4 xyzqi = xyzq[bl.i];
  float4 xyzqj = xyzq[bl.j];

  float xij = xyzqi.x + shx - xyzqj.x;
  float yij = xyzqi.y + shy - xyzqj.y;
  float zij = xyzqi.z + shz - xyzqj.z;

  float r2 = xij*xij + yij*yij + zij*zij;
  if (r2 < cutoff2) {

    float rinv = rsqrtf(r2);

    float qq;
    if (doElect) qq = one_scale14 * xyzqi.w * xyzqj.w;

    int vdwIndex = bl.vdwtypei + bl.vdwtypej*vdwCoefTableWidth;
#if __CUDA_ARCH__ >= 350
    float2 ljab = __ldg(&vdwCoefTable[vdwIndex]);
#else
    float2 ljab = tex1Dfetch<float2>(vdwCoefTableTex, vdwIndex);
#endif

    float4 fi = tex1D<float4>(forceTableTex, rinv);
    float4 ei;
    if (doEnergy) {
      ei = tex1D<float4>(energyTableTex, rinv);
      energyVdw   += (double)(ljab.x * ei.z + ljab.y * ei.y);
      if (doElect) {
        energyNbond += qq * ei.x;
        if (doSlow) energySlow  += qq * ei.w;
      }
    }

    float fNbond;
    if (doElect) {
      fNbond = -(ljab.x * fi.z + ljab.y * fi.y + qq * fi.x);
    } else {
      fNbond = -(ljab.x * fi.z + ljab.y * fi.y);
    }
    T T_fxij, T_fyij, T_fzij;
    calcComponentForces<T>(fNbond, xij, yij, zij, T_fxij, T_fyij, T_fzij);
    storeForces<T>(T_fxij, T_fyij, T_fzij, bl.i, stride, forceNbond);
    storeForces<T>(-T_fxij, -T_fyij, -T_fzij, bl.j, stride, forceNbond);

    float fSlow;
    if (doSlow && doElect) {
      fSlow = -qq * fi.w;
      T T_fxij, T_fyij, T_fzij;
      calcComponentForces<T>(fSlow, xij, yij, zij, T_fxij, T_fyij, T_fzij);
      storeForces<T>(T_fxij, T_fyij, T_fzij, bl.i, stride, forceSlow);
      storeForces<T>(-T_fxij, -T_fyij, -T_fzij, bl.j, stride, forceSlow);
    }

    // Store virial
    if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
      float fxij = fNbond*xij;
      float fyij = fNbond*yij;
      float fzij = fNbond*zij;
      virialNbond.xx = (fxij*xij);
      virialNbond.xy = (fxij*yij);
      virialNbond.xz = (fxij*zij);
      virialNbond.yx = (fyij*xij);
      virialNbond.yy = (fyij*yij);
      virialNbond.yz = (fyij*zij);
      virialNbond.zx = (fzij*xij);
      virialNbond.zy = (fzij*yij);
      virialNbond.zz = (fzij*zij);
#endif
    }

    // Store virial
    if (doVirial && doSlow && doElect) {
#ifdef WRITE_FULL_VIRIALS
      float fxij = fSlow*xij;
      float fyij = fSlow*yij;
      float fzij = fSlow*zij;
      virialSlow.xx = (fxij*xij);
      virialSlow.xy = (fxij*yij);
      virialSlow.xz = (fxij*zij);
      virialSlow.yx = (fyij*xij);
      virialSlow.yy = (fyij*yij);
      virialSlow.yz = (fyij*zij);
      virialSlow.zx = (fzij*xij);
      virialSlow.zy = (fzij*yij);
      virialSlow.zz = (fzij*zij);
#endif
    }

  }
}

//
// Calculates exclusions. Here doSlow = true
//
template <typename T, bool doEnergy, bool doVirial>
__device__ void exclusionForce(
  const int index,
  const CudaExclusion* __restrict__ exclusionList,
  const float r2_delta, const int r2_delta_expc,
#if __CUDA_ARCH__ >= 350
  const float* __restrict__ r2_table,
  const float4* __restrict__ exclusionTable,
#else
  cudaTextureObject_t r2_table_tex,
  cudaTextureObject_t exclusionTableTex,
#endif
  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  const float cutoff2,
  T* __restrict__ forceSlow, double &energySlow,
#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double>& virialSlow
#else
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virialSlow
#endif
  ) {

  CudaExclusion bl = exclusionList[index];

  float shx = bl.ioffsetXYZ.x*lata.x + bl.ioffsetXYZ.y*latb.x + bl.ioffsetXYZ.z*latc.x;
  float shy = bl.ioffsetXYZ.x*lata.y + bl.ioffsetXYZ.y*latb.y + bl.ioffsetXYZ.z*latc.y;
  float shz = bl.ioffsetXYZ.x*lata.z + bl.ioffsetXYZ.y*latb.z + bl.ioffsetXYZ.z*latc.z;

  float4 xyzqi = xyzq[bl.i];
  float4 xyzqj = xyzq[bl.j];

  float xij = xyzqi.x + shx - xyzqj.x;
  float yij = xyzqi.y + shy - xyzqj.y;
  float zij = xyzqi.z + shz - xyzqj.z;

  float r2 = xij*xij + yij*yij + zij*zij;
  if (r2 < cutoff2) {
    r2 += r2_delta;

    union { float f; int i; } r2i;
    r2i.f = r2;
    int table_i = (r2i.i >> 17) + r2_delta_expc;  // table_i >= 0

#if __CUDA_ARCH__ >= 350
    float r2_table_val = __ldg(&r2_table[table_i]);
#else
    float r2_table_val = tex1Dfetch<float>(r2_table_tex, table_i);
#endif
    float diffa = r2 - r2_table_val;
    float qq = xyzqi.w * xyzqj.w;

#if __CUDA_ARCH__ >= 350
    float4 slow = __ldg(&exclusionTable[table_i]);
#else
    float4 slow = tex1Dfetch<float4>(exclusionTableTex, table_i);
#endif

    if (doEnergy) {
      energySlow += (double)(qq*(((diffa * (1.0f/6.0f)*slow.x + 0.25f*slow.y ) * diffa + 0.5f*slow.z ) * diffa + slow.w));
    }

    float fSlow = -qq*((diffa * slow.x + slow.y) * diffa + slow.z);

    T T_fxij, T_fyij, T_fzij;
    calcComponentForces<T>(fSlow, xij, yij, zij, T_fxij, T_fyij, T_fzij);
    storeForces<T>(T_fxij, T_fyij, T_fzij, bl.i, stride, forceSlow);
    storeForces<T>(-T_fxij, -T_fyij, -T_fzij, bl.j, stride, forceSlow);

    // Store virial
    if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
      float fxij = fSlow*xij;
      float fyij = fSlow*yij;
      float fzij = fSlow*zij;
      virialSlow.xx = (fxij*xij);
      virialSlow.xy = (fxij*yij);
      virialSlow.xz = (fxij*zij);
      virialSlow.yx = (fyij*xij);
      virialSlow.yy = (fyij*yij);
      virialSlow.yz = (fyij*zij);
      virialSlow.zx = (fzij*xij);
      virialSlow.zy = (fzij*yij);
      virialSlow.zz = (fzij*zij);
#endif
    }
  }
}

template <typename T, bool doEnergy, bool doVirial>
__device__ void angleForce(const int index,
  const CudaAngle* __restrict__ angleList,
  const CudaAngleValue* __restrict__ angleValues,
  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  T* __restrict__ force, double &energy,
#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double>& virial
#else
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virial
#endif
  ) {

  CudaAngle al = angleList[index];

  float ishx = al.ioffsetXYZ.x*lata.x + al.ioffsetXYZ.y*latb.x + al.ioffsetXYZ.z*latc.x;
  float ishy = al.ioffsetXYZ.x*lata.y + al.ioffsetXYZ.y*latb.y + al.ioffsetXYZ.z*latc.y;
  float ishz = al.ioffsetXYZ.x*lata.z + al.ioffsetXYZ.y*latb.z + al.ioffsetXYZ.z*latc.z;

  float kshx = al.koffsetXYZ.x*lata.x + al.koffsetXYZ.y*latb.x + al.koffsetXYZ.z*latc.x;
  float kshy = al.koffsetXYZ.x*lata.y + al.koffsetXYZ.y*latb.y + al.koffsetXYZ.z*latc.y;
  float kshz = al.koffsetXYZ.x*lata.z + al.koffsetXYZ.y*latb.z + al.koffsetXYZ.z*latc.z;

  float xij = xyzq[al.i].x + ishx - xyzq[al.j].x;
  float yij = xyzq[al.i].y + ishy - xyzq[al.j].y;
  float zij = xyzq[al.i].z + ishz - xyzq[al.j].z;

  float xkj = xyzq[al.k].x + kshx - xyzq[al.j].x;
  float ykj = xyzq[al.k].y + kshy - xyzq[al.j].y;
  float zkj = xyzq[al.k].z + kshz - xyzq[al.j].z;

  float rij_inv = rsqrtf(xij*xij + yij*yij + zij*zij);
  float rkj_inv = rsqrtf(xkj*xkj + ykj*ykj + zkj*zkj);

  float xijr = xij*rij_inv;
  float yijr = yij*rij_inv;
  float zijr = zij*rij_inv;
  float xkjr = xkj*rkj_inv;
  float ykjr = ykj*rkj_inv;
  float zkjr = zkj*rkj_inv;
  float cos_theta = xijr*xkjr + yijr*ykjr + zijr*zkjr;

  CudaAngleValue angleValue = angleValues[al.itype];
  angleValue.k *= al.scale;

  float diff;
  if (angleValue.normal == 1) {
    // Restrict values of cos_theta to the interval [-0.999 ... 0.999]
    cos_theta = min(0.999f, max(-0.999f, cos_theta));
    float theta = acosf(cos_theta);
    diff = theta - angleValue.theta0;
  } else {
    diff = cos_theta - angleValue.theta0;
  }

  if (doEnergy) {
    energy += (double)(angleValue.k * diff * diff);
  }
  if (angleValue.normal == 1) {
    float inv_sin_theta = rsqrtf(1.0f - cos_theta*cos_theta);
    if (inv_sin_theta > 1.0e6) {
      diff = (diff < 0.0f) ? 1.0f : -1.0f;
    } else {
      diff *= -inv_sin_theta;
    }
  }
  diff *= -2.0f*angleValue.k;

  float dtxi = rij_inv*(xkjr - cos_theta*xijr);
  float dtxj = rkj_inv*(xijr - cos_theta*xkjr);
  float dtyi = rij_inv*(ykjr - cos_theta*yijr);
  float dtyj = rkj_inv*(yijr - cos_theta*ykjr);
  float dtzi = rij_inv*(zkjr - cos_theta*zijr);
  float dtzj = rkj_inv*(zijr - cos_theta*zkjr);

  T T_dtxi, T_dtyi, T_dtzi;
  T T_dtxj, T_dtyj, T_dtzj;
  calcComponentForces<T>(diff, dtxi, dtyi, dtzi, T_dtxi, T_dtyi, T_dtzi);
  calcComponentForces<T>(diff, dtxj, dtyj, dtzj, T_dtxj, T_dtyj, T_dtzj);
  T T_dtxk = -T_dtxi - T_dtxj;
  T T_dtyk = -T_dtyi - T_dtyj;
  T T_dtzk = -T_dtzi - T_dtzj;
  storeForces<T>(T_dtxk, T_dtyk, T_dtzk, al.j, stride, force);

  if (angleValue.k_ub) {
    float xik = xij - xkj;
    float yik = yij - ykj;
    float zik = zij - zkj;
    float rik_inv = rsqrtf(xik*xik + yik*yik + zik*zik);
    float rik = 1.0f/rik_inv;
    float diff_ub = rik - angleValue.r_ub;
    if (doEnergy) {
      energy += (double)(angleValue.k_ub * diff_ub * diff_ub);
    }
    diff_ub *= -2.0f*angleValue.k_ub*rik_inv;
    T T_dxik, T_dyik, T_dzik;
    calcComponentForces<T>(diff_ub, xik, yik, zik, T_dxik, T_dyik, T_dzik);
    T_dtxi += T_dxik;
    T_dtyi += T_dyik;
    T_dtzi += T_dzik;
    T_dtxj -= T_dxik;
    T_dtyj -= T_dyik;
    T_dtzj -= T_dzik;
  }

  storeForces<T>(T_dtxi, T_dtyi, T_dtzi, al.i, stride, force);
  storeForces<T>(T_dtxj, T_dtyj, T_dtzj, al.k, stride, force);

  // Store virial
  if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
    float fxi = ((float)T_dtxi)*force_to_float;
    float fyi = ((float)T_dtyi)*force_to_float;
    float fzi = ((float)T_dtzi)*force_to_float;
    float fxk = ((float)T_dtxj)*force_to_float;
    float fyk = ((float)T_dtyj)*force_to_float;
    float fzk = ((float)T_dtzj)*force_to_float;
    virial.xx = (fxi*xij) + (fxk*xkj);
    virial.xy = (fxi*yij) + (fxk*ykj);
    virial.xz = (fxi*zij) + (fxk*zkj);
    virial.yx = (fyi*xij) + (fyk*xkj);
    virial.yy = (fyi*yij) + (fyk*ykj);
    virial.yz = (fyi*zij) + (fyk*zkj);
    virial.zx = (fzi*xij) + (fzk*xkj);
    virial.zy = (fzi*yij) + (fzk*ykj);
    virial.zz = (fzi*zij) + (fzk*zkj);
#endif
  }
}

//
// Dihedral computation
//
// Out: df, e
//
template <bool doEnergy>
__forceinline__ __device__
void diheComp(const CudaDihedralValue* dihedralValues, int ic,
  const float sin_phi, const float cos_phi, const float scale, float& df, double& e) {

  df = 0.0f;
  if (doEnergy) e = 0.0;

  float phi = atan2f(sin_phi, cos_phi);

  bool lrep = true;
  while (lrep) {
    CudaDihedralValue dihedralValue = dihedralValues[ic];
    dihedralValue.k *= scale;

    // Last dihedral has n low bit set to 0
    lrep = (dihedralValue.n & 1);
    dihedralValue.n >>= 1;
    if (dihedralValue.n) {
      float nf = dihedralValue.n;
      float x = nf*phi - dihedralValue.delta;
      if (doEnergy) {
        float sin_x, cos_x;
        sincosf(x, &sin_x, &cos_x);
        e += (double)(dihedralValue.k*(1.0f + cos_x));
        df += (double)(nf*dihedralValue.k*sin_x);
      } else {
        float sin_x = sinf(x);
        df += (double)(nf*dihedralValue.k*sin_x);
      }
    } else {
      float diff = phi - dihedralValue.delta;
      if (diff < -(float)(M_PI)) diff += (float)(2.0*M_PI);
      if (diff >  (float)(M_PI)) diff -= (float)(2.0*M_PI);
      if (doEnergy) {
        e += (double)(dihedralValue.k*diff*diff);
      }
      df -= 2.0f*dihedralValue.k*diff;
    }
    ic++;
  }

}


template <typename T, bool doEnergy, bool doVirial>
__device__ void diheForce(const int index,
  const CudaDihedral* __restrict__ diheList,
  const CudaDihedralValue* __restrict__ dihedralValues,
  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  T* __restrict__ force, double &energy,
#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double>& virial
#else
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virial
#endif
  ) {

  CudaDihedral dl = diheList[index];

  float ishx = dl.ioffsetXYZ.x*lata.x + dl.ioffsetXYZ.y*latb.x + dl.ioffsetXYZ.z*latc.x;
  float ishy = dl.ioffsetXYZ.x*lata.y + dl.ioffsetXYZ.y*latb.y + dl.ioffsetXYZ.z*latc.y;
  float ishz = dl.ioffsetXYZ.x*lata.z + dl.ioffsetXYZ.y*latb.z + dl.ioffsetXYZ.z*latc.z;

  float jshx = dl.joffsetXYZ.x*lata.x + dl.joffsetXYZ.y*latb.x + dl.joffsetXYZ.z*latc.x;
  float jshy = dl.joffsetXYZ.x*lata.y + dl.joffsetXYZ.y*latb.y + dl.joffsetXYZ.z*latc.y;
  float jshz = dl.joffsetXYZ.x*lata.z + dl.joffsetXYZ.y*latb.z + dl.joffsetXYZ.z*latc.z;

  float lshx = dl.loffsetXYZ.x*lata.x + dl.loffsetXYZ.y*latb.x + dl.loffsetXYZ.z*latc.x;
  float lshy = dl.loffsetXYZ.x*lata.y + dl.loffsetXYZ.y*latb.y + dl.loffsetXYZ.z*latc.y;
  float lshz = dl.loffsetXYZ.x*lata.z + dl.loffsetXYZ.y*latb.z + dl.loffsetXYZ.z*latc.z;

  float xij = xyzq[dl.i].x + ishx - xyzq[dl.j].x;
  float yij = xyzq[dl.i].y + ishy - xyzq[dl.j].y;
  float zij = xyzq[dl.i].z + ishz - xyzq[dl.j].z;

  float xjk = xyzq[dl.j].x + jshx - xyzq[dl.k].x;
  float yjk = xyzq[dl.j].y + jshy - xyzq[dl.k].y;
  float zjk = xyzq[dl.j].z + jshz - xyzq[dl.k].z;

  float xlk = xyzq[dl.l].x + lshx - xyzq[dl.k].x;
  float ylk = xyzq[dl.l].y + lshy - xyzq[dl.k].y;
  float zlk = xyzq[dl.l].z + lshz - xyzq[dl.k].z;

  // A=F^G, B=H^G.
  float ax = yij*zjk - zij*yjk;
  float ay = zij*xjk - xij*zjk;
  float az = xij*yjk - yij*xjk;
  float bx = ylk*zjk - zlk*yjk;
  float by = zlk*xjk - xlk*zjk;
  float bz = xlk*yjk - ylk*xjk;

  float ra2 = ax*ax + ay*ay + az*az;
  float rb2 = bx*bx + by*by + bz*bz;
  float rg = sqrtf(xjk*xjk + yjk*yjk + zjk*zjk);

  //    if((ra2 <= rxmin2) .or. (rb2 <= rxmin2) .or. (rg <= rxmin)) then
  //          nlinear = nlinear + 1
  //       endif

  float rgr = 1.0f / rg;
  float ra2r = 1.0f / ra2;
  float rb2r = 1.0f / rb2;
  float rabr = sqrtf(ra2r*rb2r);

  float cos_phi = (ax*bx + ay*by + az*bz)*rabr;
  //
  // Note that sin(phi).G/|G|=B^A/(|A|.|B|)
  // which can be simplify to sin(phi)=|G|H.A/(|A|.|B|)
  float sin_phi = rg*rabr*(ax*xlk + ay*ylk + az*zlk);
  //
  //     Energy and derivative contributions.

  float df;
  double e;
  diheComp<doEnergy>(dihedralValues, dl.itype, sin_phi, cos_phi, dl.scale, df, e);

  if (doEnergy) energy += e;

  //
  //     Compute derivatives wrt catesian coordinates.
  //
  // GAA=dE/dphi.|G|/A^2, GBB=dE/dphi.|G|/B^2, FG=F.G, HG=H.G
  //  FGA=dE/dphi*F.G/(|G|A^2), HGB=dE/dphi*H.G/(|G|B^2)

  float fg = xij*xjk + yij*yjk + zij*zjk;
  float hg = xlk*xjk + ylk*yjk + zlk*zjk;
  ra2r *= df;
  rb2r *= df;
  float fga = fg*ra2r*rgr;
  float hgb = hg*rb2r*rgr;
  float gaa = ra2r*rg;
  float gbb = rb2r*rg;
  // DFi=dE/dFi, DGi=dE/dGi, DHi=dE/dHi.

  // Store forces
  T T_fix, T_fiy, T_fiz;
  calcComponentForces<T>(-gaa, ax, ay, az, T_fix, T_fiy, T_fiz);
  storeForces<T>(T_fix, T_fiy, T_fiz, dl.i, stride, force);

  T dgx, dgy, dgz;
  calcComponentForces<T>(fga, ax, ay, az, -hgb, bx, by, bz, dgx, dgy, dgz);
  T T_fjx = dgx - T_fix;
  T T_fjy = dgy - T_fiy;
  T T_fjz = dgz - T_fiz;
  storeForces<T>(T_fjx, T_fjy, T_fjz, dl.j, stride, force);

  T dhx, dhy, dhz;
  calcComponentForces<T>(gbb, bx, by, bz, dhx, dhy, dhz);
  T T_fkx = -dhx - dgx;
  T T_fky = -dhy - dgy;
  T T_fkz = -dhz - dgz;
  storeForces<T>(T_fkx, T_fky, T_fkz, dl.k, stride, force);
  storeForces<T>(dhx, dhy, dhz, dl.l, stride, force);

  // Store virial
  if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
    float fix = ((float)T_fix)*force_to_float;
    float fiy = ((float)T_fiy)*force_to_float;
    float fiz = ((float)T_fiz)*force_to_float;
    float fjx = ((float)dgx)*force_to_float;
    float fjy = ((float)dgy)*force_to_float;
    float fjz = ((float)dgz)*force_to_float;
    float flx = ((float)dhx)*force_to_float;
    float fly = ((float)dhy)*force_to_float;
    float flz = ((float)dhz)*force_to_float;
    virial.xx = (fix*xij) + (fjx*xjk) + (flx*xlk);
    virial.xy = (fix*yij) + (fjx*yjk) + (flx*ylk);
    virial.xz = (fix*zij) + (fjx*zjk) + (flx*zlk);
    virial.yx = (fiy*xij) + (fjy*xjk) + (fly*xlk);
    virial.yy = (fiy*yij) + (fjy*yjk) + (fly*ylk);
    virial.yz = (fiy*zij) + (fjy*zjk) + (fly*zlk);
    virial.zx = (fiz*xij) + (fjz*xjk) + (flz*xlk);
    virial.zy = (fiz*yij) + (fjz*yjk) + (flz*ylk);
    virial.zz = (fiz*zij) + (fjz*zjk) + (flz*zlk);
#endif
  }

}

__device__ __forceinline__ float3 cross(const float3 v1, const float3 v2) {
 return make_float3(
  v1.y*v2.z-v2.y*v1.z,
  v2.x*v1.z-v1.x*v2.z,
  v1.x*v2.y-v2.x*v1.y
  );
}

__device__ __forceinline__ float dot(const float3 v1, const float3 v2) {
  return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

//
// Calculates crossterms
//
template <typename T, bool doEnergy, bool doVirial>
__device__ void crosstermForce(
  const int index,
  const CudaCrossterm* __restrict__ crosstermList,
  const CudaCrosstermValue* __restrict__ crosstermValues,
  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  T* __restrict__ force, double &energy,
#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double>& virial
#else
  ComputeBondedCUDAKernel::BondedVirial* __restrict__ virial
#endif
  ) {

  CudaCrossterm cl = crosstermList[index];

  // ----------------------------------------------------------------------------
  // Angle between 1 - 2 - 3 - 4
  //
  // 1 - 2
  float3 sh12 = make_float3(
    cl.offset12XYZ.x*lata.x + cl.offset12XYZ.y*latb.x + cl.offset12XYZ.z*latc.x,
    cl.offset12XYZ.x*lata.y + cl.offset12XYZ.y*latb.y + cl.offset12XYZ.z*latc.y,
    cl.offset12XYZ.x*lata.z + cl.offset12XYZ.y*latb.z + cl.offset12XYZ.z*latc.z);

  float4 xyzq1 = xyzq[cl.i1];
  float4 xyzq2 = xyzq[cl.i2];

  float3 r12 = make_float3(
    xyzq1.x + sh12.x - xyzq2.x,
    xyzq1.y + sh12.y - xyzq2.y,
    xyzq1.z + sh12.z - xyzq2.z);

  // 2 - 3
  float3 sh23 = make_float3(
    cl.offset23XYZ.x*lata.x + cl.offset23XYZ.y*latb.x + cl.offset23XYZ.z*latc.x,
    cl.offset23XYZ.x*lata.y + cl.offset23XYZ.y*latb.y + cl.offset23XYZ.z*latc.y,
    cl.offset23XYZ.x*lata.z + cl.offset23XYZ.y*latb.z + cl.offset23XYZ.z*latc.z);

  float4 xyzq3 = xyzq[cl.i3];

  float3 r23 = make_float3(
    xyzq2.x + sh23.x - xyzq3.x,
    xyzq2.y + sh23.y - xyzq3.y,
    xyzq2.z + sh23.z - xyzq3.z);

  // 3 - 4
  float3 sh34 = make_float3(
    cl.offset34XYZ.x*lata.x + cl.offset34XYZ.y*latb.x + cl.offset34XYZ.z*latc.x,
    cl.offset34XYZ.x*lata.y + cl.offset34XYZ.y*latb.y + cl.offset34XYZ.z*latc.y,
    cl.offset34XYZ.x*lata.z + cl.offset34XYZ.y*latb.z + cl.offset34XYZ.z*latc.z);

  float4 xyzq4 = xyzq[cl.i4];

  float3 r34 = make_float3(
    xyzq3.x + sh34.x - xyzq4.x,
    xyzq3.y + sh34.y - xyzq4.y,
    xyzq3.z + sh34.z - xyzq4.z);

  // Calculate the cross products
  float3 A = cross(r12, r23);
  float3 B = cross(r23, r34);
  float3 C = cross(r23, A);

  // Calculate the inverse distances
  float inv_rA = rsqrtf(dot(A, A));
  float inv_rB = rsqrtf(dot(B, B));
  float inv_rC = rsqrtf(dot(C, C));

  //  Calculate the sin and cos
  float cos_phi = dot(A, B)*(inv_rA*inv_rB);
  float sin_phi = dot(C, B)*(inv_rC*inv_rB);

  float phi = -atan2f(sin_phi,cos_phi);

  // ----------------------------------------------------------------------------
  // Angle between 5 - 6 - 7 - 8
  //

  // 5 - 6
  float3 sh56 = make_float3(
    cl.offset56XYZ.x*lata.x + cl.offset56XYZ.y*latb.x + cl.offset56XYZ.z*latc.x,
    cl.offset56XYZ.x*lata.y + cl.offset56XYZ.y*latb.y + cl.offset56XYZ.z*latc.y,
    cl.offset56XYZ.x*lata.z + cl.offset56XYZ.y*latb.z + cl.offset56XYZ.z*latc.z);

  float4 xyzq5 = xyzq[cl.i5];
  float4 xyzq6 = xyzq[cl.i6];

  float3 r56 = make_float3(
    xyzq5.x + sh56.x - xyzq6.x,
    xyzq5.y + sh56.y - xyzq6.y,
    xyzq5.z + sh56.z - xyzq6.z);

  // 6 - 7
  float3 sh67 = make_float3(
    cl.offset67XYZ.x*lata.x + cl.offset67XYZ.y*latb.x + cl.offset67XYZ.z*latc.x,
    cl.offset67XYZ.x*lata.y + cl.offset67XYZ.y*latb.y + cl.offset67XYZ.z*latc.y,
    cl.offset67XYZ.x*lata.z + cl.offset67XYZ.y*latb.z + cl.offset67XYZ.z*latc.z);

  float4 xyzq7 = xyzq[cl.i7];

  float3 r67 = make_float3(
    xyzq6.x + sh67.x - xyzq7.x,
    xyzq6.y + sh67.y - xyzq7.y,
    xyzq6.z + sh67.z - xyzq7.z);

  // 7 - 8
  float3 sh78 = make_float3(
    cl.offset78XYZ.x*lata.x + cl.offset78XYZ.y*latb.x + cl.offset78XYZ.z*latc.x,
    cl.offset78XYZ.x*lata.y + cl.offset78XYZ.y*latb.y + cl.offset78XYZ.z*latc.y,
    cl.offset78XYZ.x*lata.z + cl.offset78XYZ.y*latb.z + cl.offset78XYZ.z*latc.z);

  float4 xyzq8 = xyzq[cl.i8];

  float3 r78 = make_float3(
    xyzq7.x + sh78.x - xyzq8.x,
    xyzq7.y + sh78.y - xyzq8.y,
    xyzq7.z + sh78.z - xyzq8.z);

  // Calculate the cross products
  float3 D = cross(r56, r67);
  float3 E = cross(r67, r78);
  float3 F = cross(r67, D);
  
  // Calculate the inverse distances
  float inv_rD = rsqrtf(dot(D, D));
  float inv_rE = rsqrtf(dot(E, E));
  float inv_rF = rsqrtf(dot(F, F));

  //  Calculate the sin and cos
  float cos_psi = dot(D, E)*(inv_rD*inv_rE);
  float sin_psi = dot(F, E)*(inv_rF*inv_rE);

  float psi = -atan2f(sin_psi,cos_psi);

  // ----------------------------------------------------------------------------
  // Calculate the energy

  const float inv_h = (float)( (CudaCrosstermValue::dim) / (2.0*M_PI) );

  // Shift angles
  phi = fmod(phi + (float)M_PI, 2.0f*(float)M_PI);
  psi = fmod(psi + (float)M_PI, 2.0f*(float)M_PI);

  // distance measured in grid points between angle and smallest value
  float phi_h = phi * inv_h;
  float psi_h = psi * inv_h;

  // find smallest numbered grid point in stencil
  int iphi = (int)floor(phi_h);
  int ipsi = (int)floor(psi_h);

  const CudaCrosstermValue& cp = crosstermValues[cl.itype];

  // Load coefficients
  float4 c[4];
  c[0] = cp.c[iphi][ipsi][0];
  c[1] = cp.c[iphi][ipsi][1];
  c[2] = cp.c[iphi][ipsi][2];
  c[3] = cp.c[iphi][ipsi][3];

  float dphi = phi_h - (float)iphi;
  float dpsi = psi_h - (float)ipsi;

  if (doEnergy) {
    float energyf =          c[3].x + dphi*( c[3].y + dphi*( c[3].z + dphi*c[3].w ) );
    energyf = energyf*dpsi + c[2].x + dphi*( c[2].y + dphi*( c[2].z + dphi*c[2].w ) );
    energyf = energyf*dpsi + c[1].x + dphi*( c[1].y + dphi*( c[1].z + dphi*c[1].w ) );
    energyf = energyf*dpsi + c[0].x + dphi*( c[0].y + dphi*( c[0].z + dphi*c[0].w ) );
    energy += energyf*cl.scale;
  }

  float dEdphi =         3.0f*(c[0].w + dpsi*( c[1].w + dpsi*( c[2].w + dpsi*c[3].w ) ));
  dEdphi = dEdphi*dphi + 2.0f*(c[0].z + dpsi*( c[1].z + dpsi*( c[2].z + dpsi*c[3].z ) ));
  dEdphi = dEdphi*dphi +      (c[0].y + dpsi*( c[1].y + dpsi*( c[2].y + dpsi*c[3].y ) ));  // 13 muls
  dEdphi *= cl.scale*inv_h;

  float dEdpsi =         3.0f*(c[3].x + dphi*( c[3].y + dphi*( c[3].z + dphi*c[3].w ) ));
  dEdpsi = dEdpsi*dpsi + 2.0f*(c[2].x + dphi*( c[2].y + dphi*( c[2].z + dphi*c[2].w ) ));
  dEdpsi = dEdpsi*dpsi +      (c[1].x + dphi*( c[1].y + dphi*( c[1].z + dphi*c[1].w ) ));  // 13 muls
  dEdpsi *= cl.scale*inv_h;

  // float normCross1 = dot(A, A);
  float square_r23 = dot(r23, r23);
  float norm_r23 = sqrtf(square_r23);
  float inv_square_r23 = 1.0f/square_r23;
  float ff1 = dEdphi*norm_r23*inv_rA*inv_rA;
  float ff2 = -dot(r12, r23)*inv_square_r23;
  float ff3 = -dot(r34, r23)*inv_square_r23;
  float ff4 = -dEdphi*norm_r23*inv_rB*inv_rB;

  float3 f1 = make_float3(ff1*A.x, ff1*A.y, ff1*A.z);
  float3 f4 = make_float3(ff4*B.x, ff4*B.y, ff4*B.z);
  float3 t1 = make_float3( ff2*f1.x - ff3*f4.x, ff2*f1.y - ff3*f4.y, ff2*f1.z - ff3*f4.z );
  float3 f2 = make_float3(  t1.x - f1.x,  t1.y - f1.y,  t1.z - f1.z);
  float3 f3 = make_float3( -t1.x - f4.x, -t1.y - f4.y, -t1.z - f4.z);

  T T_f1x, T_f1y, T_f1z;
  T T_f2x, T_f2y, T_f2z;
  T T_f3x, T_f3y, T_f3z;
  T T_f4x, T_f4y, T_f4z;
  convertForces<T>(f1.x, f1.y, f1.z, T_f1x, T_f1y, T_f1z);
  convertForces<T>(f2.x, f2.y, f2.z, T_f2x, T_f2y, T_f2z);
  convertForces<T>(f3.x, f3.y, f3.z, T_f3x, T_f3y, T_f3z);
  convertForces<T>(f4.x, f4.y, f4.z, T_f4x, T_f4y, T_f4z);
  storeForces<T>(T_f1x, T_f1y, T_f1z, cl.i1, stride, force);
  storeForces<T>(T_f2x, T_f2y, T_f2z, cl.i2, stride, force);
  storeForces<T>(T_f3x, T_f3y, T_f3z, cl.i3, stride, force);
  storeForces<T>(T_f4x, T_f4y, T_f4z, cl.i4, stride, force);

  float square_r67 = dot(r67, r67);
  float norm_r67 = sqrtf(square_r67);
  float inv_square_r67 = 1.0f/(square_r67);
  ff1 = dEdpsi*norm_r67*inv_rD*inv_rD;
  ff2 = -dot(r56, r67)*inv_square_r67;
  ff3 = -dot(r78, r67)*inv_square_r67;
  ff4 = -dEdpsi*norm_r67*inv_rE*inv_rE;
  
  float3 f5 = make_float3( ff1*D.x, ff1*D.y, ff1*D.z );
  float3 f8 = make_float3( ff4*E.x, ff4*E.y, ff4*E.z );
  float3 t2 = make_float3( ff2*f5.x - ff3*f8.x, ff2*f5.y - ff3*f8.y, ff2*f5.z - ff3*f8.z );
  float3 f6 = make_float3( t2.x - f5.x,  t2.y - f5.y,  t2.z - f5.z );
  float3 f7 = make_float3(-t2.x - f8.x, -t2.y - f8.y, -t2.z - f8.z );

  T T_f5x, T_f5y, T_f5z;
  T T_f6x, T_f6y, T_f6z;
  T T_f7x, T_f7y, T_f7z;
  T T_f8x, T_f8y, T_f8z;
  convertForces<T>(f5.x, f5.y, f5.z, T_f5x, T_f5y, T_f5z);
  convertForces<T>(f6.x, f6.y, f6.z, T_f6x, T_f6y, T_f6z);
  convertForces<T>(f7.x, f7.y, f7.z, T_f7x, T_f7y, T_f7z);
  convertForces<T>(f8.x, f8.y, f8.z, T_f8x, T_f8y, T_f8z);
  storeForces<T>(T_f5x, T_f5y, T_f5z, cl.i5, stride, force);
  storeForces<T>(T_f6x, T_f6y, T_f6z, cl.i6, stride, force);
  storeForces<T>(T_f7x, T_f7y, T_f7z, cl.i7, stride, force);
  storeForces<T>(T_f8x, T_f8y, T_f8z, cl.i8, stride, force);

  // Store virial
  if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
    float3 s12 = make_float3( f1.x + f2.x, f1.y + f2.y, f1.z + f2.z );
    float3 s56 = make_float3( f5.x + f6.x, f5.y + f6.y, f5.z + f6.z );
    virial.xx = f1.x*r12.x + s12.x*r23.x - f4.x*r34.x + f5.x*r56.x + s56.x*r67.x - f8.x*r78.x;
    virial.xy = f1.x*r12.y + s12.x*r23.y - f4.x*r34.y + f5.x*r56.y + s56.x*r67.y - f8.x*r78.y;
    virial.xz = f1.x*r12.z + s12.x*r23.z - f4.x*r34.z + f5.x*r56.z + s56.x*r67.z - f8.x*r78.z;
    virial.yx = f1.y*r12.x + s12.y*r23.x - f4.y*r34.x + f5.y*r56.x + s56.y*r67.x - f8.y*r78.x;
    virial.yy = f1.y*r12.y + s12.y*r23.y - f4.y*r34.y + f5.y*r56.y + s56.y*r67.y - f8.y*r78.y;
    virial.yz = f1.y*r12.z + s12.y*r23.z - f4.y*r34.z + f5.y*r56.z + s56.y*r67.z - f8.y*r78.z;
    virial.zx = f1.z*r12.x + s12.z*r23.x - f4.z*r34.x + f5.z*r56.x + s56.z*r67.x - f8.z*r78.x;
    virial.zy = f1.z*r12.y + s12.z*r23.y - f4.z*r34.y + f5.z*r56.y + s56.z*r67.y - f8.z*r78.y;
    virial.zz = f1.z*r12.z + s12.z*r23.z - f4.z*r34.z + f5.z*r56.z + s56.z*r67.z - f8.z*r78.z;
  }
#endif

}

#define BONDEDFORCESKERNEL_NUM_WARP 4
//
// Calculates all forces in a single kernel call
//
template <typename T, bool doEnergy, bool doVirial>
__global__ void bondedForcesKernel(
  const int start,

  const int numBonds,
  const CudaBond* __restrict__ bonds,
  const CudaBondValue* __restrict__ bondValues,

  const int numAngles,
  const CudaAngle* __restrict__ angles,
  const CudaAngleValue* __restrict__ angleValues,

  const int numDihedrals,
  const CudaDihedral* __restrict__ dihedrals,
  const CudaDihedralValue* __restrict__ dihedralValues,

  const int numImpropers,
  const CudaDihedral* __restrict__ impropers,
  const CudaDihedralValue* __restrict__ improperValues,

  const int numExclusions,
  const CudaExclusion* __restrict__ exclusions,

  const int numCrossterms,
  const CudaCrossterm* __restrict__ crossterms,
  const CudaCrosstermValue* __restrict__ crosstermValues,

  const float cutoff2,
  const float r2_delta, const int r2_delta_expc,

  const float* __restrict__ r2_table,
  const float4* __restrict__ exclusionTable,
  cudaTextureObject_t r2_table_tex,
  cudaTextureObject_t exclusionTableTex,

  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  T* __restrict__ force,
  T* __restrict__ forceSlow,
  double* __restrict__ energies_virials) {

  // Thread-block index
  int indexTB = start + blockIdx.x;

  const int numBondsTB     = (numBonds + blockDim.x - 1)/blockDim.x;
  const int numAnglesTB    = (numAngles + blockDim.x - 1)/blockDim.x;
  const int numDihedralsTB = (numDihedrals + blockDim.x - 1)/blockDim.x;
  const int numImpropersTB = (numImpropers + blockDim.x - 1)/blockDim.x;
  const int numExclusionsTB= (numExclusions + blockDim.x - 1)/blockDim.x;
  const int numCrosstermsTB= (numCrossterms + blockDim.x - 1)/blockDim.x;

  // Each thread computes single bonded interaction.
  // Each thread block computes single bonded type
  double energy;
  int energyIndex;

  if (doEnergy) {
    energy = 0.0;
    energyIndex = -1;
  }

#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double> virial;
  int virialIndex;
  if (doVirial) {
    virial.xx = 0.0;
    virial.xy = 0.0;
    virial.xz = 0.0;
    virial.yx = 0.0;
    virial.yy = 0.0;
    virial.yz = 0.0;
    virial.zx = 0.0;
    virial.zy = 0.0;
    virial.zz = 0.0;
    virialIndex = ComputeBondedCUDAKernel::normalVirialIndex_XX;
  }
#endif

  if (indexTB < numBondsTB) {
    int i = threadIdx.x + indexTB*blockDim.x;
    if (doEnergy) {
      energyIndex = ComputeBondedCUDAKernel::energyIndex_BOND;
    }
    if (i < numBonds) {
      bondForce<T, doEnergy, doVirial>
      (i, bonds, bondValues, xyzq,
        stride, lata, latb, latc,
        force, energy, virial);
    }
    goto reduce;
  }
  indexTB -= numBondsTB;

  if (indexTB < numAnglesTB) {
    int i = threadIdx.x + indexTB*blockDim.x;
    if (doEnergy) {
      energyIndex = ComputeBondedCUDAKernel::energyIndex_ANGLE;
    }
    if (i < numAngles) {
      angleForce<T, doEnergy, doVirial>
      (i, angles, angleValues, xyzq, stride,
        lata, latb, latc,
        force, energy, virial);
    }
    goto reduce;
  }
  indexTB -= numAnglesTB;

  if (indexTB < numDihedralsTB) {
    int i = threadIdx.x + indexTB*blockDim.x;
    if (doEnergy) {
      energyIndex = ComputeBondedCUDAKernel::energyIndex_DIHEDRAL;
    }
    if (doVirial) {
      virialIndex = ComputeBondedCUDAKernel::amdDiheVirialIndex_XX;
    }
    if (i < numDihedrals) {
      diheForce<T, doEnergy, doVirial>
      (i, dihedrals, dihedralValues, xyzq, stride,
        lata, latb, latc,
        force, energy, virial);
    }
    goto reduce;
  }
  indexTB -= numDihedralsTB;

  if (indexTB < numImpropersTB) {
    int i = threadIdx.x + indexTB*blockDim.x;
    if (doEnergy) {
      energyIndex = ComputeBondedCUDAKernel::energyIndex_IMPROPER;
    }
    if (i < numImpropers) {
      diheForce<T, doEnergy, doVirial>
      (i, impropers, improperValues, xyzq, stride,
        lata, latb, latc,
        force, energy, virial);
    }
    goto reduce;
  }
  indexTB -= numImpropersTB;

  if (indexTB < numExclusionsTB) {
    int i = threadIdx.x + indexTB*blockDim.x;
    if (doEnergy) {
      energyIndex = ComputeBondedCUDAKernel::energyIndex_ELECT_SLOW;
    }
    if (doVirial) {
      virialIndex = ComputeBondedCUDAKernel::slowVirialIndex_XX;
    }
    if (i < numExclusions) {
      exclusionForce<T, doEnergy, doVirial>
      (i, exclusions, r2_delta, r2_delta_expc,
#if __CUDA_ARCH__ >= 350
        r2_table, exclusionTable,
#else
        r2_table_tex, exclusionTableTex,
#endif
        xyzq, stride, lata, latb, latc, cutoff2,
        forceSlow, energy, virial);
    }
    goto reduce;
  }
  indexTB -= numExclusionsTB;

  if (indexTB < numCrosstermsTB) {
    int i = threadIdx.x + indexTB*blockDim.x;
    if (doEnergy) {
      energyIndex = ComputeBondedCUDAKernel::energyIndex_CROSSTERM;
    }
    if (doVirial) {
      virialIndex = ComputeBondedCUDAKernel::amdDiheVirialIndex_XX;
    }
    if (i < numCrossterms) {
      crosstermForce<T, doEnergy, doVirial>
      (i, crossterms, crosstermValues,
        xyzq, stride, lata, latb, latc,
        force, energy, virial);
    }
    goto reduce;
  }
  // indexTB -= numCrosstermsTB;

  reduce:

  // Write energies to global memory
  if (doEnergy) {
    // energyIndex is constant within thread-block
    __shared__ double shEnergy[BONDEDFORCESKERNEL_NUM_WARP];
#pragma unroll
    for (int i=16;i >= 1;i/=2) {
      energy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energy, i, 32);
    }
    int laneID = (threadIdx.x & (WARPSIZE - 1));
    int warpID = threadIdx.x / WARPSIZE;
    if (laneID == 0) {
      shEnergy[warpID] = energy;
    }
    BLOCK_SYNC;
    if (warpID == 0) {
      energy = (laneID < BONDEDFORCESKERNEL_NUM_WARP) ? shEnergy[laneID] : 0.0;
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        energy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energy, i, 32);
      }
      if (laneID == 0) {
        atomicAdd(&energies_virials[energyIndex], energy);
      }
    }
  }

  // Write virials to global memory
#ifdef WRITE_FULL_VIRIALS
  if (doVirial) {
#pragma unroll
    for (int i=16;i >= 1;i/=2) {
      virial.xx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.xx, i, 32);
      virial.xy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.xy, i, 32);
      virial.xz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.xz, i, 32);
      virial.yx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.yx, i, 32);
      virial.yy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.yy, i, 32);
      virial.yz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.yz, i, 32);
      virial.zx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.zx, i, 32);
      virial.zy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.zy, i, 32);
      virial.zz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.zz, i, 32);
    }
    __shared__ ComputeBondedCUDAKernel::BondedVirial<double> shVirial[BONDEDFORCESKERNEL_NUM_WARP];
    int laneID = (threadIdx.x & (WARPSIZE - 1));
    int warpID = threadIdx.x / WARPSIZE;
    if (laneID == 0) {
      shVirial[warpID] = virial;
    }
    BLOCK_SYNC;

    if (warpID == 0) {
      virial.xx = 0.0;
      virial.xy = 0.0;
      virial.xz = 0.0;
      virial.yx = 0.0;
      virial.yy = 0.0;
      virial.yz = 0.0;
      virial.zx = 0.0;
      virial.zy = 0.0;
      virial.zz = 0.0;
      if (laneID < BONDEDFORCESKERNEL_NUM_WARP) virial = shVirial[laneID];
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        virial.xx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.xx, i, 32);
        virial.xy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.xy, i, 32);
        virial.xz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.xz, i, 32);
        virial.yx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.yx, i, 32);
        virial.yy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.yy, i, 32);
        virial.yz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.yz, i, 32);
        virial.zx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.zx, i, 32);
        virial.zy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.zy, i, 32);
        virial.zz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virial.zz, i, 32);
      }   

      if (laneID == 0) {
#ifdef USE_FP_VIRIAL
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 0], llitoulli(virial.xx*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 1], llitoulli(virial.xy*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 2], llitoulli(virial.xz*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 3], llitoulli(virial.yx*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 4], llitoulli(virial.yy*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 5], llitoulli(virial.yz*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 6], llitoulli(virial.zx*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 7], llitoulli(virial.zy*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[virialIndex + 8], llitoulli(virial.zz*double_to_virial));
#else
        atomicAdd(&energies_virials[virialIndex + 0], virial.xx);
        atomicAdd(&energies_virials[virialIndex + 1], virial.xy);
        atomicAdd(&energies_virials[virialIndex + 2], virial.xz);
        atomicAdd(&energies_virials[virialIndex + 3], virial.yx);
        atomicAdd(&energies_virials[virialIndex + 4], virial.yy);
        atomicAdd(&energies_virials[virialIndex + 5], virial.yz);
        atomicAdd(&energies_virials[virialIndex + 6], virial.zx);
        atomicAdd(&energies_virials[virialIndex + 7], virial.zy);
        atomicAdd(&energies_virials[virialIndex + 8], virial.zz);
#endif
      }
    }
  }
#endif

}

template <typename T, bool doEnergy, bool doVirial, bool doElect>
__global__ void modifiedExclusionForcesKernel(
  const int start,

  const int numModifiedExclusions,
  const CudaExclusion* __restrict__ modifiedExclusions,

  const bool doSlow,
  const float one_scale14,                // 1 - scale14
  const float cutoff2,
  const int vdwCoefTableWidth,
  const float2* __restrict__ vdwCoefTable,
  cudaTextureObject_t vdwCoefTableTex, 
  cudaTextureObject_t modifiedExclusionForceTableTex, cudaTextureObject_t modifiedExclusionEnergyTableTex,

  const float4* __restrict__ xyzq,
  const int stride,
  const float3 lata, const float3 latb, const float3 latc,
  T* __restrict__ forceNbond, T* __restrict__ forceSlow,
  double* __restrict__ energies_virials
  ) {

  // index
  int i = threadIdx.x + (start + blockIdx.x)*blockDim.x;

  double energyVdw, energyNbond, energySlow;
  if (doEnergy) {
    energyVdw = 0.0;
    if (doElect) {
      energyNbond = 0.0;
      energySlow = 0.0;
    }
  }

#ifdef WRITE_FULL_VIRIALS
  ComputeBondedCUDAKernel::BondedVirial<double> virialNbond;
  ComputeBondedCUDAKernel::BondedVirial<double> virialSlow;
  if (doVirial) {
    virialNbond.xx = 0.0;
    virialNbond.xy = 0.0;
    virialNbond.xz = 0.0;
    virialNbond.yx = 0.0;
    virialNbond.yy = 0.0;
    virialNbond.yz = 0.0;
    virialNbond.zx = 0.0;
    virialNbond.zy = 0.0;
    virialNbond.zz = 0.0;
    if (doElect) {
      virialSlow.xx = 0.0;
      virialSlow.xy = 0.0;
      virialSlow.xz = 0.0;
      virialSlow.yx = 0.0;
      virialSlow.yy = 0.0;
      virialSlow.yz = 0.0;
      virialSlow.zx = 0.0;
      virialSlow.zy = 0.0;
      virialSlow.zz = 0.0;
    }
  }
#endif

  if (i < numModifiedExclusions)
  {
    modifiedExclusionForce<T, doEnergy, doVirial, doElect>
    (i, modifiedExclusions, doSlow, one_scale14, vdwCoefTableWidth,
#if __CUDA_ARCH__ >= 350
      vdwCoefTable,
#else
      vdwCoefTableTex,
#endif
      modifiedExclusionForceTableTex, modifiedExclusionEnergyTableTex,
      xyzq, stride, lata, latb, latc, cutoff2,
      energyVdw, forceNbond, energyNbond,
      forceSlow, energySlow, virialNbond, virialSlow);
  }

  // Write energies to global memory
  if (doEnergy) {
    __shared__ double shEnergyVdw[BONDEDFORCESKERNEL_NUM_WARP];
    __shared__ double shEnergyNbond[(doElect) ? BONDEDFORCESKERNEL_NUM_WARP : 1];
    __shared__ double shEnergySlow[(doElect) ? BONDEDFORCESKERNEL_NUM_WARP : 1];
#pragma unroll
    for (int i=16;i >= 1;i/=2) {
      energyVdw   += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energyVdw, i, 32);
      if (doElect) {
        energyNbond += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energyNbond, i, 32);
        energySlow  += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energySlow, i, 32);
      }
    }
    int laneID = (threadIdx.x & (WARPSIZE - 1));
    int warpID = threadIdx.x / WARPSIZE;
    if (laneID == 0) {
      shEnergyVdw[warpID]   = energyVdw;
      if (doElect) {
        shEnergyNbond[warpID] = energyNbond;
        shEnergySlow[warpID]  = energySlow;
      }
    }
    BLOCK_SYNC;
    if (warpID == 0) {
      energyVdw   = (laneID < BONDEDFORCESKERNEL_NUM_WARP) ? shEnergyVdw[laneID] : 0.0;
      if (doElect) {
        energyNbond = (laneID < BONDEDFORCESKERNEL_NUM_WARP) ? shEnergyNbond[laneID] : 0.0;
        energySlow  = (laneID < BONDEDFORCESKERNEL_NUM_WARP) ? shEnergySlow[laneID] : 0.0;
      }
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        energyVdw   += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energyVdw, i, 32);
        if (doElect) {
          energyNbond += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energyNbond, i, 32);
          energySlow  += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energySlow, i, 32);
        }
      }
      if (laneID == 0) {
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::energyIndex_LJ],         energyVdw);
        if (doElect) {
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::energyIndex_ELECT],      energyNbond);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::energyIndex_ELECT_SLOW], energySlow);
        }
      }
    }
  }

  // Write virials to global memory
#ifdef WRITE_FULL_VIRIALS
  if (doVirial) {
#pragma unroll
    for (int i=16;i >= 1;i/=2) {
      virialNbond.xx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.xx, i, 32);
      virialNbond.xy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.xy, i, 32);
      virialNbond.xz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.xz, i, 32);
      virialNbond.yx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.yx, i, 32);
      virialNbond.yy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.yy, i, 32);
      virialNbond.yz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.yz, i, 32);
      virialNbond.zx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.zx, i, 32);
      virialNbond.zy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.zy, i, 32);
      virialNbond.zz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.zz, i, 32);
      if (doElect && doSlow) {
        virialSlow.xx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.xx, i, 32);
        virialSlow.xy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.xy, i, 32);
        virialSlow.xz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.xz, i, 32);
        virialSlow.yx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.yx, i, 32);
        virialSlow.yy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.yy, i, 32);
        virialSlow.yz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.yz, i, 32);
        virialSlow.zx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.zx, i, 32);
        virialSlow.zy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.zy, i, 32);
        virialSlow.zz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.zz, i, 32);
      }
    }
    __shared__ ComputeBondedCUDAKernel::BondedVirial<double> shVirialNbond[BONDEDFORCESKERNEL_NUM_WARP];
    __shared__ ComputeBondedCUDAKernel::BondedVirial<double> shVirialSlow[(doElect) ? BONDEDFORCESKERNEL_NUM_WARP : 1];
    int laneID = (threadIdx.x & (WARPSIZE - 1));
    int warpID = threadIdx.x / WARPSIZE;
    if (laneID == 0) {
      shVirialNbond[warpID] = virialNbond;
      shVirialSlow[warpID] = virialSlow;
    }
    BLOCK_SYNC;

    virialNbond.xx = 0.0;
    virialNbond.xy = 0.0;
    virialNbond.xz = 0.0;
    virialNbond.yx = 0.0;
    virialNbond.yy = 0.0;
    virialNbond.yz = 0.0;
    virialNbond.zx = 0.0;
    virialNbond.zy = 0.0;
    virialNbond.zz = 0.0;
    if (doElect && doSlow) {
      virialSlow.xx = 0.0;
      virialSlow.xy = 0.0;
      virialSlow.xz = 0.0;
      virialSlow.yx = 0.0;
      virialSlow.yy = 0.0;
      virialSlow.yz = 0.0;
      virialSlow.zx = 0.0;
      virialSlow.zy = 0.0;
      virialSlow.zz = 0.0;
    }

    if (warpID == 0) {
      if (laneID < BONDEDFORCESKERNEL_NUM_WARP) {
        virialNbond = shVirialNbond[laneID];
        virialSlow = shVirialSlow[laneID];
      }
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        virialNbond.xx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.xx, i, 32);
        virialNbond.xy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.xy, i, 32);
        virialNbond.xz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.xz, i, 32);
        virialNbond.yx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.yx, i, 32);
        virialNbond.yy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.yy, i, 32);
        virialNbond.yz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.yz, i, 32);
        virialNbond.zx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.zx, i, 32);
        virialNbond.zy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.zy, i, 32);
        virialNbond.zz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialNbond.zz, i, 32);
        if (doElect && doSlow) {
          virialSlow.xx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.xx, i, 32);
          virialSlow.xy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.xy, i, 32);
          virialSlow.xz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.xz, i, 32);
          virialSlow.yx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.yx, i, 32);
          virialSlow.yy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.yy, i, 32);
          virialSlow.yz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.yz, i, 32);
          virialSlow.zx += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.zx, i, 32);
          virialSlow.zy += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.zy, i, 32);
          virialSlow.zz += WARP_SHUFFLE_XOR(WARP_FULL_MASK, virialSlow.zz, i, 32);
        }
      }

      if (laneID == 0)
      {
#ifdef USE_FP_VIRIAL
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XX], llitoulli(virialNbond.xx*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XY], llitoulli(virialNbond.xy*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XZ], llitoulli(virialNbond.xz*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_YX], llitoulli(virialNbond.yx*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_YY], llitoulli(virialNbond.yy*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_YZ], llitoulli(virialNbond.yz*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_ZX], llitoulli(virialNbond.zx*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_ZY], llitoulli(virialNbond.zy*double_to_virial));
        atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_ZZ], llitoulli(virialNbond.zz*double_to_virial));
        if (doElect && doSlow) {
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XX], llitoulli(virialSlow.xx*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XY], llitoulli(virialSlow.xy*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XZ], llitoulli(virialSlow.xz*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_YX], llitoulli(virialSlow.yx*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_YY], llitoulli(virialSlow.yy*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_YZ], llitoulli(virialSlow.yz*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_ZX], llitoulli(virialSlow.zx*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_ZY], llitoulli(virialSlow.zy*double_to_virial));
          atomicAdd((unsigned long long int *)&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_ZZ], llitoulli(virialSlow.zz*double_to_virial));
        }
#else
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XX], virialNbond.xx);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XY], virialNbond.xy);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XZ], virialNbond.xz);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_YX], virialNbond.yx);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_YY], virialNbond.yy);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_YZ], virialNbond.yz);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_ZX], virialNbond.zx);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_ZY], virialNbond.zy);
        atomicAdd(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_ZZ], virialNbond.zz);
        if (doElect && doSlow) {
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XX], virialSlow.xx);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XY], virialSlow.xy);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XZ], virialSlow.xz);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_YX], virialSlow.yx);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_YY], virialSlow.yy);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_YZ], virialSlow.yz);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_ZX], virialSlow.zx);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_ZY], virialSlow.zy);
          atomicAdd(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_ZZ], virialSlow.zz);
        }
#endif
      }
    }
  }
#endif

}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

//
// Class constructor
//
ComputeBondedCUDAKernel::ComputeBondedCUDAKernel(int deviceID, CudaNonbondedTables& cudaNonbondedTables) :
deviceID(deviceID), cudaNonbondedTables(cudaNonbondedTables) {

  cudaCheck(cudaSetDevice(deviceID));

  tupleData = NULL;
  tupleDataSize = 0;

  numBonds = 0;
  numAngles = 0;
  numDihedrals = 0;
  numImpropers = 0;
  numModifiedExclusions = 0;
  numExclusions = 0;
  numCrossterms = 0;

  bondValues = NULL;
  angleValues = NULL;
  dihedralValues = NULL;
  improperValues = NULL;
  crosstermValues = NULL;

  xyzq = NULL;
  xyzqSize = 0;

  forces = NULL;
  forcesSize = 0;

  allocate_device<double>(&energies_virials, energies_virials_SIZE);
}

//
// Class destructor
//
ComputeBondedCUDAKernel::~ComputeBondedCUDAKernel() {
  cudaCheck(cudaSetDevice(deviceID));

  deallocate_device<double>(&energies_virials);
  // deallocate_device<BondedVirial>(&virial);

  if (tupleData != NULL) deallocate_device<char>(&tupleData);
  if (xyzq != NULL) deallocate_device<float4>(&xyzq);
  if (forces != NULL) deallocate_device<FORCE_TYPE>(&forces);

  if (bondValues != NULL) deallocate_device<CudaBondValue>(&bondValues);
  if (angleValues != NULL) deallocate_device<CudaAngleValue>(&angleValues);
  if (dihedralValues != NULL) deallocate_device<CudaDihedralValue>(&dihedralValues);
  if (improperValues != NULL) deallocate_device<CudaDihedralValue>(&improperValues);
  if (crosstermValues != NULL) deallocate_device<CudaCrosstermValue>(&crosstermValues);
}

void ComputeBondedCUDAKernel::setupBondValues(int numBondValues, CudaBondValue* h_bondValues) {
  allocate_device<CudaBondValue>(&bondValues, numBondValues);
  copy_HtoD_sync<CudaBondValue>(h_bondValues, bondValues, numBondValues);
}

void ComputeBondedCUDAKernel::setupAngleValues(int numAngleValues, CudaAngleValue* h_angleValues) {
  allocate_device<CudaAngleValue>(&angleValues, numAngleValues);
  copy_HtoD_sync<CudaAngleValue>(h_angleValues, angleValues, numAngleValues);
}

void ComputeBondedCUDAKernel::setupDihedralValues(int numDihedralValues, CudaDihedralValue* h_dihedralValues) {
  allocate_device<CudaDihedralValue>(&dihedralValues, numDihedralValues);
  copy_HtoD_sync<CudaDihedralValue>(h_dihedralValues, dihedralValues, numDihedralValues);
}

void ComputeBondedCUDAKernel::setupImproperValues(int numImproperValues, CudaDihedralValue* h_improperValues) {
  allocate_device<CudaDihedralValue>(&improperValues, numImproperValues);
  copy_HtoD_sync<CudaDihedralValue>(h_improperValues, improperValues, numImproperValues);
}

void ComputeBondedCUDAKernel::setupCrosstermValues(int numCrosstermValues, CudaCrosstermValue* h_crosstermValues) {
  allocate_device<CudaCrosstermValue>(&crosstermValues, numCrosstermValues);
  copy_HtoD_sync<CudaCrosstermValue>(h_crosstermValues, crosstermValues, numCrosstermValues);
}

//
// Update bonded lists
//
void ComputeBondedCUDAKernel::update(
  const int numBondsIn,
  const int numAnglesIn,
  const int numDihedralsIn,
  const int numImpropersIn,
  const int numModifiedExclusionsIn,
  const int numExclusionsIn,
  const int numCrosstermsIn,
  const char* h_tupleData,
  cudaStream_t stream) {

  numBonds              = numBondsIn;
  numAngles             = numAnglesIn;
  numDihedrals          = numDihedralsIn;
  numImpropers          = numImpropersIn;
  numModifiedExclusions = numModifiedExclusionsIn;
  numExclusions         = numExclusionsIn;
  numCrossterms         = numCrosstermsIn;

  int numBondsWA     = warpAlign(numBonds);
  int numAnglesWA    = warpAlign(numAngles);
  int numDihedralsWA = warpAlign(numDihedrals);
  int numImpropersWA = warpAlign(numImpropers);
  int numModifiedExclusionsWA = warpAlign(numModifiedExclusions);
  int numExclusionsWA         = warpAlign(numExclusions);
  int numCrosstermsWA         = warpAlign(numCrossterms);

  int sizeTot = numBondsWA*sizeof(CudaBond) + numAnglesWA*sizeof(CudaAngle) + 
  numDihedralsWA*sizeof(CudaDihedral) + numImpropersWA*sizeof(CudaDihedral) +
  numModifiedExclusionsWA*sizeof(CudaExclusion) + numExclusionsWA*sizeof(CudaExclusion) + 
  numCrosstermsWA*sizeof(CudaCrossterm);

  reallocate_device<char>(&tupleData, &tupleDataSize, sizeTot, 1.4f);
  copy_HtoD<char>(h_tupleData, tupleData, sizeTot, stream);

  // Setup pointers
  int pos = 0;
  bonds = (CudaBond *)&tupleData[pos];
  pos += numBondsWA*sizeof(CudaBond);

  angles = (CudaAngle* )&tupleData[pos];
  pos += numAnglesWA*sizeof(CudaAngle);

  dihedrals = (CudaDihedral* )&tupleData[pos];
  pos += numDihedralsWA*sizeof(CudaDihedral);

  impropers = (CudaDihedral* )&tupleData[pos];
  pos += numImpropersWA*sizeof(CudaDihedral);

  modifiedExclusions = (CudaExclusion* )&tupleData[pos];
  pos += numModifiedExclusionsWA*sizeof(CudaExclusion);

  exclusions = (CudaExclusion* )&tupleData[pos];
  pos += numExclusionsWA*sizeof(CudaExclusion);

  crossterms = (CudaCrossterm* )&tupleData[pos];
  pos += numCrosstermsWA*sizeof(CudaCrossterm);
}

//
// Return stride for force array
//
int ComputeBondedCUDAKernel::getForceStride(const int atomStorageSize) {
#ifdef USE_STRIDED_FORCE
  // Align stride to 256 bytes
  return ((atomStorageSize*sizeof(FORCE_TYPE) - 1)/256 + 1)*256/sizeof(FORCE_TYPE);
#else
  return 1;
#endif
}

//
// Return size of single force array
//
int ComputeBondedCUDAKernel::getForceSize(const int atomStorageSize) {
#ifdef USE_STRIDED_FORCE
  return (3*getForceStride(atomStorageSize));
#else
  return (3*atomStorageSize);
#endif
}

//
// Return size of the all force arrays
//
int ComputeBondedCUDAKernel::getAllForceSize(const int atomStorageSize, const bool doSlow) {

  int forceSize = getForceSize(atomStorageSize);

  if (numModifiedExclusions > 0 || numExclusions > 0) {
    if (doSlow) {
      // All three force arrays [normal, nbond, slow]
      forceSize *= 3;
    } else {
      // Two force arrays [normal, nbond]
      forceSize *= 2;
    }
  }

  return forceSize;
}

//
// Compute bonded forces
//
void ComputeBondedCUDAKernel::bondedForce(
  const double scale14, const int atomStorageSize,
  const bool doEnergy, const bool doVirial, const bool doSlow,
  const float3 lata, const float3 latb, const float3 latc,
  const float cutoff2, const float r2_delta, const int r2_delta_expc,
  const float4* h_xyzq, FORCE_TYPE* h_forces,
  double *h_energies_virials,
  cudaStream_t stream) {

  int forceStorageSize = getAllForceSize(atomStorageSize, true);
  int forceCopySize = getAllForceSize(atomStorageSize, doSlow);
  int forceStride = getForceStride(atomStorageSize);

  int forceSize = getForceSize(atomStorageSize);
  int startNbond = forceSize;
  int startSlow = 2*forceSize;

  // Re-allocate coordinate and force arrays if neccessary
  reallocate_device<float4>(&xyzq, &xyzqSize, atomStorageSize, 1.4f);
  reallocate_device<FORCE_TYPE>(&forces, &forcesSize, forceStorageSize, 1.4f);

  // Copy coordinates to device
  copy_HtoD<float4>(h_xyzq, xyzq, atomStorageSize, stream);

  // Clear force array
  clear_device_array<FORCE_TYPE>(forces, forceCopySize, stream);
  if (doEnergy || doVirial) {
    clear_device_array<double>(energies_virials, energies_virials_SIZE, stream);
  }

  float one_scale14 = (float)(1.0 - scale14);

  // If doSlow = false, these exclusions are not calculated
  int numExclusionsDoSlow = doSlow ? numExclusions : 0;

  int nthread = BONDEDFORCESKERNEL_NUM_WARP * WARPSIZE;

  int numBondsTB     = (numBonds + nthread - 1)/nthread;
  int numAnglesTB    = (numAngles + nthread - 1)/nthread;
  int numDihedralsTB = (numDihedrals + nthread - 1)/nthread;
  int numImpropersTB = (numImpropers + nthread - 1)/nthread;
  int numExclusionsTB= (numExclusionsDoSlow + nthread - 1)/nthread;
  int numCrosstermsTB= (numCrossterms + nthread - 1)/nthread;

  int nblock = numBondsTB + numAnglesTB + numDihedralsTB + numImpropersTB + 
  numExclusionsTB + numCrosstermsTB;
  int shmem_size = 0;

  // printf("%d %d %d %d %d %d nblock %d\n",
  //   numBonds, numAngles, numDihedrals, numImpropers, numModifiedExclusions, numExclusions, nblock);

  int max_nblock = deviceCUDA->getMaxNumBlocks();

  int start = 0;
  while (start < nblock)
  {
    int nleft = nblock - start;
    int nblock_use = min(max_nblock, nleft);

#define CALL(DOENERGY, DOVIRIAL) \
  bondedForcesKernel<FORCE_TYPE, DOENERGY, DOVIRIAL> \
  <<< nblock_use, nthread, shmem_size, stream >>> \
  (start, numBonds, bonds, bondValues, \
    numAngles, angles, angleValues, \
    numDihedrals, dihedrals, dihedralValues, \
    numImpropers, impropers, improperValues, \
    numExclusionsDoSlow, exclusions, \
    numCrossterms, crossterms, crosstermValues, \
    cutoff2, \
    r2_delta, r2_delta_expc, \
    cudaNonbondedTables.get_r2_table(), cudaNonbondedTables.getExclusionTable(), \
    cudaNonbondedTables.get_r2_table_tex(), cudaNonbondedTables.getExclusionTableTex(), \
    xyzq, forceStride, \
    lata, latb, latc, \
    forces, &forces[startSlow], energies_virials);

    if (!doEnergy && !doVirial) CALL(0, 0);
    if (!doEnergy && doVirial)  CALL(0, 1);
    if (doEnergy && !doVirial)  CALL(1, 0);
    if (doEnergy && doVirial)   CALL(1, 1);

#undef CALL
    cudaCheck(cudaGetLastError());

    start += nblock_use;
  }

  nthread = BONDEDFORCESKERNEL_NUM_WARP * WARPSIZE;
  nblock = (numModifiedExclusions + nthread - 1)/nthread;

  bool doElect = (one_scale14 == 0.0f) ? false : true;

  start = 0;
  while (start < nblock)
  {
    int nleft = nblock - start;
    int nblock_use = min(max_nblock, nleft);

#define CALL(DOENERGY, DOVIRIAL, DOELECT) \
  modifiedExclusionForcesKernel<FORCE_TYPE, DOENERGY, DOVIRIAL, DOELECT> \
  <<< nblock_use, nthread, shmem_size, stream >>> \
  (start, numModifiedExclusions, modifiedExclusions, \
    doSlow, one_scale14, cutoff2, \
    cudaNonbondedTables.getVdwCoefTableWidth(), cudaNonbondedTables.getExclusionVdwCoefTable(), \
    cudaNonbondedTables.getExclusionVdwCoefTableTex(), \
    cudaNonbondedTables.getModifiedExclusionForceTableTex(), cudaNonbondedTables.getModifiedExclusionEnergyTableTex(), \
    xyzq, forceStride, lata, latb, latc, \
    &forces[startNbond], &forces[startSlow], energies_virials);

    if (!doEnergy && !doVirial && !doElect) CALL(0, 0, 0);
    if (!doEnergy && doVirial && !doElect)  CALL(0, 1, 0);
    if (doEnergy && !doVirial && !doElect)  CALL(1, 0, 0);
    if (doEnergy && doVirial && !doElect)   CALL(1, 1, 0);

    if (!doEnergy && !doVirial && doElect)  CALL(0, 0, 1);
    if (!doEnergy && doVirial && doElect)   CALL(0, 1, 1);
    if (doEnergy && !doVirial && doElect)   CALL(1, 0, 1);
    if (doEnergy && doVirial && doElect)    CALL(1, 1, 1);

#undef CALL
    cudaCheck(cudaGetLastError());

    start += nblock_use;
  }

  copy_DtoH<FORCE_TYPE>(forces, h_forces, forceCopySize, stream);
  if (doEnergy || doVirial) {
    copy_DtoH<double>(energies_virials, h_energies_virials, energies_virials_SIZE, stream);
  }

}
