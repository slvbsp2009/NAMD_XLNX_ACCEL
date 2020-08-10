#include <cuda.h>
#include <cub/cub.cuh>
#include "CudaComputeNonbondedKernel.h"
#include "CudaTileListKernel.h"
#include "DeviceCUDA.h"
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;

#define OVERALLOC 1.2f

void NAMD_die(const char *);

#define MAX_CONST_EXCLUSIONS 2048  // cache size is 8k
__constant__ unsigned int constExclusions[MAX_CONST_EXCLUSIONS];

#define NONBONDKERNEL_NUM_WARP 4

template<bool doEnergy, bool doSlow>
__device__ __forceinline__
void calcForceEnergy(const float r2, const float qi, const float qj,
  const float dx, const float dy, const float dz,
  const int vdwtypei, const int vdwtypej, const float2* __restrict__ vdwCoefTable,
  cudaTextureObject_t vdwCoefTableTex, 
  cudaTextureObject_t forceTableTex, cudaTextureObject_t energyTableTex,
  float3& iforce, float3& iforceSlow, float3& jforce, float3& jforceSlow,
  float& energyVdw, float& energyElec, float& energySlow) {

  int vdwIndex = vdwtypej + vdwtypei;
#if __CUDA_ARCH__ >= 350
  float2 ljab = __ldg(&vdwCoefTable[vdwIndex]);
#else
  float2 ljab = tex1Dfetch<float2>(vdwCoefTableTex, vdwIndex);
#endif

  float rinv = rsqrtf(r2);
  float4 ei;
  float4 fi = tex1D<float4>(forceTableTex, rinv);
  if (doEnergy) ei = tex1D<float4>(energyTableTex, rinv);

  float fSlow = qi * qj;
  float f = ljab.x * fi.z + ljab.y * fi.y + fSlow * fi.x;

  if (doEnergy) {
    energyVdw  += ljab.x * ei.z + ljab.y * ei.y;
    energyElec += fSlow * ei.x;
    if (doSlow) energySlow += fSlow * ei.w;
  }
  if (doSlow) fSlow *= fi.w;

  float fx = dx * f;
  float fy = dy * f;
  float fz = dz * f;
  iforce.x += fx;
  iforce.y += fy;
  iforce.z += fz;
  jforce.x -= fx;
  jforce.y -= fy;
  jforce.z -= fz;

  if (doSlow) {
    float fxSlow = dx * fSlow;
    float fySlow = dy * fSlow;
    float fzSlow = dz * fSlow;
    iforceSlow.x += fxSlow;
    iforceSlow.y += fySlow;
    iforceSlow.z += fzSlow;
    jforceSlow.x -= fxSlow;
    jforceSlow.y -= fySlow;
    jforceSlow.z -= fzSlow;
  }
}

template<bool doSlow>
__device__ __forceinline__
void storeForces(const int pos, const float3 force, const float3 forceSlow,
  float4* __restrict__ devForces, float4* __restrict__ devForcesSlow) {
  atomicAdd(&devForces[pos].x, force.x);
  atomicAdd(&devForces[pos].y, force.y);
  atomicAdd(&devForces[pos].z, force.z);
  if (doSlow) {
    atomicAdd(&devForcesSlow[pos].x, forceSlow.x);
    atomicAdd(&devForcesSlow[pos].y, forceSlow.y);
    atomicAdd(&devForcesSlow[pos].z, forceSlow.z);
  }
}

template<bool doSlow>
__device__ __forceinline__
void storeForces(const int pos, const float3 force, const float3 forceSlow,
                 float* __restrict__ devForces_x, 
                 float* __restrict__ devForces_y, 
                 float* __restrict__ devForces_z,
                 float* __restrict__ devForcesSlow_x, 
                 float* __restrict__ devForcesSlow_y, 
                 float* __restrict__ devForcesSlow_z)
{
  atomicAdd(&devForces_x[pos], force.x);
  atomicAdd(&devForces_y[pos], force.y);
  atomicAdd(&devForces_z[pos], force.z);
  if (doSlow) {
    atomicAdd(&devForcesSlow_x[pos], forceSlow.x);
    atomicAdd(&devForcesSlow_y[pos], forceSlow.y);
    atomicAdd(&devForcesSlow_z[pos], forceSlow.z);
  }
}

template<bool doSlow>
__device__ __forceinline__
void storeForces(const int pos, const float3 force, const float3 forceSlow,
  float3* __restrict__ forces, float3* __restrict__ forcesSlow) {
  atomicAdd(&forces[pos].x, force.x);
  atomicAdd(&forces[pos].y, force.y);
  atomicAdd(&forces[pos].z, force.z);
  if (doSlow) {
    atomicAdd(&forcesSlow[pos].x, forceSlow.x);
    atomicAdd(&forcesSlow[pos].y, forceSlow.y);
    atomicAdd(&forcesSlow[pos].z, forceSlow.z);
  }
}

template<bool doPairlist>
__device__ __forceinline__
void shuffleNext(float& xyzq_j_w, int& vdwtypej, int& jatomIndex, int& jexclMaxdiff, int& jexclIndex) {
  xyzq_j_w = WARP_SHUFFLE(WARP_FULL_MASK, xyzq_j_w, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);
  vdwtypej = WARP_SHUFFLE(WARP_FULL_MASK, vdwtypej, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);
  if (doPairlist) {
    jatomIndex   = WARP_SHUFFLE(WARP_FULL_MASK, jatomIndex, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);    
    jexclIndex   = WARP_SHUFFLE(WARP_FULL_MASK, jexclIndex, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);
    jexclMaxdiff = WARP_SHUFFLE(WARP_FULL_MASK, jexclMaxdiff, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);
  }
}

template<bool doPairlist>
__device__ __forceinline__
void shuffleNext(float& xyzq_j_w, int& vdwtypej, int& jatomIndex) {
  xyzq_j_w = WARP_SHUFFLE(WARP_FULL_MASK, xyzq_j_w, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);
  vdwtypej = WARP_SHUFFLE(WARP_FULL_MASK, vdwtypej, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);
  if (doPairlist) {
    jatomIndex   = WARP_SHUFFLE(WARP_FULL_MASK, jatomIndex, (threadIdx.x+1) & (WARPSIZE-1), WARPSIZE);    
  }
}

template<bool doSlow>
__device__ __forceinline__
void shuffleNext(float3& jforce, float3& jforceSlow) {
  jforce.x = WARP_SHUFFLE(WARP_FULL_MASK, jforce.x, (threadIdx.x+1)&(WARPSIZE-1), WARPSIZE);
  jforce.y = WARP_SHUFFLE(WARP_FULL_MASK, jforce.y, (threadIdx.x+1)&(WARPSIZE-1), WARPSIZE);
  jforce.z = WARP_SHUFFLE(WARP_FULL_MASK, jforce.z, (threadIdx.x+1)&(WARPSIZE-1), WARPSIZE);
  if (doSlow) {
    jforceSlow.x = WARP_SHUFFLE(WARP_FULL_MASK, jforceSlow.x, (threadIdx.x+1)&(WARPSIZE-1), WARPSIZE);
    jforceSlow.y = WARP_SHUFFLE(WARP_FULL_MASK, jforceSlow.y, (threadIdx.x+1)&(WARPSIZE-1), WARPSIZE);
    jforceSlow.z = WARP_SHUFFLE(WARP_FULL_MASK, jforceSlow.z, (threadIdx.x+1)&(WARPSIZE-1), WARPSIZE);
  }
}

//#define USE_NEW_EXCL_METHOD

//
// Returns the lower estimate for the distance between a bounding box and a set of atoms
//
__device__ __forceinline__ float distsq(const BoundingBox a, const float4 b) {
  float dx = max(0.0f, fabsf(a.x - b.x) - a.wx);
  float dy = max(0.0f, fabsf(a.y - b.y) - a.wy);
  float dz = max(0.0f, fabsf(a.z - b.z) - a.wz);
  float r2 = dx*dx + dy*dy + dz*dz;
  return r2;
}

#define LARGE_FLOAT (float)(1.0e10)

//
// Nonbonded force kernel
//
template <bool doEnergy, bool doVirial, bool doSlow, bool doPairlist, bool doStreaming>
__global__ void
__launch_bounds__(WARPSIZE*NONBONDKERNEL_NUM_WARP,
  doPairlist ? (10) : (doEnergy ? (10) : (10) )
  )
nonbondedForceKernel(
  const int start, const int numTileLists,
  const TileList* __restrict__ tileLists, TileExcl* __restrict__ tileExcls,
  const int* __restrict__ tileJatomStart,
  const int vdwCoefTableWidth, const float2* __restrict__ vdwCoefTable, const int* __restrict__ vdwTypes,
  const float3 lata, const float3 latb, const float3 latc,
  const float4* __restrict__ xyzq, const float cutoff2,
  cudaTextureObject_t vdwCoefTableTex,
  cudaTextureObject_t forceTableTex, cudaTextureObject_t energyTableTex,
  // ----------
  // doPairlist
  const int atomStorageSize, const float plcutoff2, const PatchPairRecord* __restrict__ patchPairs,
  const int* __restrict__ atomIndex,
  const int2* __restrict__ exclIndexMaxDiff, const unsigned int* __restrict__ overflowExclusions,
  unsigned int* __restrict__ tileListDepth, int* __restrict__ tileListOrder,
  int* __restrict__ jtiles, TileListStat* __restrict__ tileListStat,
  const BoundingBox* __restrict__ boundingBoxes,
#ifdef USE_NEW_EXCL_METHOD
  const int* __restrict__ minmaxExclAtom,
#endif
  // ----------
  float4* __restrict__ devForces, float4* __restrict__ devForcesSlow,
  float * __restrict__ devForce_x,
  float * __restrict__ devForce_y,
  float * __restrict__ devForce_z,
  float * __restrict__ devForce_w,
  float * __restrict__ devForceSlow_x,
  float * __restrict__ devForceSlow_y,
  float * __restrict__ devForceSlow_z,
  float * __restrict__ devForceSlow_w,                     
  // ---- USE_STREAMING_FORCES ----
  const int numPatches,
  unsigned int* __restrict__ patchNumCount,
  const CudaPatchRecord* __restrict__ cudaPatches,
  float4* __restrict__ mapForces, float4* __restrict__ mapForcesSlow,
  int* __restrict__ mapPatchReadyQueue,
  int* __restrict__ outputOrder,
  // ------------------------------
  TileListVirialEnergy* __restrict__ virialEnergy) {

  // Single warp takes care of one list of tiles
  // for (int itileList = (threadIdx.x + blockDim.x*blockIdx.x)/WARPSIZE;itileList < numTileLists;itileList += blockDim.x*gridDim.x/WARPSIZE)
  int itileList = start + threadIdx.x/WARPSIZE + blockDim.x/WARPSIZE*blockIdx.x;
  if (itileList < numTileLists)
  {

    float3 iforce;
    float3 iforceSlow;
    float energyVdw, energyElec, energySlow;
    int nexcluded;
    unsigned int itileListLen;
    int2 patchInd;
    int2 patchNumList;
    __shared__ float4 s_xyzq[NONBONDKERNEL_NUM_WARP][WARPSIZE];
    __shared__ int s_vdwtypej[NONBONDKERNEL_NUM_WARP][WARPSIZE];
    __shared__ float3 s_jforce[NONBONDKERNEL_NUM_WARP][WARPSIZE];
    __shared__ float3 s_jforceSlow[NONBONDKERNEL_NUM_WARP][WARPSIZE];
    __shared__ int s_jatomIndex[NONBONDKERNEL_NUM_WARP][WARPSIZE];            

    // Start computation
    {
      // Warp index (0...warpsize-1)
      const int wid = threadIdx.x % WARPSIZE;
      const int iwarp = threadIdx.x / WARPSIZE;

      TileList tmp = tileLists[itileList];
      int iatomStart = tmp.iatomStart;
      int jtileStart = tmp.jtileStart;
      int jtileEnd   = tmp.jtileEnd;
      patchInd     = tmp.patchInd;
      patchNumList = tmp.patchNumList;

      float shx = tmp.offsetXYZ.x*lata.x + tmp.offsetXYZ.y*latb.x + tmp.offsetXYZ.z*latc.x;
      float shy = tmp.offsetXYZ.x*lata.y + tmp.offsetXYZ.y*latb.y + tmp.offsetXYZ.z*latc.y;
      float shz = tmp.offsetXYZ.x*lata.z + tmp.offsetXYZ.y*latb.z + tmp.offsetXYZ.z*latc.z;

      // DH - set zeroShift flag if magnitude of shift vector is zero
      bool zeroShift = ! (shx*shx + shy*shy + shz*shz > 0);

      int iatomSize, iatomFreeSize, jatomSize, jatomFreeSize;
      if (doPairlist) {
        PatchPairRecord PPStmp = patchPairs[itileList];
        iatomSize     = PPStmp.iatomSize;
        iatomFreeSize = PPStmp.iatomFreeSize;
        jatomSize     = PPStmp.jatomSize;
        jatomFreeSize = PPStmp.jatomFreeSize;
      }

      // Write to global memory here to avoid register spilling
      if (doVirial) {
        if (wid == 0) {
          virialEnergy[itileList].shx = shx;
          virialEnergy[itileList].shy = shy;
          virialEnergy[itileList].shz = shz;
        }
      }

      // Load i-atom data (and shift coordinates)
      float4 xyzq_i = xyzq[iatomStart + wid];
      xyzq_i.x += shx;
      xyzq_i.y += shy;
      xyzq_i.z += shz;
      int vdwtypei = vdwTypes[iatomStart + wid]*vdwCoefTableWidth;

      // Load i-atom data (and shift coordinates)
      BoundingBox boundingBoxI;
      if (doPairlist) {
        boundingBoxI = boundingBoxes[iatomStart/WARPSIZE];
        boundingBoxI.x += shx;
        boundingBoxI.y += shy;
        boundingBoxI.z += shz;
      }

      // Get i-atom global index
#ifdef USE_NEW_EXCL_METHOD
      int iatomIndex, minExclAtom, maxExclAtom;
#else
      int iatomIndex;
#endif
      if (doPairlist) {
#ifdef USE_NEW_EXCL_METHOD
        iatomIndex = atomIndex[iatomStart + wid];
        int2 tmp = minmaxExclAtom[iatomStart + wid];
        minExclAtom = tmp.x;
        maxExclAtom = tmp.y;
#else
        iatomIndex = atomIndex[iatomStart + wid];
#endif
      }

      // i-forces in registers
      // float3 iforce;
      iforce.x = 0.0f;
      iforce.y = 0.0f;
      iforce.z = 0.0f;

      // float3 iforceSlow;
      if (doSlow) {
        iforceSlow.x = 0.0f;
        iforceSlow.y = 0.0f;
        iforceSlow.z = 0.0f;
      }

      // float energyVdw, energyElec, energySlow;
      if (doEnergy) {
        energyVdw = 0.0f;
        energyElec = 0.0f;
        if (doSlow) energySlow = 0.0f;
      }

      // Number of exclusions
      // NOTE: Lowest bit is used as indicator bit for tile pairs:
      //       bit 0 tile has no atoms within pairlist cutoff
      //       bit 1 tile has atoms within pairlist cutoff
      // int nexcluded;
      if (doPairlist) nexcluded = 0;

      // Number of i loops and free atoms
      int nfreei;
      if (doPairlist) {
        int nloopi = min(iatomSize - iatomStart, WARPSIZE);
        nfreei = max(iatomFreeSize - iatomStart, 0);
        if (wid >= nloopi) {
          xyzq_i.x = -LARGE_FLOAT;
          xyzq_i.y = -LARGE_FLOAT;
          xyzq_i.z = -LARGE_FLOAT;
        }
      }

      // tile list stuff
      // int itileListLen;
      // int minJatomStart;
      if (doPairlist) {
        // minJatomStart = tileJatomStart[jtileStart];
        itileListLen = 0;
      }

      // Exclusion index and maxdiff
      int iexclIndex, iexclMaxdiff;
      if (doPairlist) {
        int2 tmp = exclIndexMaxDiff[iatomStart + wid];
        iexclIndex   = tmp.x;
        iexclMaxdiff = tmp.y;
      }

      for (int jtile=jtileStart;jtile <= jtileEnd;jtile++) {

        // Load j-atom starting index and exclusion mask
        int jatomStart = tileJatomStart[jtile];

        float4 xyzq_j = xyzq[jatomStart + wid];
        WARP_SYNC(WARP_FULL_MASK);        

        // Check for early bail
        if (doPairlist) {
          float r2bb = distsq(boundingBoxI, xyzq_j);
          if (WARP_ALL(WARP_FULL_MASK, r2bb > plcutoff2)) continue;
        }
        unsigned int excl = (doPairlist) ? 0 : tileExcls[jtile].excl[wid];
        int vdwtypej = vdwTypes[jatomStart + wid];
        s_vdwtypej[iwarp][wid] = vdwtypej;

        // Get i-atom global index
        if (doPairlist) {
          s_jatomIndex[iwarp][wid] = atomIndex[jatomStart + wid];
        }

        // Number of j loops and free atoms
        int nfreej;
        if (doPairlist) {
          int nloopj = min(jatomSize - jatomStart, WARPSIZE);
          nfreej = max(jatomFreeSize - jatomStart, 0);
          //if (nfreei == 0 && nfreej == 0) continue;
          if (wid >= nloopj) {
            xyzq_j.x = LARGE_FLOAT;
            xyzq_j.y = LARGE_FLOAT;
            xyzq_j.z = LARGE_FLOAT;
          }
        }

        s_xyzq[iwarp][wid] = xyzq_j; 

        // DH - self requires that zeroShift is also set
        const bool self = zeroShift && (iatomStart == jatomStart);
        const int modval = (self) ? 2*WARPSIZE-1 : WARPSIZE-1;

        s_jforce[iwarp][wid] = make_float3(0.0f, 0.0f, 0.0f);
        if (doSlow)
          s_jforceSlow[iwarp][wid] = make_float3(0.0f, 0.0f, 0.0f); 
        WARP_SYNC(WARP_FULL_MASK);
        

        int t = (self) ? 1 : 0;

        if (doPairlist) {
          // Build pair list
          // NOTE: Pairlist update, we must also include the diagonal since this is used
          //       in GBIS phase 2.
          // Clear the lowest (indicator) bit
          nexcluded &= (~1);

          // For self tiles, do the diagonal term (t=0).
          // NOTE: No energies are computed here, since this self-diagonal term is only for GBIS phase 2
          if (self) {
            int j = (0 + wid) & modval;
            xyzq_j = s_xyzq[iwarp][j];
            float dx = xyzq_j.x - xyzq_i.x;
            float dy = xyzq_j.y - xyzq_i.y;
            float dz = xyzq_j.z - xyzq_i.z;       

            float r2 = dx*dx + dy*dy + dz*dz;

            if (j < WARPSIZE && r2 < plcutoff2) {
              // We have atom pair within the pairlist cutoff => Set indicator bit
              nexcluded |= 1;
            }
          }

          for (;t < WARPSIZE;t++) {
            int j = (t + wid) & modval;

            excl >>= 1;
            if (j < WARPSIZE ) {
              xyzq_j = s_xyzq[iwarp][j];
              float dx = xyzq_j.x - xyzq_i.x;
              float dy = xyzq_j.y - xyzq_i.y;
              float dz = xyzq_j.z - xyzq_i.z;
              float r2 = dx*dx + dy*dy + dz*dz;
              // We have atom pair within the pairlist cutoff => Set indicator bit
              if(r2 < plcutoff2){
                nexcluded |= 1;
                if (j < nfreej || wid < nfreei) {
                  bool excluded = false;
                  int indexdiff = s_jatomIndex[iwarp][j] - iatomIndex;
                  if ( abs(indexdiff) <= iexclMaxdiff) {
                    indexdiff += iexclIndex;
                    int indexword = ((unsigned int) indexdiff) >> 5;

                    if ( indexword < MAX_CONST_EXCLUSIONS ) {
                      indexword = constExclusions[indexword];
                    } else {
                      indexword = overflowExclusions[indexword];
                    }

                    excluded = ((indexword & (1<<(indexdiff&31))) != 0);
                  }
                  if (excluded) nexcluded += 2;
                  if (!excluded) excl |= 0x80000000;
                  if (!excluded && r2 < cutoff2) {
                    calcForceEnergy<doEnergy, doSlow>(
                      r2, xyzq_i.w, xyzq_j.w, dx, dy, dz,
                      vdwtypei, s_vdwtypej[iwarp][j],
                      vdwCoefTable, vdwCoefTableTex, forceTableTex, energyTableTex,
                      iforce, iforceSlow,
                      s_jforce[iwarp][j], s_jforceSlow[iwarp][j],
                      energyVdw,
                      energyElec, energySlow);
                  }
                }
              }
            }
            WARP_SYNC(WARP_FULL_MASK);
         } // t
       } else {
          // Just compute forces
          if (self) {
            excl >>= 1;
          }
          for (;t < WARPSIZE;t++) {
            if ((excl & 1)) {
              xyzq_j = s_xyzq[iwarp][(wid+t) & (WARPSIZE-1)];
              float dx = xyzq_j.x - xyzq_i.x;
              float dy = xyzq_j.y - xyzq_i.y;
              float dz = xyzq_j.z - xyzq_i.z;              

              float r2 = dx*dx + dy*dy + dz*dz;

              if (r2 < cutoff2) {
                calcForceEnergy<doEnergy, doSlow>(
                  r2, xyzq_i.w, xyzq_j.w, dx, dy, dz,
                  vdwtypei, s_vdwtypej[iwarp][(wid+t) & (WARPSIZE-1)], vdwCoefTable, 
                  vdwCoefTableTex, forceTableTex, energyTableTex,
                  iforce, iforceSlow,
                  s_jforce[iwarp][(wid+t) & (WARPSIZE-1)],
                  s_jforceSlow[iwarp][(wid+t) & (WARPSIZE-1)],
                  energyVdw, energyElec, energySlow);
              } // (r2 < cutoff2)
            } // (excl & 1)
            excl >>= 1;
            WARP_SYNC(WARP_FULL_MASK);
          } // t
        }

        // Write j-forces
        storeForces<doSlow>(jatomStart + wid, s_jforce[iwarp][wid], s_jforceSlow[iwarp][wid],
          devForce_x, devForce_y, devForce_z,
          devForceSlow_x, devForceSlow_y, devForceSlow_z);

        // Write exclusions
        if (doPairlist && WARP_ANY(WARP_FULL_MASK, nexcluded & 1)) {
          int anyexcl = (65536 | WARP_ANY(WARP_FULL_MASK, excl));
          // Mark this jtile as non-empty:
          //  VdW:      1 if tile has atom pairs within pairlist cutoff and some these atoms interact
          //  GBIS: 65536 if tile has atom pairs within pairlist cutoff but not necessary interacting (i.e. these atoms are fixed or excluded)
          if (wid == 0) jtiles[jtile] = anyexcl;
          // Store exclusions
          tileExcls[jtile].excl[wid] = excl;
          // itileListLen:
          // lower 16 bits number of tiles with atom pairs within pairlist cutoff that interact
          // upper 16 bits number of tiles with atom pairs within pairlist cutoff (but not necessary interacting)
          itileListLen += anyexcl;
          // NOTE, this minJatomStart is only stored once for the first tile list entry
          // minJatomStart = min(minJatomStart, jatomStart);
        }

      } // jtile

      // Write i-forces
      storeForces<doSlow>(iatomStart + wid, iforce, iforceSlow,
                          devForce_x, devForce_y, devForce_z,
                          devForceSlow_x, devForceSlow_y, devForceSlow_z);
    }
    // Done with computation

    // Save pairlist stuff
    if (doPairlist) {

      // Warp index (0...warpsize-1)
      const int wid = threadIdx.x % WARPSIZE;

      if (wid == 0) {
        // minJatomStart is in range [0 ... atomStorageSize-1]
        //int atom0 = (minJatomStart)/WARPSIZE;
        // int atom0 = 0;
        // int storageOffset = atomStorageSize/WARPSIZE;
        // int itileListLen = 0;
        // for (int jtile=jtileStart;jtile <= jtileEnd;jtile++) itileListLen += jtiles[jtile];
        // Store 0 if itileListLen == 0
        // tileListDepth[itileList] = (itileListLen > 0)*(itileListLen*storageOffset + atom0);
        tileListDepth[itileList] = itileListLen;
        tileListOrder[itileList] = itileList;
        // Number of active tilelists with tile with atom pairs within pairlist cutoff that interact
        if ((itileListLen & 65535) > 0) atomicAdd(&tileListStat->numTileLists, 1);
        // Number of active tilelists with tiles with atom pairs within pairlist cutoff (but not necessary interacting)
        if (itileListLen > 0) atomicAdd(&tileListStat->numTileListsGBIS, 1);
        // NOTE: always numTileListsGBIS >= numTileLists
      }

      typedef cub::WarpReduce<int> WarpReduceInt;
      __shared__ typename WarpReduceInt::TempStorage tempStorage[NONBONDKERNEL_NUM_WARP];
      int warpId = threadIdx.x / WARPSIZE;
      // Remove indicator bit
      nexcluded >>= 1;
      volatile int nexcludedWarp = WarpReduceInt(tempStorage[warpId]).Sum(nexcluded);
      if (wid == 0) atomicAdd(&tileListStat->numExcluded, nexcludedWarp);

    }

    if (doVirial) {
      // Warp index (0...warpsize-1)
      const int wid = threadIdx.x % WARPSIZE;

      typedef cub::WarpReduce<float> WarpReduce;
      __shared__ typename WarpReduce::TempStorage tempStorage[NONBONDKERNEL_NUM_WARP];
      int warpId = threadIdx.x / WARPSIZE;
      volatile float iforcexSum = WarpReduce(tempStorage[warpId]).Sum(iforce.x);
      WARP_SYNC(WARP_FULL_MASK);
      volatile float iforceySum = WarpReduce(tempStorage[warpId]).Sum(iforce.y);
      WARP_SYNC(WARP_FULL_MASK);
      volatile float iforcezSum = WarpReduce(tempStorage[warpId]).Sum(iforce.z);
      WARP_SYNC(WARP_FULL_MASK);
      if (wid == 0) {
        virialEnergy[itileList].forcex = iforcexSum;
        virialEnergy[itileList].forcey = iforceySum;
        virialEnergy[itileList].forcez = iforcezSum;
      }

      if (doSlow) {
        iforcexSum = WarpReduce(tempStorage[warpId]).Sum(iforceSlow.x);
        WARP_SYNC(WARP_FULL_MASK);
        iforceySum = WarpReduce(tempStorage[warpId]).Sum(iforceSlow.y);
        WARP_SYNC(WARP_FULL_MASK);
        iforcezSum = WarpReduce(tempStorage[warpId]).Sum(iforceSlow.z);
        WARP_SYNC(WARP_FULL_MASK);
        if (wid == 0) {
          virialEnergy[itileList].forceSlowx = iforcexSum;
          virialEnergy[itileList].forceSlowy = iforceySum;
          virialEnergy[itileList].forceSlowz = iforcezSum;
        }
      }
    }

    // Reduce energy
    if (doEnergy) {
      // NOTE: We must hand write these warp-wide reductions to avoid excess register spillage
      //       (Why does CUB suck here?)
#pragma unroll
      for (int i=16;i >= 1;i/=2) {
        energyVdw += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energyVdw, i, 32);
        energyElec += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energyElec, i, 32);
        if (doSlow) energySlow += WARP_SHUFFLE_XOR(WARP_FULL_MASK, energySlow, i, 32);
      }

      if (threadIdx.x % WARPSIZE == 0) {
        virialEnergy[itileList].energyVdw  = energyVdw;
        virialEnergy[itileList].energyElec = energyElec;
        if (doSlow) virialEnergy[itileList].energySlow = energySlow;
      }
    }

    if (doStreaming) {
      // Make sure devForces and devForcesSlow have been written into device memory
      WARP_SYNC(WARP_FULL_MASK);
      __threadfence();

      int patchDone[2] = {false, false};
      const int wid = threadIdx.x % WARPSIZE;
      if (wid == 0) {
        int patchCountOld0 = atomicInc(&patchNumCount[patchInd.x], (unsigned int)(patchNumList.x-1));
        patchDone[0] = (patchCountOld0 + 1 == patchNumList.x);
        if (patchInd.x != patchInd.y) {
          int patchCountOld1 = atomicInc(&patchNumCount[patchInd.y], (unsigned int)(patchNumList.y-1));
          patchDone[1] = (patchCountOld1 + 1 == patchNumList.y);
        }
      }

      patchDone[0] = WARP_ANY(WARP_FULL_MASK, patchDone[0]);
      patchDone[1] = WARP_ANY(WARP_FULL_MASK, patchDone[1]);

      if (patchDone[0]) {
        // Patch 1 is done, write onto host-mapped memory
        CudaPatchRecord patch = cudaPatches[patchInd.x];
        int start = patch.atomStart;
        int end   = start + patch.numAtoms;
        for (int i=start+wid;i < end;i+=WARPSIZE) {
          mapForces[i] = make_float4(devForce_x[i],
            devForce_y[i], devForce_z[i], devForce_w[i]);
          if (doSlow){
            mapForcesSlow[i] = make_float4(devForceSlow_x[i],
                                           devForceSlow_y[i], 
                                           devForceSlow_z[i], 
                                           devForceSlow_w[i]);
          }
        }
      }
      if (patchDone[1]) {
        // Patch 2 is done
        CudaPatchRecord patch = cudaPatches[patchInd.y];
        int start = patch.atomStart;
        int end   = start + patch.numAtoms;
        for (int i=start+wid;i < end;i+=WARPSIZE) {
          mapForces[i] = make_float4(devForce_x[i],  devForce_y[i], devForce_z[i], devForce_w[i]);
          if (doSlow){
            mapForcesSlow[i] = make_float4(devForceSlow_x[i],
                                          devForceSlow_y[i], 
                                          devForceSlow_z[i], 
                                          devForceSlow_w[i]);
          }
        }
      }

      if (patchDone[0] || patchDone[1]) {
        // Make sure mapForces and mapForcesSlow are up-to-date
        WARP_SYNC(WARP_FULL_MASK);
        __threadfence_system();
        // Add patch into "patchReadyQueue"
        if (wid == 0) {
          if (patchDone[0]) {
            int ind = atomicAdd(&tileListStat->patchReadyQueueCount, 1);
            // int ind = atomicInc((unsigned int *)&mapPatchReadyQueue[numPatches], numPatches-1);
            mapPatchReadyQueue[ind] = patchInd.x;
          }
          if (patchDone[1]) {
            int ind = atomicAdd(&tileListStat->patchReadyQueueCount, 1);
            // int ind = atomicInc((unsigned int *)&mapPatchReadyQueue[numPatches], numPatches-1);
            mapPatchReadyQueue[ind] = patchInd.y;
          }
        }
      }
    }

    if (doStreaming && outputOrder != NULL && threadIdx.x % WARPSIZE == 0) {
      int index = atomicAdd(&tileListStat->outputOrderIndex, 1);
      outputOrder[index] = itileList;
    }
  } // if (itileList < numTileLists)
}

//
// Finish up - reduce virials from nonbonded kernel
//
#define REDUCENONBONDEDVIRIALKERNEL_NUM_WARP 32
__global__ void reduceNonbondedVirialKernel(const bool doSlow,
  const int atomStorageSize,
  const float4* __restrict__ xyzq,
  const float4* __restrict__ devForces, const float4* __restrict__ devForcesSlow,
  VirialEnergy* __restrict__ virialEnergy) {

  for (int ibase = blockIdx.x*blockDim.x;ibase < atomStorageSize;ibase += blockDim.x*gridDim.x)
  {
    int i = ibase + threadIdx.x;

    // Set to zero to avoid nan*0
    float4 pos;
    pos.x = 0.0f;
    pos.y = 0.0f;
    pos.z = 0.0f;
    float4 force, forceSlow;
    force.x = 0.0f;
    force.y = 0.0f;
    force.z = 0.0f;
    forceSlow.x = 0.0f;
    forceSlow.y = 0.0f;
    forceSlow.z = 0.0f;
    if (i < atomStorageSize) {
      pos = xyzq[i];
      force = devForces[i];
      if (doSlow) forceSlow = devForcesSlow[i];
    }
    // Reduce across the entire thread block
    float vxxt = force.x*pos.x;
    float vxyt = force.x*pos.y;
    float vxzt = force.x*pos.z;
    float vyxt = force.y*pos.x;
    float vyyt = force.y*pos.y;
    float vyzt = force.y*pos.z;
    float vzxt = force.z*pos.x;
    float vzyt = force.z*pos.y;
    float vzzt = force.z*pos.z;
    // atomicAdd(&virialEnergy->virial[0], (double)vxx);
    // atomicAdd(&virialEnergy->virial[1], (double)vxy);
    // atomicAdd(&virialEnergy->virial[2], (double)vxz);
    // atomicAdd(&virialEnergy->virial[3], (double)vyx);
    // atomicAdd(&virialEnergy->virial[4], (double)vyy);
    // atomicAdd(&virialEnergy->virial[5], (double)vyz);
    // atomicAdd(&virialEnergy->virial[6], (double)vzx);
    // atomicAdd(&virialEnergy->virial[7], (double)vzy);
    // atomicAdd(&virialEnergy->virial[8], (double)vzz);

    typedef cub::BlockReduce<float, REDUCENONBONDEDVIRIALKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
    __shared__ typename BlockReduce::TempStorage tempStorage;
    volatile float vxx = BlockReduce(tempStorage).Sum(vxxt); BLOCK_SYNC;
    volatile float vxy = BlockReduce(tempStorage).Sum(vxyt); BLOCK_SYNC;
    volatile float vxz = BlockReduce(tempStorage).Sum(vxzt); BLOCK_SYNC;
    volatile float vyx = BlockReduce(tempStorage).Sum(vyxt); BLOCK_SYNC;
    volatile float vyy = BlockReduce(tempStorage).Sum(vyyt); BLOCK_SYNC;
    volatile float vyz = BlockReduce(tempStorage).Sum(vyzt); BLOCK_SYNC;
    volatile float vzx = BlockReduce(tempStorage).Sum(vzxt); BLOCK_SYNC;
    volatile float vzy = BlockReduce(tempStorage).Sum(vzyt); BLOCK_SYNC;
    volatile float vzz = BlockReduce(tempStorage).Sum(vzzt); BLOCK_SYNC;
    if (threadIdx.x == 0) {
      atomicAdd(&virialEnergy->virial[0], (double)vxx);
      atomicAdd(&virialEnergy->virial[1], (double)vxy);
      atomicAdd(&virialEnergy->virial[2], (double)vxz);
      atomicAdd(&virialEnergy->virial[3], (double)vyx);
      atomicAdd(&virialEnergy->virial[4], (double)vyy);
      atomicAdd(&virialEnergy->virial[5], (double)vyz);
      atomicAdd(&virialEnergy->virial[6], (double)vzx);
      atomicAdd(&virialEnergy->virial[7], (double)vzy);
      atomicAdd(&virialEnergy->virial[8], (double)vzz);
    }

    if (doSlow) {
      // if (isnan(forceSlow.x) || isnan(forceSlow.y) || isnan(forceSlow.z))
      float vxxSlowt = forceSlow.x*pos.x;
      float vxySlowt = forceSlow.x*pos.y;
      float vxzSlowt = forceSlow.x*pos.z;
      float vyxSlowt = forceSlow.y*pos.x;
      float vyySlowt = forceSlow.y*pos.y;
      float vyzSlowt = forceSlow.y*pos.z;
      float vzxSlowt = forceSlow.z*pos.x;
      float vzySlowt = forceSlow.z*pos.y;
      float vzzSlowt = forceSlow.z*pos.z;
      // atomicAdd(&virialEnergy->virialSlow[0], (double)vxxSlow);
      // atomicAdd(&virialEnergy->virialSlow[1], (double)vxySlow);
      // atomicAdd(&virialEnergy->virialSlow[2], (double)vxzSlow);
      // atomicAdd(&virialEnergy->virialSlow[3], (double)vyxSlow);
      // atomicAdd(&virialEnergy->virialSlow[4], (double)vyySlow);
      // atomicAdd(&virialEnergy->virialSlow[5], (double)vyzSlow);
      // atomicAdd(&virialEnergy->virialSlow[6], (double)vzxSlow);
      // atomicAdd(&virialEnergy->virialSlow[7], (double)vzySlow);
      // atomicAdd(&virialEnergy->virialSlow[8], (double)vzzSlow);
      volatile float vxxSlow = BlockReduce(tempStorage).Sum(vxxSlowt); BLOCK_SYNC;
      volatile float vxySlow = BlockReduce(tempStorage).Sum(vxySlowt); BLOCK_SYNC;
      volatile float vxzSlow = BlockReduce(tempStorage).Sum(vxzSlowt); BLOCK_SYNC;
      volatile float vyxSlow = BlockReduce(tempStorage).Sum(vyxSlowt); BLOCK_SYNC;
      volatile float vyySlow = BlockReduce(tempStorage).Sum(vyySlowt); BLOCK_SYNC;
      volatile float vyzSlow = BlockReduce(tempStorage).Sum(vyzSlowt); BLOCK_SYNC;
      volatile float vzxSlow = BlockReduce(tempStorage).Sum(vzxSlowt); BLOCK_SYNC;
      volatile float vzySlow = BlockReduce(tempStorage).Sum(vzySlowt); BLOCK_SYNC;
      volatile float vzzSlow = BlockReduce(tempStorage).Sum(vzzSlowt); BLOCK_SYNC;
      if (threadIdx.x == 0) {
        atomicAdd(&virialEnergy->virialSlow[0], (double)vxxSlow);
        atomicAdd(&virialEnergy->virialSlow[1], (double)vxySlow);
        atomicAdd(&virialEnergy->virialSlow[2], (double)vxzSlow);
        atomicAdd(&virialEnergy->virialSlow[3], (double)vyxSlow);
        atomicAdd(&virialEnergy->virialSlow[4], (double)vyySlow);
        atomicAdd(&virialEnergy->virialSlow[5], (double)vyzSlow);
        atomicAdd(&virialEnergy->virialSlow[6], (double)vzxSlow);
        atomicAdd(&virialEnergy->virialSlow[7], (double)vzySlow);
        atomicAdd(&virialEnergy->virialSlow[8], (double)vzzSlow);
      }
    }
  
  }

}

#define REDUCEVIRIALENERGYKERNEL_NUM_WARP 32
__global__ void reduceVirialEnergyKernel(
  const bool doEnergy, const bool doVirial, const bool doSlow,
  const int numTileLists,
  const TileListVirialEnergy* __restrict__ tileListVirialEnergy,
  VirialEnergy* __restrict__ virialEnergy) {

  for (int ibase = blockIdx.x*blockDim.x;ibase < numTileLists;ibase += blockDim.x*gridDim.x)
  {
    int itileList = ibase + threadIdx.x;
    TileListVirialEnergy ve;
    if (itileList < numTileLists) {
      ve = tileListVirialEnergy[itileList];
    } else {
      // Set to zero to avoid nan*0
      if (doVirial) {
        ve.shx = 0.0f;
        ve.shy = 0.0f;
        ve.shz = 0.0f;
        ve.forcex = 0.0f;
        ve.forcey = 0.0f;
        ve.forcez = 0.0f;
        ve.forceSlowx = 0.0f;
        ve.forceSlowy = 0.0f;
        ve.forceSlowz = 0.0f;
      }
      if (doEnergy) {
        ve.energyVdw = 0.0;
        ve.energyElec = 0.0;
        ve.energySlow = 0.0;
        // ve.energyGBIS = 0.0;
      }
    }

    if (doVirial) {
      typedef cub::BlockReduce<float, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
      __shared__ typename BlockReduce::TempStorage tempStorage;
      float vxxt = ve.forcex*ve.shx;
      float vxyt = ve.forcex*ve.shy;
      float vxzt = ve.forcex*ve.shz;
      float vyxt = ve.forcey*ve.shx;
      float vyyt = ve.forcey*ve.shy;
      float vyzt = ve.forcey*ve.shz;
      float vzxt = ve.forcez*ve.shx;
      float vzyt = ve.forcez*ve.shy;
      float vzzt = ve.forcez*ve.shz;
      volatile float vxx = BlockReduce(tempStorage).Sum(vxxt); BLOCK_SYNC;
      volatile float vxy = BlockReduce(tempStorage).Sum(vxyt); BLOCK_SYNC;
      volatile float vxz = BlockReduce(tempStorage).Sum(vxzt); BLOCK_SYNC;
      volatile float vyx = BlockReduce(tempStorage).Sum(vyxt); BLOCK_SYNC;
      volatile float vyy = BlockReduce(tempStorage).Sum(vyyt); BLOCK_SYNC;
      volatile float vyz = BlockReduce(tempStorage).Sum(vyzt); BLOCK_SYNC;
      volatile float vzx = BlockReduce(tempStorage).Sum(vzxt); BLOCK_SYNC;
      volatile float vzy = BlockReduce(tempStorage).Sum(vzyt); BLOCK_SYNC;
      volatile float vzz = BlockReduce(tempStorage).Sum(vzzt); BLOCK_SYNC;
      if (threadIdx.x == 0) {
        atomicAdd(&virialEnergy->virial[0], (double)vxx);
        atomicAdd(&virialEnergy->virial[1], (double)vxy);
        atomicAdd(&virialEnergy->virial[2], (double)vxz);
        atomicAdd(&virialEnergy->virial[3], (double)vyx);
        atomicAdd(&virialEnergy->virial[4], (double)vyy);
        atomicAdd(&virialEnergy->virial[5], (double)vyz);
        atomicAdd(&virialEnergy->virial[6], (double)vzx);
        atomicAdd(&virialEnergy->virial[7], (double)vzy);
        atomicAdd(&virialEnergy->virial[8], (double)vzz);
      }

      if (doSlow) {
        typedef cub::BlockReduce<float, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
        __shared__ typename BlockReduce::TempStorage tempStorage;
        float vxxt = ve.forceSlowx*ve.shx;
        float vxyt = ve.forceSlowx*ve.shy;
        float vxzt = ve.forceSlowx*ve.shz;
        float vyxt = ve.forceSlowy*ve.shx;
        float vyyt = ve.forceSlowy*ve.shy;
        float vyzt = ve.forceSlowy*ve.shz;
        float vzxt = ve.forceSlowz*ve.shx;
        float vzyt = ve.forceSlowz*ve.shy;
        float vzzt = ve.forceSlowz*ve.shz;
        volatile float vxx = BlockReduce(tempStorage).Sum(vxxt); BLOCK_SYNC;
        volatile float vxy = BlockReduce(tempStorage).Sum(vxyt); BLOCK_SYNC;
        volatile float vxz = BlockReduce(tempStorage).Sum(vxzt); BLOCK_SYNC;
        volatile float vyx = BlockReduce(tempStorage).Sum(vyxt); BLOCK_SYNC;
        volatile float vyy = BlockReduce(tempStorage).Sum(vyyt); BLOCK_SYNC;
        volatile float vyz = BlockReduce(tempStorage).Sum(vyzt); BLOCK_SYNC;
        volatile float vzx = BlockReduce(tempStorage).Sum(vzxt); BLOCK_SYNC;
        volatile float vzy = BlockReduce(tempStorage).Sum(vzyt); BLOCK_SYNC;
        volatile float vzz = BlockReduce(tempStorage).Sum(vzzt); BLOCK_SYNC;
        if (threadIdx.x == 0) {
          atomicAdd(&virialEnergy->virialSlow[0], (double)vxx);
          atomicAdd(&virialEnergy->virialSlow[1], (double)vxy);
          atomicAdd(&virialEnergy->virialSlow[2], (double)vxz);
          atomicAdd(&virialEnergy->virialSlow[3], (double)vyx);
          atomicAdd(&virialEnergy->virialSlow[4], (double)vyy);
          atomicAdd(&virialEnergy->virialSlow[5], (double)vyz);
          atomicAdd(&virialEnergy->virialSlow[6], (double)vzx);
          atomicAdd(&virialEnergy->virialSlow[7], (double)vzy);
          atomicAdd(&virialEnergy->virialSlow[8], (double)vzz);
        }
      }
    }

    if (doEnergy) {
      typedef cub::BlockReduce<double, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
      __shared__ typename BlockReduce::TempStorage tempStorage;
      volatile double energyVdw  = BlockReduce(tempStorage).Sum(ve.energyVdw); BLOCK_SYNC;
      volatile double energyElec = BlockReduce(tempStorage).Sum(ve.energyElec); BLOCK_SYNC;
      if (threadIdx.x == 0) {
          atomicAdd(&virialEnergy->energyVdw, (double)energyVdw);
          atomicAdd(&virialEnergy->energyElec, (double)energyElec);
      }
      if (doSlow) {
        volatile double energySlow = BlockReduce(tempStorage).Sum(ve.energySlow); BLOCK_SYNC;
        if (threadIdx.x == 0) atomicAdd(&virialEnergy->energySlow, (double)energySlow);
      }
      // if (doGBIS) {
      //   double energyGBIS = BlockReduce(tempStorage).Sum(ve.energyGBIS); BLOCK_SYNC;
      //   if (threadIdx.x == 0) atomicAdd(&virialEnergy->energyGBIS, (double)energyGBIS);
      // }
    }

  }

}

#define REDUCEGBISENERGYKERNEL_NUM_WARP 32
__global__ void reduceGBISEnergyKernel(const int numTileLists,
  const TileListVirialEnergy* __restrict__ tileListVirialEnergy,
  VirialEnergy* __restrict__ virialEnergy) {

  for (int ibase = blockIdx.x*blockDim.x;ibase < numTileLists;ibase += blockDim.x*gridDim.x)
  {
    int itileList = ibase + threadIdx.x;
    double energyGBISt = 0.0;
    if (itileList < numTileLists) {
      energyGBISt = tileListVirialEnergy[itileList].energyGBIS;
    }

    typedef cub::BlockReduce<double, REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE> BlockReduce;
    __shared__ typename BlockReduce::TempStorage tempStorage;
    volatile double energyGBIS = BlockReduce(tempStorage).Sum(energyGBISt); BLOCK_SYNC;
    if (threadIdx.x == 0) atomicAdd(&virialEnergy->energyGBIS, (double)energyGBIS);
  }

}

// ##############################################################################################
// ##############################################################################################
// ##############################################################################################

CudaComputeNonbondedKernel::CudaComputeNonbondedKernel(int deviceID, CudaNonbondedTables& cudaNonbondedTables,
  bool doStreaming) : deviceID(deviceID), cudaNonbondedTables(cudaNonbondedTables), doStreaming(doStreaming) {
  
  cudaCheck(cudaSetDevice(deviceID));

  overflowExclusions = NULL;
  overflowExclusionsSize = 0;

  exclIndexMaxDiff = NULL;
  exclIndexMaxDiffSize = 0;

  atomIndex = NULL;
  atomIndexSize = 0;

  vdwTypes = NULL;
  vdwTypesSize = 0;

  patchNumCount = NULL;
  patchNumCountSize = 0;

  patchReadyQueue = NULL;
  patchReadyQueueSize = 0;

  force_x = force_y = force_z = force_w = NULL;
  forceSize = 0;
  forceSlow_x = forceSlow_y = forceSlow_z = forceSlow_w = NULL;
  forceSlowSize = 0;
}

void CudaComputeNonbondedKernel::reallocate_forceSOA(int atomStorageSize)
{
  reallocate_device<float>(&force_x, &forceSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&force_y, &forceSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&force_z, &forceSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&force_w, &forceSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&forceSlow_x, &forceSlowSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&forceSlow_y, &forceSlowSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&forceSlow_z, &forceSlowSize, atomStorageSize, 1.4f);
  reallocate_device<float>(&forceSlow_w, &forceSlowSize, atomStorageSize, 1.4f);  
}

CudaComputeNonbondedKernel::~CudaComputeNonbondedKernel() {
  cudaCheck(cudaSetDevice(deviceID));
  if (overflowExclusions != NULL) deallocate_device<unsigned int>(&overflowExclusions);
  if (exclIndexMaxDiff != NULL) deallocate_device<int2>(&exclIndexMaxDiff);
  if (atomIndex != NULL) deallocate_device<int>(&atomIndex);
  if (vdwTypes != NULL) deallocate_device<int>(&vdwTypes);
  if (patchNumCount != NULL) deallocate_device<unsigned int>(&patchNumCount);
  if (patchReadyQueue != NULL) deallocate_host<int>(&patchReadyQueue);
  if (force_x != NULL) deallocate_device<float>(&force_x);
  if (force_y != NULL) deallocate_device<float>(&force_y);
  if (force_z != NULL) deallocate_device<float>(&force_z);
  if (force_w != NULL) deallocate_device<float>(&force_w);
  if (forceSlow_x != NULL) deallocate_device<float>(&forceSlow_x);
  if (forceSlow_y != NULL) deallocate_device<float>(&forceSlow_y);
  if (forceSlow_z != NULL) deallocate_device<float>(&forceSlow_z);
  if (forceSlow_w != NULL) deallocate_device<float>(&forceSlow_w);  
}

void CudaComputeNonbondedKernel::updateVdwTypesExcl(const int atomStorageSize, const int* h_vdwTypes,
  const int2* h_exclIndexMaxDiff, const int* h_atomIndex, cudaStream_t stream) {

  reallocate_device<int>(&vdwTypes, &vdwTypesSize, atomStorageSize, OVERALLOC);
  reallocate_device<int2>(&exclIndexMaxDiff, &exclIndexMaxDiffSize, atomStorageSize, OVERALLOC);
  reallocate_device<int>(&atomIndex, &atomIndexSize, atomStorageSize, OVERALLOC);

  copy_HtoD<int>(h_vdwTypes, vdwTypes, atomStorageSize, stream);
  copy_HtoD<int2>(h_exclIndexMaxDiff, exclIndexMaxDiff, atomStorageSize, stream);
  copy_HtoD<int>(h_atomIndex, atomIndex, atomStorageSize, stream);
}

int* CudaComputeNonbondedKernel::getPatchReadyQueue() {
  if (!doStreaming) {
    NAMD_die("CudaComputeNonbondedKernel::getPatchReadyQueue() called on non-streaming kernel");
  }
  return patchReadyQueue;
}

template <int doSlow>
__global__ void transposeForcesKernel(float4 *f, float4 *fSlow,
                                      float *fx, float *fy, float *fz, float *fw,
                                      float *fSlowx, float *fSlowy, float *fSlowz, float *fSloww,
                                      int n)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if (tid < n) {
    f[tid] = make_float4(fx[tid], fy[tid], fz[tid], fw[tid]);
    if (doSlow) {
      fSlow[tid] = make_float4(fSlowx[tid], fSlowy[tid], fSlowz[tid], fSloww[tid]);
    }
  }
}



void CudaComputeNonbondedKernel::nonbondedForce(CudaTileListKernel& tlKernel,
  const int atomStorageSize, const bool doPairlist,
  const bool doEnergy, const bool doVirial, const bool doSlow,
  const float3 lata, const float3 latb, const float3 latc,
  const float4* h_xyzq, const float cutoff2, 
  float4* d_forces, float4* d_forcesSlow,
  float4* h_forces, float4* h_forcesSlow,
  cudaStream_t stream) {

  if (!doPairlist) copy_HtoD<float4>(h_xyzq, tlKernel.get_xyzq(), atomStorageSize, stream);

  // clear_device_array<float4>(d_forces, atomStorageSize, stream);
  // if (doSlow) clear_device_array<float4>(d_forcesSlow, atomStorageSize, stream);

  
  // XXX TODO: Clear all of these
  if(1){
     // two clears
     tlKernel.clearTileListStat(stream);
     clear_device_array<float>(force_x, atomStorageSize, stream);
     clear_device_array<float>(force_y, atomStorageSize, stream);
     clear_device_array<float>(force_z, atomStorageSize, stream);
     clear_device_array<float>(force_w, atomStorageSize, stream);
     if (doSlow) {
       clear_device_array<float>(forceSlow_x, atomStorageSize, stream);
       clear_device_array<float>(forceSlow_y, atomStorageSize, stream);
       clear_device_array<float>(forceSlow_z, atomStorageSize, stream);
       clear_device_array<float>(forceSlow_w, atomStorageSize, stream);
     }
  }

  // --- streaming ----
  float4* m_forces = NULL;
  float4* m_forcesSlow = NULL;
  int* m_patchReadyQueue = NULL;
  int numPatches = 0;
  unsigned int* patchNumCountPtr = NULL;
  if (doStreaming) {
    numPatches = tlKernel.getNumPatches();
    if (reallocate_device<unsigned int>(&patchNumCount, &patchNumCountSize, numPatches)) {
      // If re-allocated, clear array
      clear_device_array<unsigned int>(patchNumCount, numPatches, stream);
    }
    patchNumCountPtr = patchNumCount;
    bool re = reallocate_host<int>(&patchReadyQueue, &patchReadyQueueSize, numPatches, cudaHostAllocMapped);
    if (re) {
      // If re-allocated, re-set to "-1"
      for (int i=0;i < numPatches;i++) patchReadyQueue[i] = -1;
    }
    cudaCheck(cudaHostGetDevicePointer(&m_patchReadyQueue, patchReadyQueue, 0));
    cudaCheck(cudaHostGetDevicePointer(&m_forces, h_forces, 0));
    cudaCheck(cudaHostGetDevicePointer(&m_forcesSlow, h_forcesSlow, 0));
  }
  // -----------------

  if (doVirial || doEnergy) {
    tlKernel.setTileListVirialEnergyLength(tlKernel.getNumTileLists());
  }

  int shMemSize = 0;

  int* outputOrderPtr = tlKernel.getOutputOrder();

  int nwarp = NONBONDKERNEL_NUM_WARP;
  int nthread = WARPSIZE*nwarp;
  int start = 0;
  while (start < tlKernel.getNumTileLists())
  {

    int nleft = tlKernel.getNumTileLists() - start;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (nleft-1)/nwarp+1);

#define CALL(DOENERGY, DOVIRIAL, DOSLOW, DOPAIRLIST, DOSTREAMING) \
    nonbondedForceKernel<DOENERGY, DOVIRIAL, DOSLOW, DOPAIRLIST, DOSTREAMING> \
  <<< nblock, nthread, shMemSize, stream >>>  \
  (start, tlKernel.getNumTileLists(), tlKernel.getTileLists(), tlKernel.getTileExcls(), tlKernel.getTileJatomStart(), \
    cudaNonbondedTables.getVdwCoefTableWidth(), cudaNonbondedTables.getVdwCoefTable(), \
    vdwTypes, lata, latb, latc, tlKernel.get_xyzq(), cutoff2, \
    cudaNonbondedTables.getVdwCoefTableTex(), cudaNonbondedTables.getForceTableTex(), cudaNonbondedTables.getEnergyTableTex(), \
    atomStorageSize, tlKernel.get_plcutoff2(), tlKernel.getPatchPairs(), atomIndex, exclIndexMaxDiff, overflowExclusions, \
    tlKernel.getTileListDepth(), tlKernel.getTileListOrder(), tlKernel.getJtiles(), tlKernel.getTileListStatDevPtr(), \
    tlKernel.getBoundingBoxes(), d_forces, d_forcesSlow, \
    force_x, force_y, force_z, force_w, \
    forceSlow_x, forceSlow_y, forceSlow_z, forceSlow_w, \
    numPatches, patchNumCountPtr, tlKernel.getCudaPatches(), m_forces, m_forcesSlow, m_patchReadyQueue, \
    outputOrderPtr, tlKernel.getTileListVirialEnergy()); called=true

    bool called = false;

    if (doStreaming) {
      if (!doEnergy && !doVirial && !doSlow && !doPairlist) CALL(0, 0, 0, 0, 1);
      if (!doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(0, 0, 1, 0, 1);
      if (!doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(0, 1, 0, 0, 1);
      if (!doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(0, 1, 1, 0, 1);
      if ( doEnergy && !doVirial && !doSlow && !doPairlist) CALL(1, 0, 0, 0, 1);
      if ( doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(1, 0, 1, 0, 1);
      if ( doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(1, 1, 0, 0, 1);
      if ( doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(1, 1, 1, 0, 1);

      if (!doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(0, 0, 0, 1, 1);
      if (!doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(0, 0, 1, 1, 1);
      if (!doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(0, 1, 0, 1, 1);
      if (!doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(0, 1, 1, 1, 1);
      if ( doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(1, 0, 0, 1, 1);
      if ( doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(1, 0, 1, 1, 1);
      if ( doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(1, 1, 0, 1, 1);
      if ( doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(1, 1, 1, 1, 1);
    } else {
      if (!doEnergy && !doVirial && !doSlow && !doPairlist) CALL(0, 0, 0, 0, 0);
      if (!doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(0, 0, 1, 0, 0);
      if (!doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(0, 1, 0, 0, 0);
      if (!doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(0, 1, 1, 0, 0);
      if ( doEnergy && !doVirial && !doSlow && !doPairlist) CALL(1, 0, 0, 0, 0);
      if ( doEnergy && !doVirial &&  doSlow && !doPairlist) CALL(1, 0, 1, 0, 0);
      if ( doEnergy &&  doVirial && !doSlow && !doPairlist) CALL(1, 1, 0, 0, 0);
      if ( doEnergy &&  doVirial &&  doSlow && !doPairlist) CALL(1, 1, 1, 0, 0);

      if (!doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(0, 0, 0, 1, 0);
      if (!doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(0, 0, 1, 1, 0);
      if (!doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(0, 1, 0, 1, 0);
      if (!doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(0, 1, 1, 1, 0);
      if ( doEnergy && !doVirial && !doSlow &&  doPairlist) CALL(1, 0, 0, 1, 0);
      if ( doEnergy && !doVirial &&  doSlow &&  doPairlist) CALL(1, 0, 1, 1, 0);
      if ( doEnergy &&  doVirial && !doSlow &&  doPairlist) CALL(1, 1, 0, 1, 0);
      if ( doEnergy &&  doVirial &&  doSlow &&  doPairlist) CALL(1, 1, 1, 1, 0);
    }

    if (!called) {
      NAMD_die("CudaComputeNonbondedKernel::nonbondedForce, none of the kernels called");
    }

    {
      int block = 128;
      int grid = (atomStorageSize + block - 1)/block;
      if (doSlow) 
        transposeForcesKernel<1><<<grid, block, 0, stream>>>(d_forces, d_forcesSlow,
                       force_x, force_y, force_z, force_w,
                       forceSlow_x, forceSlow_y, forceSlow_z, forceSlow_w,
                       atomStorageSize);
      else
        transposeForcesKernel<0><<<grid, block, 0, stream>>>(d_forces, d_forcesSlow,
                       force_x, force_y, force_z, force_w,
                       forceSlow_x, forceSlow_y, forceSlow_z, forceSlow_w,
                       atomStorageSize);        
    }

#undef CALL
    cudaCheck(cudaGetLastError());

    start += nblock*nwarp;
  }

}

//
// Perform virial and energy reductions for non-bonded force calculation
//
void CudaComputeNonbondedKernel::reduceVirialEnergy(CudaTileListKernel& tlKernel,
  const int atomStorageSize, const bool doEnergy, const bool doVirial, const bool doSlow, const bool doGBIS,
  float4* d_forces, float4* d_forcesSlow,
  VirialEnergy* d_virialEnergy, cudaStream_t stream) {

  if (doEnergy || doVirial) {
    clear_device_array<VirialEnergy>(d_virialEnergy, 1, stream);
  }

  if (doVirial)
  {
    int nthread = REDUCENONBONDEDVIRIALKERNEL_NUM_WARP*WARPSIZE;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (atomStorageSize-1)/nthread+1);
    reduceNonbondedVirialKernel <<< nblock, nthread, 0, stream >>>
    (doSlow, atomStorageSize, tlKernel.get_xyzq(), d_forces, d_forcesSlow, d_virialEnergy);
    cudaCheck(cudaGetLastError());
  }

  if (doVirial || doEnergy)
  {
    int nthread = REDUCEVIRIALENERGYKERNEL_NUM_WARP*WARPSIZE;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getTileListVirialEnergyLength()-1)/nthread+1);
    reduceVirialEnergyKernel <<< nblock, nthread, 0, stream >>>
    (doEnergy, doVirial, doSlow, tlKernel.getTileListVirialEnergyLength(), tlKernel.getTileListVirialEnergy(), d_virialEnergy);
    cudaCheck(cudaGetLastError());
  }  

  if (doGBIS && doEnergy)
  {
    int nthread = REDUCEGBISENERGYKERNEL_NUM_WARP*WARPSIZE;
    int nblock = min(deviceCUDA->getMaxNumBlocks(), (tlKernel.getTileListVirialEnergyGBISLength()-1)/nthread+1);
    reduceGBISEnergyKernel <<< nblock, nthread, 0, stream >>>
    (tlKernel.getTileListVirialEnergyGBISLength(), tlKernel.getTileListVirialEnergy(), d_virialEnergy);
    cudaCheck(cudaGetLastError());
  }

}

void CudaComputeNonbondedKernel::bindExclusions(int numExclusions, unsigned int* exclusion_bits) {
	int nconst = ( numExclusions < MAX_CONST_EXCLUSIONS ? numExclusions : MAX_CONST_EXCLUSIONS );
	cudaCheck(cudaMemcpyToSymbol(constExclusions, exclusion_bits, nconst*sizeof(unsigned int), 0));

  reallocate_device<unsigned int>(&overflowExclusions, &overflowExclusionsSize, numExclusions);
  copy_HtoD_sync<unsigned int>(exclusion_bits, overflowExclusions, numExclusions);
}
