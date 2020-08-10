#ifndef CUDACOMPUTENONBONDEDKERNEL_H
#define CUDACOMPUTENONBONDEDKERNEL_H
#include "CudaUtils.h"
#include "CudaTileListKernel.h"
#include "CudaNonbondedTables.h"
#ifdef NAMD_CUDA

class CudaComputeNonbondedKernel {
private:

  const int deviceID;
  CudaNonbondedTables& cudaNonbondedTables;
  const bool doStreaming;

  // Exclusions
  unsigned int* overflowExclusions;
  int overflowExclusionsSize;

  int2* exclIndexMaxDiff;
  int exclIndexMaxDiffSize;

  // Atom indices
  int* atomIndex;
  int atomIndexSize;

  // VdW types
  int* vdwTypes;
  int vdwTypesSize;

  unsigned int* patchNumCount;
  int patchNumCountSize;

  int* patchReadyQueue;
  int patchReadyQueueSize;

  float *force_x, *force_y, *force_z, *force_w;
  int forceSize;
  float *forceSlow_x, *forceSlow_y, *forceSlow_z, *forceSlow_w;
  int forceSlowSize;
public:
  CudaComputeNonbondedKernel(int deviceID, CudaNonbondedTables& cudaNonbondedTables, bool doStreaming);
  ~CudaComputeNonbondedKernel();

  void updateVdwTypesExcl(const int atomStorageSize, const int* h_vdwTypes,
    const int2* h_exclIndexMaxDiff, const int* h_atomIndex, cudaStream_t stream);

  void nonbondedForce(CudaTileListKernel& tlKernel,
    const int atomStorageSize, const bool doPairlist,
    const bool doEnergy, const bool doVirial, const bool doSlow,
    const float3 lata, const float3 latb, const float3 latc,
    const float4* h_xyzq, const float cutoff2, 
    float4* d_forces, float4* d_forcesSlow,
    float4* h_forces, float4* h_forcesSlow,
    cudaStream_t stream);

  void reduceVirialEnergy(CudaTileListKernel& tlKernel,
    const int atomStorageSize, const bool doEnergy, const bool doVirial, const bool doSlow, const bool doGBIS,
    float4* d_forces, float4* d_forcesSlow,
    VirialEnergy* d_virialEnergy, cudaStream_t stream);

  void getVirialEnergy(VirialEnergy* h_virialEnergy, cudaStream_t stream);

  void bindExclusions(int numExclusions, unsigned int* exclusion_bits);

  int* getPatchReadyQueue();

  void reallocate_forceSOA(int atomStorageSize); 
};

#endif // NAMD_CUDA
#endif // CUDACOMPUTENONBONDEDKERNEL_H