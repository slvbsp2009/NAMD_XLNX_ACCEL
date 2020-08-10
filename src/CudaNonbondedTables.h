#ifndef CUDANONBONDEDTABLES_H
#define CUDANONBONDEDTABLES_H

#ifdef NAMD_CUDA
#include <cuda.h>

class CudaNonbondedTables {
private:
  const int deviceID;

  float2 *vdwCoefTable;
  int vdwCoefTableWidth;
  cudaTextureObject_t vdwCoefTableTex;

  // Non-bonded
  cudaArray_t forceArray;
  cudaTextureObject_t forceTableTex;

  cudaArray_t energyArray;
  cudaTextureObject_t energyTableTex;

  // Modified exclusions
  cudaArray_t modifiedExclusionForceArray;
  cudaTextureObject_t modifiedExclusionForceTableTex;
  
  cudaArray_t modifiedExclusionEnergyArray;
  cudaTextureObject_t modifiedExclusionEnergyTableTex;

  // Exclusions
  float2 *exclusionVdwCoefTable;
  cudaTextureObject_t exclusionVdwCoefTableTex;

  // cudaArray_t exclusionForceArray;
  // cudaTextureObject_t exclusionForceTableTex;
  
  // cudaArray_t exclusionEnergyArray;
  // cudaTextureObject_t exclusionEnergyTableTex;

  float4* exclusionTable;
  float* r2_table;
  cudaTextureObject_t exclusionTableTex;
  cudaTextureObject_t r2_table_tex;

  void buildVdwCoefTable(bool update=false);
  void buildForceAndEnergyTables(int tableSize);

public:
  CudaNonbondedTables(const int deviceID);
  ~CudaNonbondedTables();

  float2* getVdwCoefTable() {return vdwCoefTable;}
  int getVdwCoefTableWidth() {return vdwCoefTableWidth;}
  cudaTextureObject_t getVdwCoefTableTex() {return vdwCoefTableTex;}
  cudaTextureObject_t getForceTableTex() {return forceTableTex;}
  cudaTextureObject_t getEnergyTableTex() {return energyTableTex;}

  void updateTables();

  float2* getExclusionVdwCoefTable() {return exclusionVdwCoefTable;}
  cudaTextureObject_t getExclusionVdwCoefTableTex() {return exclusionVdwCoefTableTex;}
  // cudaTextureObject_t getExclusionForceTableTex() {return exclusionForceTableTex;}
  // cudaTextureObject_t getExclusionEnergyTableTex() {return exclusionEnergyTableTex;}
  float4* getExclusionTable() {return exclusionTable;}
  float* get_r2_table() {return r2_table;}
  cudaTextureObject_t getExclusionTableTex() {return exclusionTableTex;}
  cudaTextureObject_t get_r2_table_tex() {return r2_table_tex;}

  cudaTextureObject_t getModifiedExclusionForceTableTex() {return modifiedExclusionForceTableTex;}
  cudaTextureObject_t getModifiedExclusionEnergyTableTex() {return modifiedExclusionEnergyTableTex;}

};

#endif // NAMD_CUDA
#endif // CUDANONBONDEDTABLES_H
