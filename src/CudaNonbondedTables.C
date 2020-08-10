#include "charm++.h"
#include "NamdTypes.h"
#include "ComputeNonbondedUtil.h"
#include "LJTable.h"
#include "CudaUtils.h"
#include "CudaNonbondedTables.h"

#ifdef NAMD_CUDA

CudaNonbondedTables::CudaNonbondedTables(const int deviceID) : deviceID(deviceID) {

  vdwCoefTable = NULL;
  vdwCoefTableWidth = 0;
  vdwCoefTableTex = 0;  

  forceTableTex = 0;
  energyTableTex = 0;

  exclusionTable = NULL;
  r2_table = NULL;
  exclusionTableTex = 0;
  r2_table_tex = 0;

  modifiedExclusionForceTableTex = 0;
  modifiedExclusionEnergyTableTex = 0;

  cudaCheck(cudaSetDevice(deviceID));
  buildForceAndEnergyTables(4096);
  buildVdwCoefTable();
}

CudaNonbondedTables::~CudaNonbondedTables() {
  cudaCheck(cudaSetDevice(deviceID));
  if (vdwCoefTable != NULL) deallocate_device<float2>(&vdwCoefTable);
  if (exclusionTable != NULL) deallocate_device<float4>(&exclusionTable);
  if (r2_table != NULL) deallocate_device<float>(&r2_table);

  cudaCheck(cudaFreeArray(forceArray));
  cudaCheck(cudaFreeArray(energyArray));

  cudaCheck(cudaDestroyTextureObject(vdwCoefTableTex));
  cudaCheck(cudaDestroyTextureObject(forceTableTex));
  cudaCheck(cudaDestroyTextureObject(energyTableTex));

  cudaCheck(cudaDestroyTextureObject(exclusionTableTex));
  cudaCheck(cudaDestroyTextureObject(r2_table_tex));
}

template <typename T>
void bindTextureObject(int size, T* h_table, T*& d_table, cudaTextureObject_t& tex, bool update=false) {
  // Copy to device
  if ( ! update) {
    allocate_device<T>(&d_table, size);
  }
  else {
    cudaCheck(cudaDestroyTextureObject(tex));
  }
  copy_HtoD_sync<T>(h_table, d_table, size);

  // Create texture object
  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeLinear;
  resDesc.res.linear.devPtr = d_table;
  resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
  resDesc.res.linear.desc.x = sizeof(float)*8; // bits per channel
  if (sizeof(T) >= sizeof(float)*2) resDesc.res.linear.desc.y = sizeof(float)*8; // bits per channel
  if (sizeof(T) >= sizeof(float)*3) resDesc.res.linear.desc.z = sizeof(float)*8; // bits per channel
  if (sizeof(T) >= sizeof(float)*4) resDesc.res.linear.desc.w = sizeof(float)*8; // bits per channel
  resDesc.res.linear.sizeInBytes = size*sizeof(T);

  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.readMode = cudaReadModeElementType;

  cudaCheck(cudaCreateTextureObject(&tex, &resDesc, &texDesc, NULL));
}

//
// Builds the VdW Lennard-Jones coefficient table. Swiped from ComputeNonbondedCUDA.C
// NOTE: Can only be called once
//
void CudaNonbondedTables::buildVdwCoefTable(bool update/*=false*/) {

  const LJTable* const ljTable = ComputeNonbondedUtil:: ljTable;
  const int dim = ljTable->get_table_dim();

  // round dim up to odd multiple of 16
  int tsize = (((dim+16+31)/32)*32)-16;
  if ( tsize < dim ) NAMD_bug("CudaNonbondedTables::buildVdwCoefTable bad tsize");

  float2 *h_vdwCoefTable = new float2[tsize*tsize];
  float2 *h_exclusionVdwCoefTable = new float2[tsize*tsize];
  float2 *row = h_vdwCoefTable;
  float2 *exclusionRow = h_exclusionVdwCoefTable;
  for ( int i=0; i<dim; ++i, row += tsize, exclusionRow += tsize ) {
    for ( int j=0; j<dim; ++j ) {
      const LJTable::TableEntry *e = ljTable->table_val(i,j);
      row[j].x = e->A * ComputeNonbondedUtil::scaling;
      row[j].y = e->B * ComputeNonbondedUtil::scaling;
      exclusionRow[j].x = ((e+1)->A - e->A)* ComputeNonbondedUtil::scaling;
      exclusionRow[j].y = ((e+1)->B - e->B)* ComputeNonbondedUtil::scaling;
    }
  }

  vdwCoefTableWidth = tsize;

  bindTextureObject<float2>(tsize*tsize, h_vdwCoefTable, vdwCoefTable, vdwCoefTableTex, update);
  bindTextureObject<float2>(tsize*tsize, h_exclusionVdwCoefTable, exclusionVdwCoefTable, exclusionVdwCoefTableTex, update);

  delete [] h_vdwCoefTable;
  delete [] h_exclusionVdwCoefTable;

  if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Updated CUDA LJ table with %d x %d elements.\n", dim, dim);
  }
}

// Update tables from newer CPU values
void CudaNonbondedTables::updateTables() {
  buildVdwCoefTable(true);
}

//
// Builds force and energy tables. Swiped from ComputeNonbondedCUDA.C
//
template <typename T>
void buildForceAndEnergyTable(const int tableSize, const double* r2list, const BigReal* src_table, const bool flip,
  const BigReal prefac, const int dst_stride, T* dst_force, T* dst_energy) {

  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  union { double f; int32 i[2]; } byte_order_test;
  byte_order_test.f = 1.0;  // should occupy high-order bits only
  int32 *r2iilist = (int32*)r2list + ( byte_order_test.i[0] ? 0 : 1 );

  for ( int i=1; i<tableSize; ++i ) {
    double r = ((double) tableSize) / ( (double) i + 0.5 );
    int table_i = (r2iilist[2*i] >> 14) + r2_delta_expc;  // table_i >= 0

    if ( r > ComputeNonbondedUtil::cutoff ) {
      dst_force[i*dst_stride] = 0.;
      dst_energy[i*dst_stride] = 0.;
      continue;
    }

    BigReal diffa = r2list[i] - ComputeNonbondedUtil::r2_table[table_i];

    BigReal table_a, table_b, table_c, table_d;
    if (flip) {
      table_a = src_table[4*table_i+3];
      table_b = src_table[4*table_i+2];
      table_c = src_table[4*table_i+1];
      table_d = src_table[4*table_i];
    } else {
      table_a = src_table[4*table_i];
      table_b = src_table[4*table_i+1];
      table_c = src_table[4*table_i+2];
      table_d = src_table[4*table_i+3];
    }

    BigReal grad = ( 3. * table_d * diffa + 2. * table_c ) * diffa + table_b;
    dst_force[i*dst_stride] = prefac * 2. * grad;
    BigReal ener = table_a + diffa * ( ( table_d * diffa + table_c ) * diffa + table_b);
    dst_energy[i*dst_stride] = prefac * ener;
  }

  dst_force[0] = 0.;
  dst_energy[0] = dst_energy[1*dst_stride];
}

template <typename T>
void bindTextureObject(int tableSize, int tableWidth, T* h_table,
  cudaArray_t& array, cudaTextureObject_t& tableTex) {

  cudaChannelFormatDesc desc;
  memset(&desc, 0, sizeof(desc));
  desc.x = sizeof(T)*8;
  if (tableWidth >= 2) desc.y = sizeof(T)*8;
  if (tableWidth >= 3) desc.z = sizeof(T)*8;
  if (tableWidth >= 4) desc.w = sizeof(T)*8;
  desc.f = cudaChannelFormatKindFloat;
  cudaCheck(cudaMallocArray(&array, &desc, tableSize, 1));
  cudaCheck(cudaMemcpyToArray(array, 0, 0, h_table, tableSize*sizeof(T)*tableWidth, cudaMemcpyHostToDevice));

  cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = array;

  cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeClamp;
  texDesc.filterMode = cudaFilterModeLinear;
  texDesc.normalizedCoords = 1;

  cudaCheck(cudaCreateTextureObject(&tableTex, &resDesc, &texDesc, NULL));
}

void CudaNonbondedTables::buildForceAndEnergyTables(int tableSize) {

  // Build r2list
  const BigReal r2_delta = ComputeNonbondedUtil:: r2_delta;
  const int r2_delta_exp = ComputeNonbondedUtil:: r2_delta_exp;
  const int r2_delta_expc = 64 * (r2_delta_exp - 1023);

  double* r2list = new double[tableSize];  // double to match cpu code
  for ( int i=1; i<tableSize; ++i ) {
    double r = ((double) tableSize) / ( (double) i + 0.5 );
    r2list[i] = r*r + r2_delta;
  }

  // Non-bonded
  {
    float4* t = new float4[tableSize];
    float4* et = new float4[tableSize];

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::fast_table, false, 1.0,
      4, &t[0].x, &et[0].x);

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::vdwb_table, false, -1.0,
      4, &t[0].y, &et[0].y);

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::vdwa_table, false, 1.0,
      4, &t[0].z, &et[0].z);

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::scor_table, false, 1.0,
      4, &t[0].w, &et[0].w);

    bindTextureObject<float>(tableSize, 4, (float *)t, forceArray, forceTableTex);
    bindTextureObject<float>(tableSize, 4, (float *)et, energyArray, energyTableTex);
    delete [] t;
    delete [] et;
  }

  // Modified exclusions
  {
    float4* t = new float4[tableSize];
    float4* et = new float4[tableSize];

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::fast_table, false, -1.0,
      4, &t[0].x, &et[0].x);

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::vdwb_table, false, -1.0,
      4, &t[0].y, &et[0].y);

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::vdwa_table, false, 1.0,
      4, &t[0].z, &et[0].z);

    buildForceAndEnergyTable<float>(tableSize, r2list, ComputeNonbondedUtil::slow_table, true, -1.0,
      4, &t[0].w, &et[0].w);
    // delete [] flip_slow_table;

    bindTextureObject<float>(tableSize, 4, (float *)t, modifiedExclusionForceArray, modifiedExclusionForceTableTex);
    bindTextureObject<float>(tableSize, 4, (float *)et, modifiedExclusionEnergyArray, modifiedExclusionEnergyTableTex);
    delete [] t;
    delete [] et;
  }

  // Exclusions
  {
    BigReal* corr_full_table = new BigReal[ComputeNonbondedUtil::table_length*4];
    for (int i=0;i < ComputeNonbondedUtil::table_length*4;i++) {
      corr_full_table[i] =  ComputeNonbondedUtil::corr_table[i] - ComputeNonbondedUtil::full_table[i];
    }

    float4* h_exclusionTable = new float4[ComputeNonbondedUtil::table_length];
    float* h_r2_table = new float[ComputeNonbondedUtil::table_length];
    for (int i=0;i < ComputeNonbondedUtil::table_length;i++) {
      h_exclusionTable[i].x = 6.0*corr_full_table[4*i + 3];
      h_exclusionTable[i].y = 4.0*corr_full_table[4*i + 2];
      h_exclusionTable[i].z = 2.0*corr_full_table[4*i + 1];
      h_exclusionTable[i].w = 1.0*corr_full_table[4*i + 0];
      h_r2_table[i] = ComputeNonbondedUtil::r2_table[i];
    }

    bindTextureObject<float>(ComputeNonbondedUtil::table_length, h_r2_table, r2_table, r2_table_tex);
    bindTextureObject<float4>(ComputeNonbondedUtil::table_length, h_exclusionTable, exclusionTable, exclusionTableTex);

    delete [] h_exclusionTable;
    delete [] h_r2_table;

  }

  if ( ! CmiPhysicalNodeID(CkMyPe()) ) {
    CkPrintf("Info: Updated CUDA force table with %d elements.\n", tableSize);
  }

  delete [] r2list;
}

#endif // NAMD_CUDA
