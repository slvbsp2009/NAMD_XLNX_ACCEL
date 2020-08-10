#ifndef COMPUTEBONDEDCUDAKERNEL_H
#define COMPUTEBONDEDCUDAKERNEL_H
#include "CudaUtils.h"
#include "TupleTypesCUDA.h"
#include "CudaNonbondedTables.h"

#ifdef NAMD_CUDA

// Use Fixed point (24.40) force?
#define USE_FP_FORCE
#define FORCE_TYPE long long int
#define USE_STRIDED_FORCE

#ifndef USE_STRIDED_FORCE
#error "Non-USE_STRIDED_FORCE not implemented"
#endif

// Use Fixed point (34.30) virial?
#define USE_FP_VIRIAL
#ifdef USE_FP_VIRIAL
#define VIRIAL_TYPE long long int
#else
#define VIRIAL_TYPE double
#endif

#define WRITE_FULL_VIRIALS

// Scaling factors for 24.40 fixed point
#ifdef USE_FP_FORCE
static __constant__ const float float_to_force = (float)(1ll << 40);
static __constant__ const float force_to_float = (float)1.0/(float)(1ll << 40);
static __constant__ const double force_to_double = (double)1.0/(double)(1ll << 40);
#else
static __constant__ const float float_to_force = 1.0f;
static __constant__ const float force_to_float = 1.0f;
static __constant__ const double force_to_double = 1.0;
#endif

#ifdef USE_FP_VIRIAL
static __constant__ const float float_to_virial = (float)(1ll << 30);
static __constant__ const double double_to_virial = (double)(1ll << 30);
static __constant__ const double virial_to_double = (double)1.0/(double)(1ll << 30);
static __constant__ const long long int CONVERT_TO_VIR = (1ll << 10);
#endif

class ComputeBondedCUDAKernel {
public:

  // Enumeration for energies_virials[]
  enum {energyIndex_BOND=0, energyIndex_ANGLE, energyIndex_DIHEDRAL, energyIndex_IMPROPER,
    energyIndex_ELECT, energyIndex_LJ, energyIndex_ELECT_SLOW, energyIndex_CROSSTERM,
    normalVirialIndex_XX, normalVirialIndex_XY, normalVirialIndex_XZ,
    normalVirialIndex_YX, normalVirialIndex_YY, normalVirialIndex_YZ,
    normalVirialIndex_ZX, normalVirialIndex_ZY, normalVirialIndex_ZZ,
    nbondVirialIndex_XX, nbondVirialIndex_XY, nbondVirialIndex_XZ,
    nbondVirialIndex_YX, nbondVirialIndex_YY, nbondVirialIndex_YZ,
    nbondVirialIndex_ZX, nbondVirialIndex_ZY, nbondVirialIndex_ZZ,
    slowVirialIndex_XX, slowVirialIndex_XY, slowVirialIndex_XZ,
    slowVirialIndex_YX, slowVirialIndex_YY, slowVirialIndex_YZ,
    slowVirialIndex_ZX, slowVirialIndex_ZY, slowVirialIndex_ZZ,
    amdDiheVirialIndex_XX, amdDiheVirialIndex_XY, amdDiheVirialIndex_XZ,
    amdDiheVirialIndex_YX, amdDiheVirialIndex_YY, amdDiheVirialIndex_YZ,
    amdDiheVirialIndex_ZX, amdDiheVirialIndex_ZY, amdDiheVirialIndex_ZZ,
    energies_virials_SIZE};

  template <typename T>
  struct BondedVirial {
#ifdef WRITE_FULL_VIRIALS
    T xx;
    T xy;
    T xz;
    T yx;
    T yy;
    T yz;
    T zx;
    T zy;
    T zz;
#else
#error "non-WRITE_FULL_VIRIALS not implemented yet"
    union {
      double sforce_dp[27][3];
      long long int sforce_fp[27][3];
    };
#endif
  };

private:
  const int deviceID;
  CudaNonbondedTables& cudaNonbondedTables;

  // This stores all bonds, angles, dihedrals, and impropers in a single 
  // contigious memory array.
  char* tupleData;
  int tupleDataSize;

  // ---------------------------------------------------------------------------------
  // NOTE: bonds, angles, dihedrals, impropers, etc. - pointers below are 
  // computed pointers pointing to tupleData -array
  // DO NOT DEALLOCATE THESE!
  int numBonds;
  CudaBond* bonds;

  int numAngles;
  CudaAngle* angles;

  int numDihedrals;
  CudaDihedral* dihedrals;

  int numImpropers;
  CudaDihedral* impropers;

  int numModifiedExclusions;
  CudaExclusion* modifiedExclusions;

  int numExclusions;
  CudaExclusion* exclusions;

  int numCrossterms;
  CudaCrossterm* crossterms;
  // ---------------------------------------------------------------------------------
  
  // Device memory for coordinates
  float4* xyzq;
  int xyzqSize;

  // Device memory for forces:
  // [normal, nbond, slow]
  FORCE_TYPE* forces;
  int forcesSize;

  CudaBondValue* bondValues;
  CudaAngleValue* angleValues;
  CudaDihedralValue* dihedralValues;
  CudaDihedralValue* improperValues;
  CudaCrosstermValue* crosstermValues;

  // Accumulated energy values for every bonded type
  double* energies_virials;

public:

  ComputeBondedCUDAKernel(int deviceID, CudaNonbondedTables& cudaNonbondedTables);
  ~ComputeBondedCUDAKernel();

  static int warpAlign(const int n) {return ((n + WARPSIZE - 1)/WARPSIZE)*WARPSIZE;} 

  void update(
    const int numBondsIn,
    const int numAnglesIn,
    const int numDihedralsIn,
    const int numImpropersIn,
    const int numModifiedExclusionsIn,
    const int numExclusionsIn,
    const int numCrosstermsIn,
    const char* h_tupleData,
    cudaStream_t stream);

  void setupBondValues(int numBondValues, CudaBondValue* h_bondValues);
  void setupAngleValues(int numAngleValues, CudaAngleValue* h_angleValues);
  void setupDihedralValues(int numDihedralValues, CudaDihedralValue* h_dihedralValues);
  void setupImproperValues(int numImproperValues, CudaDihedralValue* h_improperValues);
  void setupCrosstermValues(int numCrosstermValues, CudaCrosstermValue* h_crosstermValues);

  int getForceStride(const int atomStorageSize);
  int getForceSize(const int atomStorageSize);
  int getAllForceSize(const int atomStorageSize, const bool doSlow);

  void bondedForce(
    const double scale14, const int atomStorageSize,
    const bool doEnergy, const bool doVirial, const bool doSlow,
    const float3 lata, const float3 latb, const float3 latc,
    const float cutoff2, const float r2_delta, const int r2_delta_expc,
    const float4* h_xyzq, FORCE_TYPE* h_forces,
    double *h_energies,
    cudaStream_t stream);

};

#endif

#endif // COMPUTEBONDEDCUDAKERNEL_H