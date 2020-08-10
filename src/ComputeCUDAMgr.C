#include "NamdTypes.h"
#include "common.h"
#include "Node.h"
#include "ComputeCUDAMgr.h"

#include "DeviceCUDA.h"
#ifdef NAMD_CUDA
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;

//
// Class constructor
//
ComputeCUDAMgr::ComputeCUDAMgr() {
	// __sdag_init();
  numDevices = 0;
  // numNodesContributed = 0;
  // numDevicesMax = 0;
}

//
// Class constructor
//
ComputeCUDAMgr::ComputeCUDAMgr(CkMigrateMessage *) {
	// __sdag_init();
  NAMD_bug("ComputeCUDAMgr cannot be migrated");
  numDevices = 0;
  // numNodesContributed = 0;
  // numDevicesMax = 0;
}

//
// Class destructor
//
ComputeCUDAMgr::~ComputeCUDAMgr() {
  for (int i=0;i < numDevices;i++) {
    if (cudaNonbondedTablesList[i] != NULL) delete cudaNonbondedTablesList[i];
    if (cudaComputeNonbondedList[i] != NULL) delete cudaComputeNonbondedList[i];
#ifdef BONDED_CUDA
    if (computeBondedCUDAList[i] != NULL) delete computeBondedCUDAList[i];
#endif
  }
}

//
// Initialize manager
// This gets called on rank 0 of each node
//
void ComputeCUDAMgr::initialize(CkQdMsg *msg) {
	if (msg != NULL) delete msg;

	numDevices = deviceCUDA->getDeviceCount();

  // Create pointers to devices
  cudaNonbondedTablesList.resize(numDevices, NULL);
  cudaComputeNonbondedList.resize(numDevices, NULL);
#ifdef BONDED_CUDA
  computeBondedCUDAList.resize(numDevices, NULL);
#endif

  // Create CUDA non-bonded tables for all devices that are used for computation
  for (int i=0;i < deviceCUDA->getNumDevice();i++) {
    int deviceID = deviceCUDA->getDeviceIDbyRank(i);
    cudaNonbondedTablesList[deviceID] = new CudaNonbondedTables(deviceID);
  }
}

//
// Update nonbonded tables
// Should be called only on rank 0 of each node
//
void ComputeCUDAMgr::update() {
  if ( CkMyRank() ) NAMD_bug("ComputeCUDAMgr::update() should be called only by rank 0");
  for (int i=0;  i < deviceCUDA->getNumDevice();  i++) {
    int deviceID = deviceCUDA->getDeviceIDbyRank(i);
    // calls update function from CudaNonbondedTables
    cudaNonbondedTablesList[deviceID]->updateTables();
  }
}

ComputeCUDAMgr* ComputeCUDAMgr::getComputeCUDAMgr() {
  // Get pointer to ComputeCUDAMgr on this node
  CProxy_ComputeCUDAMgr computeCUDAMgrProxy = CkpvAccess(BOCclass_group).computeCUDAMgr;
  ComputeCUDAMgr* computeCUDAMgr = computeCUDAMgrProxy.ckLocalBranch();
  if (computeCUDAMgr == NULL)
    NAMD_bug("getComputeCUDAMgr, unable to locate local branch of BOC entry ComputeCUDAMgr");
  return computeCUDAMgr;
}

//
// Creates CudaComputeNonbonded object
//
CudaComputeNonbonded* ComputeCUDAMgr::createCudaComputeNonbonded(ComputeID c) {
  int deviceID = deviceCUDA->getDeviceID();
  if (cudaComputeNonbondedList.at(deviceID) != NULL)
    NAMD_bug("ComputeCUDAMgr::createCudaComputeNonbonded called twice");
  if (cudaNonbondedTablesList.at(deviceID) == NULL)
    NAMD_bug("ComputeCUDAMgr::createCudaComputeNonbonded, non-bonded CUDA tables not created");
  bool doStreaming = !deviceCUDA->getNoStreaming() && !Node::Object()->simParameters->GBISOn;
  cudaComputeNonbondedList[deviceID] = new CudaComputeNonbonded(c, deviceID, *cudaNonbondedTablesList[deviceID], doStreaming);
  return cudaComputeNonbondedList[deviceID];
}

//
// Returns CudaComputeNonbonded for this Pe
//
CudaComputeNonbonded* ComputeCUDAMgr::getCudaComputeNonbonded() {
  // Get device ID for this Pe
  int deviceID = deviceCUDA->getDeviceID();
  CudaComputeNonbonded* p = cudaComputeNonbondedList[deviceID];
  if (p == NULL)
    NAMD_bug("ComputeCUDAMgr::getCudaComputeNonbonded(), device not created yet");
  return p;
}

#ifdef BONDED_CUDA
//
// Creates ComputeBondedCUDA object
//
ComputeBondedCUDA* ComputeCUDAMgr::createComputeBondedCUDA(ComputeID c, ComputeMgr* computeMgr) {
  int deviceID = deviceCUDA->getDeviceID();
  if (computeBondedCUDAList.at(deviceID) != NULL)
    NAMD_bug("ComputeCUDAMgr::createComputeBondedCUDA called twice");
  if (cudaNonbondedTablesList.at(deviceID) == NULL)
    NAMD_bug("ComputeCUDAMgr::createCudaComputeNonbonded, non-bonded CUDA tables not created");
  computeBondedCUDAList[deviceID] = new ComputeBondedCUDA(c, computeMgr, deviceID, *cudaNonbondedTablesList[deviceID]);
  return computeBondedCUDAList[deviceID];
}

//
// Returns ComputeBondedCUDA for this Pe
//
ComputeBondedCUDA* ComputeCUDAMgr::getComputeBondedCUDA() {
  // Get device ID for this Pe
  int deviceID = deviceCUDA->getDeviceID();
  ComputeBondedCUDA* p = computeBondedCUDAList[deviceID];
  if (p == NULL)
    NAMD_bug("ComputeCUDAMgr::getComputeBondedCUDA(), device not created yet");
  return p;
}
#endif // BONDED_CUDA

#endif // NAMD_CUDA

#include "ComputeCUDAMgr.def.h"
