#include "charm++.h"
#include "NamdTypes.h"
#include "ComputeMgr.h"
#include "WorkDistrib.h"
#include "ProxyMgr.h"
#include "CudaUtils.h"
#include "DeviceCUDA.h"
#include "ComputeBonds.h"
#include "ComputeAngles.h"
#include "ComputeDihedrals.h"
#include "ComputeImpropers.h"
#include "ComputeCrossterms.h"
#include "ComputeNonbondedCUDAExcl.h"
#include "ComputeBondedCUDA.h"

#ifdef NAMD_CUDA
#ifdef BONDED_CUDA

#include <algorithm>   // std::find

const int ComputeBondedCUDA::CudaTupleTypeSize[Tuples::NUM_TUPLE_TYPES] = {
  sizeof(CudaBond),          // Bonds
  sizeof(CudaAngle),         // Angles
  sizeof(CudaDihedral),      // Dihedrals
  sizeof(CudaDihedral),      // Impropers
  sizeof(CudaExclusion),     // Exclusions
  sizeof(CudaCrossterm)      // Crossterms
};

extern "C" void CcdCallBacksReset(void *ignored, double curWallTime);  // fix Charm++

//
// Class constructor
//
ComputeBondedCUDA::ComputeBondedCUDA(ComputeID c, ComputeMgr* computeMgr, int deviceID,
  CudaNonbondedTables& cudaNonbondedTables) :
Compute(c), computeMgr(computeMgr), deviceID(deviceID), masterPe(CkMyPe()),
bondedKernel(deviceID, cudaNonbondedTables) 
{

  computes.resize(CkMyNodeSize());
  patchIDsPerRank.resize(CkMyNodeSize());
  numExclPerRank.resize(CkMyNodeSize());
  for (int i=0;i < numExclPerRank.size();i++) {
    numExclPerRank[i].numModifiedExclusions = 0;
    numExclPerRank[i].numExclusions = 0;
  }

  atomMap.allocateMap(Node::Object()->molecule->numAtoms);

  flags = NULL;

  tupleData = NULL;
  tupleDataSize = 0;

  atoms = NULL;
  atomsSize = 0;

  forces = NULL;
  forcesSize = 0;

  energies_virials = NULL;

  initializeCalled = false;

  SimParameters *params = Node::Object()->simParameters;
  accelMDdoDihe = false;
  if (params->accelMDOn) {
     if (params->accelMDdihe || params->accelMDdual) accelMDdoDihe=true;
  }
}

//
// Class destructor
//
ComputeBondedCUDA::~ComputeBondedCUDA() {
  cudaCheck(cudaSetDevice(deviceID));

  if (atoms != NULL) deallocate_host<CudaAtom>(&atoms);
  if (forces != NULL) deallocate_host<FORCE_TYPE>(&forces);
  if (energies_virials != NULL) deallocate_host<double>(&energies_virials);
  if (tupleData != NULL) deallocate_host<char>(&tupleData);

  if (initializeCalled) {
    cudaCheck(cudaStreamDestroy(stream));
    cudaCheck(cudaEventDestroy(forceDoneEvent));
    CmiDestroyLock(lock);
    delete reduction;
  }

  // NOTE: unregistering happens in [sync] -entry method
  computeMgr->sendUnregisterBoxesOnPe(pes, this);
}

void ComputeBondedCUDA::unregisterBoxesOnPe() {
  for (int i=0;i < patchIDsPerRank[CkMyRank()].size();i++) {
    PatchID patchID = patchIDsPerRank[CkMyRank()][i];
    TuplePatchElem* tpe = tuplePatchList.find(TuplePatchElem(patchID));
    if (tpe == NULL || tpe->p == NULL) {
      NAMD_bug("ComputeBondedCUDA::unregisterBoxesOnPe, TuplePatchElem not found or setup incorrectly");
    }
    Patch* patch = tpe->p;
    if (tpe->positionBox != NULL) patch->unregisterPositionPickup(this, &tpe->positionBox);
    if (tpe->avgPositionBox != NULL) patch->unregisterAvgPositionPickup(this, &tpe->avgPositionBox);
    if (tpe->forceBox != NULL) patch->unregisterForceDeposit(this, &tpe->forceBox);
  }
}

//
// Register compute for a given PE. pids is a list of patches the PE has
// Called by master PE
//
void ComputeBondedCUDA::registerCompute(int pe, int type, PatchIDList& pids) {

  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::registerCompute() called on non master PE");

  int rank = CkRankOf(pe);

  HomeCompute& homeCompute = computes[rank].homeCompute;
  if (homeCompute.patchIDs.size() == 0) {
    homeCompute.isBasePatch.resize(PatchMap::Object()->numPatches(), 0);
    homeCompute.patchIDs.resize(pids.size());
    for (int i=0;i < pids.size();i++) {
      homeCompute.patchIDs[i] = pids[i];
      homeCompute.isBasePatch[pids[i]] = 1;
    }
  } else {
    if (homeCompute.patchIDs.size() != pids.size()) {
      NAMD_bug("ComputeBondedCUDA::registerCompute(), homeComputes, patch IDs do not match (1)");
    }
    for (int i=0;i < pids.size();i++) {
      if (homeCompute.patchIDs[i] != pids[i]) {
        NAMD_bug("ComputeBondedCUDA::registerCompute(), homeComputes, patch IDs do not match (2)");
      }
    }
  }

  switch(type) {
    case computeBondsType:
    homeCompute.tuples.push_back(new HomeTuples<BondElem, Bond, BondValue>(Tuples::BOND));
    break;

    case computeAnglesType:
    homeCompute.tuples.push_back(new HomeTuples<AngleElem, Angle, AngleValue>(Tuples::ANGLE));
    break;

    case computeDihedralsType:
    homeCompute.tuples.push_back(new HomeTuples<DihedralElem, Dihedral, DihedralValue>(Tuples::DIHEDRAL));
    break;

    case computeImpropersType:
    homeCompute.tuples.push_back(new HomeTuples<ImproperElem, Improper, ImproperValue>(Tuples::IMPROPER));
    break;

    case computeExclsType:
    homeCompute.tuples.push_back(new HomeTuples<ExclElem, Exclusion, int>(Tuples::EXCLUSION));
    break;

    case computeCrosstermsType:
    homeCompute.tuples.push_back(new HomeTuples<CrosstermElem, Crossterm, CrosstermValue>(Tuples::CROSSTERM));
    break;

    default:
    NAMD_bug("ComputeBondedCUDA::registerCompute(), Unsupported compute type");
    break;
  }

}

//
// Register self compute for a given PE
// Called by master PE
//
void ComputeBondedCUDA::registerSelfCompute(int pe, int type, int pid) {

  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::registerSelfCompute() called on non master PE");

  int rank = CkRankOf(pe);

  std::vector< SelfCompute >& selfComputes = computes[rank].selfComputes;
  auto it = find(selfComputes.begin(), selfComputes.end(), SelfCompute(type));
  if (it == selfComputes.end()) {
    // Type not found, add new one
    selfComputes.push_back(SelfCompute(type));
    it = selfComputes.begin() + (selfComputes.size() - 1);

    switch(type) {
      case computeSelfBondsType:
      it->tuples = new SelfTuples<BondElem, Bond, BondValue>(Tuples::BOND);
      break;

      case computeSelfAnglesType:
      it->tuples = new SelfTuples<AngleElem, Angle, AngleValue>(Tuples::ANGLE);
      break;

      case computeSelfDihedralsType:
      it->tuples = new SelfTuples<DihedralElem, Dihedral, DihedralValue>(Tuples::DIHEDRAL);
      break;

      case computeSelfImpropersType:
      it->tuples = new SelfTuples<ImproperElem, Improper, ImproperValue>(Tuples::IMPROPER);
      break;

      case computeSelfExclsType:
      it->tuples = new SelfTuples<ExclElem, Exclusion, int>(Tuples::EXCLUSION);
      break;

      case computeSelfCrosstermsType:
      it->tuples = new SelfTuples<CrosstermElem, Crossterm, CrosstermValue>(Tuples::CROSSTERM);
      break;

      default:
      NAMD_bug("ComputeBondedCUDA::registerSelfCompute(), Unsupported compute type");
      break;
    }

  }

  // Add patch ID for this type
  it->patchIDs.push_back(pid);
}

void ComputeBondedCUDA::assignPatchesOnPe() {

  PatchMap* patchMap = PatchMap::Object();
  for (int i=0;i < patchIDsPerRank[CkMyRank()].size();i++) {
    PatchID patchID = patchIDsPerRank[CkMyRank()][i];
    ProxyMgr::Object()->createProxy(patchID);
    Patch* patch = patchMap->patch(patchID);
    if (patch == NULL)
      NAMD_bug("ComputeBondedCUDA::assignPatchesOnPe, patch not found");
    if (flags == NULL) flags = &patchMap->patch(patchID)->flags;
    TuplePatchElem* tpe = tuplePatchList.find(TuplePatchElem(patchID));
    if (tpe == NULL) {
      NAMD_bug("ComputeBondedCUDA::assignPatchesOnPe, TuplePatchElem not found");
    }
    if (tpe->p != NULL) {
      NAMD_bug("ComputeBondedCUDA::assignPatchesOnPe, TuplePatchElem already registered");
    }
    // Assign patch and register coordinates and forces manually
    tpe->p = patch;
    tpe->positionBox = patch->registerPositionPickup(this);
    tpe->avgPositionBox = patch->registerAvgPositionPickup(this);
    tpe->forceBox = patch->registerForceDeposit(this);
  }
}

//
// atomUpdate() can be called by any Pe
//
void ComputeBondedCUDA::atomUpdate() {
  atomsChangedIn = true;
}

//
// Enqueu doWork on masterPe and return "no work"
// Can be called by any Pe
//
int ComputeBondedCUDA::noWork() {
  computeMgr->sendMessageEnqueueWork(masterPe, this);
  return 1;
}

void ComputeBondedCUDA::messageEnqueueWork() {
  if (masterPe != CkMyPe())
    NAMD_bug("ComputeBondedCUDA::messageEnqueueWork() must be called from master PE");
  WorkDistrib::messageEnqueueWork(this);
}

//
// Sends open-box commands to PEs
// Called on master PE
//
void ComputeBondedCUDA::doWork() {

  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::doWork() called on non master PE");

  // Read value of atomsChangedIn, which is set in atomUpdate(), and reset it.
  // atomsChangedIn can be set to true by any Pe
  // atomsChanged can only be set by masterPe
  // This use of double varibles makes sure we don't have race condition
  atomsChanged = atomsChangedIn;
  atomsChangedIn = false;

  if (getNumPatches() == 0) return;  // No work do to

  if (flags == NULL)
    NAMD_bug("ComputeBondedCUDA::doWork(), no flags set");

  // Read flags
  lattice = flags->lattice;
  doEnergy = flags->doEnergy;
  doVirial = flags->doVirial;
  doSlow   = flags->doFullElectrostatics;
  doMolly  = flags->doMolly;

  if (atomsChanged) {
    // Re-calculate patch atom numbers and storage
    updatePatches();
  }

  // Open boxes on Pes and launch work to masterPe
  computeMgr->sendOpenBoxesOnPe(pes, this);
}

//
// This gets called when patch finishes on a PE
//
void ComputeBondedCUDA::patchReady(PatchID pid, int doneMigration, int seq) {
  if (doneMigration) {
    // auto it = patchIndex.find(pid);
    // if (it == patchIndex.end())
    //   NAMD_bug("ComputeBondedCUDA::patchReady, Patch ID not found");
    patches[patchIndex[pid]].numAtoms = PatchMap::Object()->patch(pid)->getNumAtoms();
  }
  CmiLock(lock);
  Compute::patchReady(pid, doneMigration, seq);
  CmiUnlock(lock);
}

//
//
//
void ComputeBondedCUDA::updatePatches() {
  int atomStart = 0;
  for (int i=0;i < patches.size();i++) {
    patches[i].atomStart = atomStart;
    atomStart += patches[i].numAtoms;
  }
  atomStorageSize = atomStart;

  // Re-allocate atoms
  reallocate_host<CudaAtom>(&atoms, &atomsSize, atomStorageSize, 1.4f);
}

//
// Map atoms GPU-wide
//
void ComputeBondedCUDA::mapAtoms() {

  for (int i=0;i < getNumPatches();i++) {
    TuplePatchElem* tpe = tuplePatchList.find(TuplePatchElem(allPatchIDs[i]));
    atomMappers[i]->registerIDsCompAtomExt(tpe->xExt, tpe->xExt + tpe->p->getNumAtoms());
  }

}

//
// Unmap atoms GPU-wide
//
void ComputeBondedCUDA::unmapAtoms() {

  for (int i=0;i < getNumPatches();i++) {
    TuplePatchElem* tpe = tuplePatchList.find(TuplePatchElem(allPatchIDs[i]));
    atomMappers[i]->unregisterIDsCompAtomExt(tpe->xExt, tpe->xExt + tpe->p->getNumAtoms());
  }

}

//
// Open all patches that have been assigned to this Pe
//
void ComputeBondedCUDA::openBoxesOnPe() {

  std::vector<int>& patchIDs = patchIDsPerRank[CkMyRank()];

  for (auto it=patchIDs.begin();it != patchIDs.end();it++) {
    PatchID patchID = *it;
    TuplePatchElem* tpe = tuplePatchList.find(TuplePatchElem(patchID));
    tpe->x = tpe->positionBox->open();
    tpe->xExt = tpe->p->getCompAtomExtInfo();
    if ( doMolly ) tpe->x_avg = tpe->avgPositionBox->open();
    tpe->r = tpe->forceBox->open();
    tpe->f = tpe->r->f[Results::normal];
    if (accelMDdoDihe) tpe->af = tpe->r->f[Results::amdf]; // for dihedral-only or dual-boost accelMD

    // Copy atoms
    int pi = patchIndex[patchID];
    int atomStart = patches[pi].atomStart;
    int numAtoms = patches[pi].numAtoms;
    CompAtom* compAtom = tpe->x;
    const CompAtomExt *aExt = tpe->p->getCompAtomExtInfo();
    const CudaAtom *src = tpe->p->getCudaAtomList();
    for (int i=0;i < numAtoms;i++) {
      int j = aExt[i].sortOrder;
      atoms[atomStart + j] = src[i];
    }

  }

  bool done = false;
  CmiLock(lock);
  patchesCounter -= patchIDs.size();
  if (patchesCounter == 0) {
    patchesCounter = getNumPatches();
    done = true;
  }
  CmiUnlock(lock);
  if (done) {
    if (atomsChanged) {
      mapAtoms();
      computeMgr->sendLoadTuplesOnPe(pes, this);
    } else {
      computeMgr->sendLaunchWork(masterPe, this);
    }
  }
}

void countNumExclusions(Tuples* tuples, int& numModifiedExclusions, int& numExclusions) {
  numModifiedExclusions = 0;
  int ntuples = tuples->getNumTuples();
  ExclElem* src = (ExclElem *)(tuples->getTupleList());
  for (int ituple=0;ituple < ntuples;ituple++) {
    if (src[ituple].modified) numModifiedExclusions++;
  }
  numExclusions = ntuples - numModifiedExclusions;
}

//
// Load tuples on PE. Note: this can only after boxes on all PEs have been opened
//
void ComputeBondedCUDA::loadTuplesOnPe() {

  int numModifiedExclusions = 0;
  int numExclusions = 0;

  std::vector< SelfCompute >& selfComputes = computes[CkMyRank()].selfComputes;
  // Loop over self compute types
  for (auto it=selfComputes.begin();it != selfComputes.end();it++) {
    it->tuples->loadTuples(tuplePatchList, NULL, &atomMap, it->patchIDs);
    // For exclusions, we must count the number of modified and non-modified exclusions
    if (it->tuples->getType() == Tuples::EXCLUSION) {
      int tmp1, tmp2;
      countNumExclusions(it->tuples, tmp1, tmp2);
      numModifiedExclusions += tmp1;
      numExclusions += tmp2;
    }
  }

  HomeCompute& homeCompute = computes[CkMyRank()].homeCompute;
  for (int i=0;i < homeCompute.tuples.size();i++) {
    homeCompute.tuples[i]->loadTuples(tuplePatchList,
      homeCompute.isBasePatch.data(), &atomMap,
      homeCompute.patchIDs);
    // For exclusions, we must count the number of modified and non-modified exclusions
    if (homeCompute.tuples[i]->getType() == Tuples::EXCLUSION) {
      int tmp1, tmp2;
      countNumExclusions(homeCompute.tuples[i], tmp1, tmp2);
      numModifiedExclusions += tmp1;
      numExclusions += tmp2;
    }
  }

  // Store number of exclusions
  numExclPerRank[CkMyRank()].numModifiedExclusions = numModifiedExclusions;
  numExclPerRank[CkMyRank()].numExclusions         = numExclusions;

  bool done = false;
  CmiLock(lock);
  patchesCounter -= patchIDsPerRank[CkMyRank()].size();
  if (patchesCounter == 0) {
    patchesCounter = getNumPatches();
    done = true;
  }
  CmiUnlock(lock);
  if (done) {
    computeMgr->sendLaunchWork(masterPe, this);
  }
}

void ComputeBondedCUDA::copyBondData(const int ntuples, const BondElem* __restrict__ src,
  const BondValue* __restrict__ bond_array, CudaBond* __restrict__ dst) {

  PatchMap* patchMap = PatchMap::Object();
  for (int ituple=0;ituple < ntuples;ituple++) {
    CudaBond dstval;
    auto p0 = src[ituple].p[0];
    auto p1 = src[ituple].p[1];
    int pi0 = patchIndex[p0->patchID];
    int pi1 = patchIndex[p1->patchID];
    int l0 = src[ituple].localIndex[0];
    int l1 = src[ituple].localIndex[1];
    dstval.i = l0 + patches[pi0].atomStart;
    dstval.j = l1 + patches[pi1].atomStart;
    dstval.itype = (src[ituple].value - bond_array);
    Position position1 = p0->x[l0].position;
    Position position2 = p1->x[l1].position;
    Vector shiftVec = lattice.wrap_delta_scaled(position1, position2);
    shiftVec += patchMap->center(p0->patchID) - patchMap->center(p1->patchID);
    dstval.ioffsetXYZ = make_float3((float)shiftVec.x, (float)shiftVec.y, (float)shiftVec.z);
    dstval.scale = src[ituple].scale;
    dst[ituple] = dstval;
  }
}

void ComputeBondedCUDA::copyAngleData(const int ntuples, const AngleElem* __restrict__ src,
  const AngleValue* __restrict__ angle_array, CudaAngle* __restrict__ dst) {

  PatchMap* patchMap = PatchMap::Object();
  for (int ituple=0;ituple < ntuples;ituple++) {
    CudaAngle dstval;
    auto p0 = src[ituple].p[0];
    auto p1 = src[ituple].p[1];
    auto p2 = src[ituple].p[2];
    int pi0 = patchIndex[p0->patchID];
    int pi1 = patchIndex[p1->patchID];
    int pi2 = patchIndex[p2->patchID];
    int l0 = src[ituple].localIndex[0];
    int l1 = src[ituple].localIndex[1];
    int l2 = src[ituple].localIndex[2];
    dstval.i = l0 + patches[pi0].atomStart;
    dstval.j = l1 + patches[pi1].atomStart;
    dstval.k = l2 + patches[pi2].atomStart;
    dstval.itype = (src[ituple].value - angle_array);
    Position position1 = p0->x[l0].position;
    Position position2 = p1->x[l1].position;
    Position position3 = p2->x[l2].position;
    Vector shiftVec12 = lattice.wrap_delta_scaled(position1, position2);
    Vector shiftVec32 = lattice.wrap_delta_scaled(position3, position2);
    shiftVec12 += patchMap->center(p0->patchID) - patchMap->center(p1->patchID);
    shiftVec32 += patchMap->center(p2->patchID) - patchMap->center(p1->patchID);
    dstval.ioffsetXYZ = make_float3((float)shiftVec12.x, (float)shiftVec12.y, (float)shiftVec12.z);
    dstval.koffsetXYZ = make_float3((float)shiftVec32.x, (float)shiftVec32.y, (float)shiftVec32.z);
    dstval.scale = src[ituple].scale;
    dst[ituple] = dstval;
  }
}

//
// Used for both dihedrals and impropers
//
template <bool doDihedral, typename T, typename P>
void ComputeBondedCUDA::copyDihedralData(const int ntuples, const T* __restrict__ src,
  const P* __restrict__ p_array, CudaDihedral* __restrict__ dst) {

  PatchMap* patchMap = PatchMap::Object();
  for (int ituple=0;ituple < ntuples;ituple++) {
    CudaDihedral dstval;
    auto p0 = src[ituple].p[0];
    auto p1 = src[ituple].p[1];
    auto p2 = src[ituple].p[2];
    auto p3 = src[ituple].p[3];
    int pi0 = patchIndex[p0->patchID];
    int pi1 = patchIndex[p1->patchID];
    int pi2 = patchIndex[p2->patchID];
    int pi3 = patchIndex[p3->patchID];
    int l0 = src[ituple].localIndex[0];
    int l1 = src[ituple].localIndex[1];
    int l2 = src[ituple].localIndex[2];
    int l3 = src[ituple].localIndex[3];
    dstval.i = l0 + patches[pi0].atomStart;
    dstval.j = l1 + patches[pi1].atomStart;
    dstval.k = l2 + patches[pi2].atomStart;
    dstval.l = l3 + patches[pi3].atomStart;
    if (doDihedral) {
      dstval.itype = dihedralMultMap[(src[ituple].value - p_array)];
    } else {
      dstval.itype = improperMultMap[(src[ituple].value - p_array)];
    }
    Position position1 = p0->x[l0].position;
    Position position2 = p1->x[l1].position;
    Position position3 = p2->x[l2].position;
    Position position4 = p3->x[l3].position;
    Vector shiftVec12 = lattice.wrap_delta_scaled(position1, position2);
    Vector shiftVec23 = lattice.wrap_delta_scaled(position2, position3);
    Vector shiftVec43 = lattice.wrap_delta_scaled(position4, position3);
    shiftVec12 += patchMap->center(p0->patchID) - patchMap->center(p1->patchID);
    shiftVec23 += patchMap->center(p1->patchID) - patchMap->center(p2->patchID);
    shiftVec43 += patchMap->center(p3->patchID) - patchMap->center(p2->patchID);
    dstval.ioffsetXYZ = make_float3((float)shiftVec12.x, (float)shiftVec12.y, (float)shiftVec12.z);
    dstval.joffsetXYZ = make_float3((float)shiftVec23.x, (float)shiftVec23.y, (float)shiftVec23.z);
    dstval.loffsetXYZ = make_float3((float)shiftVec43.x, (float)shiftVec43.y, (float)shiftVec43.z);
    dstval.scale = src[ituple].scale;
    dst[ituple] = dstval;
  }
}

void ComputeBondedCUDA::copyExclusionData(const int ntuples, const ExclElem* __restrict__ src, const int typeSize,
  CudaExclusion* __restrict__ dst1, CudaExclusion* __restrict__ dst2, int& pos, int& pos2) {

  PatchMap* patchMap = PatchMap::Object();
  for (int ituple=0;ituple < ntuples;ituple++) {
    auto p0 = src[ituple].p[0];
    auto p1 = src[ituple].p[1];
    int pi0 = patchIndex[p0->patchID];
    int pi1 = patchIndex[p1->patchID];
    int l0 = src[ituple].localIndex[0];
    int l1 = src[ituple].localIndex[1];
    CompAtom& ca1 = p0->x[l0];
    CompAtom& ca2 = p1->x[l1];
    Position position1 = ca1.position;
    Position position2 = ca2.position;
    Vector shiftVec = lattice.wrap_delta_scaled(position1, position2);
    shiftVec += patchMap->center(p0->patchID) - patchMap->center(p1->patchID);
    CudaExclusion ce;
    ce.i            = l0 + patches[pi0].atomStart;
    ce.j            = l1 + patches[pi1].atomStart;
    ce.vdwtypei    = ca1.vdwType;
    ce.vdwtypej    = ca2.vdwType;
    ce.ioffsetXYZ = make_float3((float)shiftVec.x, (float)shiftVec.y, (float)shiftVec.z);
    //
    if (src[ituple].modified) {
      *dst1 = ce;
      dst1++;
      pos += typeSize;
    } else {
      *dst2 = ce;
      dst2++;
      pos2 += typeSize;
    }
  }
}

void ComputeBondedCUDA::copyCrosstermData(const int ntuples, const CrosstermElem* __restrict__ src,
  const CrosstermValue* __restrict__ crossterm_array, CudaCrossterm* __restrict__ dst) {

  PatchMap* patchMap = PatchMap::Object();
  for (int ituple=0;ituple < ntuples;ituple++) {
    auto p0 = src[ituple].p[0];
    auto p1 = src[ituple].p[1];
    auto p2 = src[ituple].p[2];
    auto p3 = src[ituple].p[3];
    auto p4 = src[ituple].p[4];
    auto p5 = src[ituple].p[5];
    auto p6 = src[ituple].p[6];
    auto p7 = src[ituple].p[7];
    int pi0 = patchIndex[p0->patchID];
    int pi1 = patchIndex[p1->patchID];
    int pi2 = patchIndex[p2->patchID];
    int pi3 = patchIndex[p3->patchID];
    int pi4 = patchIndex[p4->patchID];
    int pi5 = patchIndex[p5->patchID];
    int pi6 = patchIndex[p6->patchID];
    int pi7 = patchIndex[p7->patchID];
    int l0 = src[ituple].localIndex[0];
    int l1 = src[ituple].localIndex[1];
    int l2 = src[ituple].localIndex[2];
    int l3 = src[ituple].localIndex[3];
    int l4 = src[ituple].localIndex[4];
    int l5 = src[ituple].localIndex[5];
    int l6 = src[ituple].localIndex[6];
    int l7 = src[ituple].localIndex[7];
    dst[ituple].i1 = l0 + patches[pi0].atomStart;
    dst[ituple].i2 = l1 + patches[pi1].atomStart;
    dst[ituple].i3 = l2 + patches[pi2].atomStart;
    dst[ituple].i4 = l3 + patches[pi3].atomStart;
    dst[ituple].i5 = l4 + patches[pi4].atomStart;
    dst[ituple].i6 = l5 + patches[pi5].atomStart;
    dst[ituple].i7 = l6 + patches[pi6].atomStart;
    dst[ituple].i8 = l7 + patches[pi7].atomStart;
    dst[ituple].itype = (src[ituple].value - crossterm_array);
    Position position1 = p0->x[l0].position;
    Position position2 = p1->x[l1].position;
    Position position3 = p2->x[l2].position;
    Position position4 = p3->x[l3].position;
    Position position5 = p4->x[l4].position;
    Position position6 = p5->x[l5].position;
    Position position7 = p6->x[l6].position;
    Position position8 = p7->x[l7].position;
    Vector shiftVec12 = lattice.wrap_delta_scaled(position1, position2);
    Vector shiftVec23 = lattice.wrap_delta_scaled(position2, position3);
    Vector shiftVec34 = lattice.wrap_delta_scaled(position3, position4);
    Vector shiftVec56 = lattice.wrap_delta_scaled(position5, position6);
    Vector shiftVec67 = lattice.wrap_delta_scaled(position6, position7);
    Vector shiftVec78 = lattice.wrap_delta_scaled(position7, position8);
    shiftVec12 += patchMap->center(p0->patchID) - patchMap->center(p1->patchID);
    shiftVec23 += patchMap->center(p1->patchID) - patchMap->center(p2->patchID);
    shiftVec34 += patchMap->center(p2->patchID) - patchMap->center(p3->patchID);
    shiftVec56 += patchMap->center(p4->patchID) - patchMap->center(p5->patchID);
    shiftVec67 += patchMap->center(p5->patchID) - patchMap->center(p6->patchID);
    shiftVec78 += patchMap->center(p6->patchID) - patchMap->center(p7->patchID);
    dst[ituple].offset12XYZ = make_float3( (float)shiftVec12.x, (float)shiftVec12.y, (float)shiftVec12.z);
    dst[ituple].offset23XYZ = make_float3( (float)shiftVec23.x, (float)shiftVec23.y, (float)shiftVec23.z);
    dst[ituple].offset34XYZ = make_float3( (float)shiftVec34.x, (float)shiftVec34.y, (float)shiftVec34.z);
    dst[ituple].offset56XYZ = make_float3( (float)shiftVec56.x, (float)shiftVec56.y, (float)shiftVec56.z);
    dst[ituple].offset67XYZ = make_float3( (float)shiftVec67.x, (float)shiftVec67.y, (float)shiftVec67.z);
    dst[ituple].offset78XYZ = make_float3( (float)shiftVec78.x, (float)shiftVec78.y, (float)shiftVec78.z);
    dst[ituple].scale = src[ituple].scale;
  }
}

void ComputeBondedCUDA::tupleCopyWorker(int first, int last, void *result, int paraNum, void *param) {
  ComputeBondedCUDA* c = (ComputeBondedCUDA *)param;
  c->tupleCopyWorker(first, last);
}

void ComputeBondedCUDA::tupleCopyWorker(int first, int last) {
  if (first == -1) {
    // Separate exclusions into modified, and non-modified
    int pos = exclusionStartPos;
    int pos2 = exclusionStartPos2;
    for (auto it = tupleList[Tuples::EXCLUSION].begin();it != tupleList[Tuples::EXCLUSION].end();it++) {
      int ntuples = (*it)->getNumTuples();
      copyExclusionData(ntuples, (ExclElem *)(*it)->getTupleList(), CudaTupleTypeSize[Tuples::EXCLUSION],
        (CudaExclusion *)&tupleData[pos], (CudaExclusion *)&tupleData[pos2], pos, pos2);
    }
    first = 0;
  }
  for (int i=first;i <= last;i++) {
    switch (tupleCopyWorkList[i].tupletype) {

      case Tuples::BOND:
      {
        copyBondData(tupleCopyWorkList[i].ntuples, (BondElem *)tupleCopyWorkList[i].tupleElemList,
          Node::Object()->parameters->bond_array, (CudaBond *)&tupleData[tupleCopyWorkList[i].tupleDataPos]);
      }
      break;

      case Tuples::ANGLE:
      {
        copyAngleData(tupleCopyWorkList[i].ntuples, (AngleElem *)tupleCopyWorkList[i].tupleElemList,
          Node::Object()->parameters->angle_array, (CudaAngle *)&tupleData[tupleCopyWorkList[i].tupleDataPos]);
      }
      break;

      case Tuples::DIHEDRAL:
      {
        copyDihedralData<true, DihedralElem, DihedralValue>(tupleCopyWorkList[i].ntuples,
          (DihedralElem *)tupleCopyWorkList[i].tupleElemList, Node::Object()->parameters->dihedral_array,
          (CudaDihedral *)&tupleData[tupleCopyWorkList[i].tupleDataPos]);
      }
      break;

      case Tuples::IMPROPER:
      {
        copyDihedralData<false, ImproperElem, ImproperValue>(tupleCopyWorkList[i].ntuples,
          (ImproperElem *)tupleCopyWorkList[i].tupleElemList, Node::Object()->parameters->improper_array,
          (CudaDihedral *)&tupleData[tupleCopyWorkList[i].tupleDataPos]);
      }
      break;

      case Tuples::CROSSTERM:
      {
        copyCrosstermData(tupleCopyWorkList[i].ntuples, (CrosstermElem *)tupleCopyWorkList[i].tupleElemList,
          Node::Object()->parameters->crossterm_array, (CudaCrossterm *)&tupleData[tupleCopyWorkList[i].tupleDataPos]);
      }
      break;

      default:
      NAMD_bug("ComputeBondedCUDA::tupleCopyWorker, Unsupported tuple type");
      break;
    }
  }
}

//
// Copies tuple data form individual buffers to a single contigious buffer
// NOTE: This is done on the master PE
//
void ComputeBondedCUDA::copyTupleData() {

  PatchMap* patchMap = PatchMap::Object();

  // Count the number of exclusions
  int numModifiedExclusions = 0;
  int numExclusions = 0;
  for (int i=0;i < numExclPerRank.size();i++) {
    numModifiedExclusions += numExclPerRank[i].numModifiedExclusions;
    numExclusions         += numExclPerRank[i].numExclusions;
  }
  int numModifiedExclusionsWA = ComputeBondedCUDAKernel::warpAlign(numModifiedExclusions);
  int numExclusionsWA         = ComputeBondedCUDAKernel::warpAlign(numExclusions);

  // Count the number of tuples for each type
  int posWA = 0;
  exclusionStartPos = 0;
  exclusionStartPos2 = 0;
  tupleCopyWorkList.clear();
  for (int tupletype=0;tupletype < Tuples::NUM_TUPLE_TYPES;tupletype++) {
    // Take temporary position
    int pos = posWA;
    if (tupletype == Tuples::EXCLUSION) {
      exclusionStartPos = pos;
      exclusionStartPos2 = pos + numModifiedExclusionsWA*CudaTupleTypeSize[Tuples::EXCLUSION];
    }
    // Count for total number of tuples for this tupletype
    int num = 0;
    for (auto it = tupleList[tupletype].begin();it != tupleList[tupletype].end();it++) {
      int ntuples = (*it)->getNumTuples();
      num += ntuples;
      if (tupletype != Tuples::EXCLUSION) {
        TupleCopyWork tupleCopyWork;
        tupleCopyWork.tupletype     = tupletype;
        tupleCopyWork.ntuples       = ntuples;
        tupleCopyWork.tupleElemList = (*it)->getTupleList();
        tupleCopyWork.tupleDataPos  = pos;
        tupleCopyWorkList.push_back(tupleCopyWork);
        pos += ntuples*CudaTupleTypeSize[tupletype];
      }
    }
    numTuplesPerType[tupletype] = num;
    //
    if (tupletype == Tuples::EXCLUSION) {
      // Warp-align exclusions separately
      posWA += (numModifiedExclusionsWA + numExclusionsWA)*CudaTupleTypeSize[tupletype];
    } else {
      posWA += ComputeBondedCUDAKernel::warpAlign(num)*CudaTupleTypeSize[tupletype];
    }
  }
  if (numModifiedExclusions + numExclusions != numTuplesPerType[Tuples::EXCLUSION]) {
    NAMD_bug("ComputeBondedCUDA::copyTupleData, invalid number of exclusions");
  }

  // Set flags for finishPatchesOnPe
  hasExclusions = (numExclusions > 0);
  hasModifiedExclusions = (numModifiedExclusions > 0);

  // Re-allocate storage as needed
  // reallocate_host<char>(&tupleData, &tupleDataSize, size, 1.2f);
  reallocate_host<char>(&tupleData, &tupleDataSize, posWA, 1.2f);

#if CMK_SMP && USE_CKLOOP
  int useCkLoop = Node::Object()->simParameters->useCkLoop;
  if (useCkLoop >= 1) {
    CkLoop_Parallelize(tupleCopyWorker, 1, (void *)this, CkMyNodeSize(), -1, tupleCopyWorkList.size() - 1);
  } else
#endif
  {
    tupleCopyWorker(-1, tupleCopyWorkList.size() - 1);
  }

  bondedKernel.update(numTuplesPerType[Tuples::BOND], numTuplesPerType[Tuples::ANGLE],
    numTuplesPerType[Tuples::DIHEDRAL], numTuplesPerType[Tuples::IMPROPER],
    numModifiedExclusions, numExclusions, numTuplesPerType[Tuples::CROSSTERM],
    tupleData, stream);

  // Re-allocate forces
  int forceStorageSize = bondedKernel.getAllForceSize(atomStorageSize, true);
  reallocate_host<FORCE_TYPE>(&forces, &forcesSize, forceStorageSize, 1.4f);
}

//
// Launch work on GPU
//
void ComputeBondedCUDA::launchWork() {
  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::launchWork() called on non master PE");

  cudaCheck(cudaSetDevice(deviceID));

  if (atomsChanged) {
    copyTupleData();
  }

  float3 lata = make_float3(lattice.a().x, lattice.a().y, lattice.a().z);
  float3 latb = make_float3(lattice.b().x, lattice.b().y, lattice.b().z);
  float3 latc = make_float3(lattice.c().x, lattice.c().y, lattice.c().z);

  int r2_delta_expc = 64 * (ComputeNonbondedUtil::r2_delta_exp - 127);

  // Calculate forces
  bondedKernel.bondedForce(
    ComputeNonbondedUtil::scale14,
    atomStorageSize,
    doEnergy, doVirial, doSlow,
    lata, latb, latc,
    (float)ComputeNonbondedUtil::cutoff2,
    (float)ComputeNonbondedUtil::r2_delta, r2_delta_expc,
    (const float4*)atoms, forces,
    energies_virials,
    stream);

  forceDoneSetCallback();
}

void ComputeBondedCUDA::forceDoneCheck(void *arg, double walltime) {
  ComputeBondedCUDA* c = (ComputeBondedCUDA *)arg;

  if (CkMyPe() != c->masterPe)
    NAMD_bug("ComputeBondedCUDA::forceDoneCheck called on non masterPe");

  cudaCheck(cudaSetDevice(c->deviceID));

  cudaError_t err = cudaEventQuery(c->forceDoneEvent);
  if (err == cudaSuccess) {
    // Event has occurred
    c->checkCount = 0;
    traceUserBracketEvent(CUDA_BONDED_KERNEL_EVENT, c->beforeForceCompute, walltime);
    c->finishPatches();
    return;
  } else if (err != cudaErrorNotReady) {
    // Anything else is an error
    char errmsg[256];
    sprintf(errmsg,"in ComputeBondedCUDA::forceDoneCheck after polling %d times over %f s",
            c->checkCount, walltime - c->beforeForceCompute);
    cudaDie(errmsg,err);
  }

  // Event has not occurred
  c->checkCount++;
  if (c->checkCount >= 1000000) {
    char errmsg[256];
    sprintf(errmsg,"ComputeBondedCUDA::forceDoneCheck polled %d times over %f s",
            c->checkCount, walltime - c->beforeForceCompute);
    cudaDie(errmsg,err);
  }

  // Call again 
  CcdCallBacksReset(0, walltime);
  CcdCallFnAfter(forceDoneCheck, arg, 0.1);
}

//
// Set call back for all the work in the stream at this point
//
void ComputeBondedCUDA::forceDoneSetCallback() {
  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::forceDoneSetCallback called on non masterPe");
  cudaCheck(cudaSetDevice(deviceID));
  cudaCheck(cudaEventRecord(forceDoneEvent, stream));
  checkCount = 0;
  CcdCallBacksReset(0, CmiWallTimer());
  // Start timer for CUDA kernel
  beforeForceCompute = CkWallTimer();
  // Set the call back at 0.1ms
  CcdCallFnAfter(forceDoneCheck, this, 0.1);
}

inline void convertForceToDouble(const FORCE_TYPE *af, const int forceStride, double& afx, double& afy, double& afz) {
#ifdef USE_STRIDED_FORCE
  FORCE_TYPE afxt = af[0];
  FORCE_TYPE afyt = af[forceStride];
  FORCE_TYPE afzt = af[forceStride*2];
#else
  FORCE_TYPE afxt = af->x;
  FORCE_TYPE afyt = af->y;
  FORCE_TYPE afzt = af->z;
#endif
#ifdef USE_FP_FORCE
  afx = ((double)afxt)*force_to_double;
  afy = ((double)afyt)*force_to_double;
  afz = ((double)afzt)*force_to_double;
#else
  afx = afxt;
  afy = afyt;
  afz = afzt;
#endif
}

template <bool sumNbond, bool sumSlow>
void finishForceLoop(const int numAtoms, const int forceStride,
  const FORCE_TYPE* __restrict__ af, const FORCE_TYPE* __restrict__ af_nbond, const FORCE_TYPE* __restrict__ af_slow, 
  Force* __restrict__ f, Force* __restrict__ f_nbond, Force* __restrict__ f_slow) {

  for (int j=0;j < numAtoms;j++) {
    {
      double afx, afy, afz;
      convertForceToDouble(af + j, forceStride, afx, afy, afz);
      f[j].x += afx;
      f[j].y += afy;
      f[j].z += afz;
    }
    if (sumNbond)
    {
      double afx, afy, afz;
      convertForceToDouble(af_nbond + j, forceStride, afx, afy, afz);
      f_nbond[j].x += afx;
      f_nbond[j].y += afy;
      f_nbond[j].z += afz;
    }
    if (sumSlow)
    {
      double afx, afy, afz;
      convertForceToDouble(af_slow + j, forceStride, afx, afy, afz);
      f_slow[j].x += afx;
      f_slow[j].y += afy;
      f_slow[j].z += afz;
    }
  }

}

//
// Finish all patches that are on this pe
//
void ComputeBondedCUDA::finishPatchesOnPe() {

  PatchMap* patchMap = PatchMap::Object();
  int myRank = CkMyRank();

  const int forceStride = bondedKernel.getForceStride(atomStorageSize);
  const int forceSize = bondedKernel.getForceSize(atomStorageSize);
  const bool sumNbond = hasModifiedExclusions;
  const bool sumSlow = (hasModifiedExclusions || hasExclusions) && doSlow;

  for (int i=0;i < patchIDsPerRank[myRank].size();i++) {
    PatchID patchID = patchIDsPerRank[myRank][i];
    Patch* patch = patchMap->patch(patchID);
    TuplePatchElem* tpe = tuplePatchList.find(TuplePatchElem(patchID));
    if (tpe == NULL) {
      NAMD_bug("ComputeBondedCUDA::finishPatchesOnPe, TuplePatchElem not found");
    }

    int pi = patchIndex[patchID];
    int numAtoms = patches[pi].numAtoms;
    int atomStart = patches[pi].atomStart;

    Force *f = tpe->f;
    Force *f_nbond = tpe->r->f[Results::nbond];
    Force *f_slow = tpe->r->f[Results::slow];

    FORCE_TYPE *af       = forces + atomStart;
    FORCE_TYPE *af_nbond = forces + forceSize + atomStart;
    FORCE_TYPE *af_slow  = forces + 2*forceSize + atomStart;

    if (!sumNbond && !sumSlow) {
      finishForceLoop<false, false>(numAtoms, forceStride, af, af_nbond, af_slow, f, f_nbond, f_slow);
    } else if (sumNbond && !sumSlow) {
      finishForceLoop<true, false>(numAtoms, forceStride, af, af_nbond, af_slow, f, f_nbond, f_slow);
    } else if (!sumNbond && sumSlow) {
      finishForceLoop<false, true>(numAtoms, forceStride, af, af_nbond, af_slow, f, f_nbond, f_slow);
    } else if (sumNbond && sumSlow) {
      finishForceLoop<true, true>(numAtoms, forceStride, af, af_nbond, af_slow, f, f_nbond, f_slow);
    } else {
      NAMD_bug("ComputeBondedCUDA::finishPatchesOnPe, logically impossible choice");
    }

    // for (int j=0;j < numAtoms;j++) {
    //   {
    //     double afx, afy, afz;
    //     convertForceToDouble(af + j, forceStride, afx, afy, afz);
    //     f[j].x += afx;
    //     f[j].y += afy;
    //     f[j].z += afz;
    //   }
    //   if (sumNbond)
    //   {
    //     double afx, afy, afz;
    //     convertForceToDouble(af_nbond + j, forceStride, afx, afy, afz);
    //     f_nbond[j].x += afx;
    //     f_nbond[j].y += afy;
    //     f_nbond[j].z += afz;
    //   }
    //   if (sumSlow)
    //   {
    //     double afx, afy, afz;
    //     convertForceToDouble(af_slow + j, forceStride, afx, afy, afz);
    //     f_slow[j].x += afx;
    //     f_slow[j].y += afy;
    //     f_slow[j].z += afz;
    //   }
    // }

    tpe->forceBox->close(&tpe->r);
    tpe->positionBox->close(&tpe->x);
    if ( doMolly ) tpe->avgPositionBox->close(&tpe->x_avg);
  }

  bool done = false;
  CmiLock(lock);
  patchesCounter -= patchIDsPerRank[CkMyRank()].size();
  if (patchesCounter == 0) {
    patchesCounter = getNumPatches();
    done = true;
  }
  CmiUnlock(lock);
  if (done) {
    computeMgr->sendFinishReductions(masterPe, this);
  }

}

void ComputeBondedCUDA::finishPatches() {

  if (atomsChanged) {
    unmapAtoms();
  }

  computeMgr->sendFinishPatchesOnPe(pes, this);
}

#ifdef WRITE_FULL_VIRIALS
#ifdef USE_FP_VIRIAL
void convertVirial(double *virial) {
  long long int *virial_lli = (long long int *)virial;
  for (int i=0;i < 9;i++) {
    virial[i] = ((double)virial_lli[i])*virial_to_double;
  }
}
#endif
#endif

//
// Finish & submit reductions
//
void ComputeBondedCUDA::finishReductions() {

  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::finishReductions() called on non masterPe");

  // static int ncall = 0;
  // ncall++;

  int pos = 0;
  for (int tupletype=0;tupletype < Tuples::NUM_TUPLE_TYPES;tupletype++) {
    if (numTuplesPerType[tupletype] > 0) {

      if (doEnergy) {
        switch (tupletype) {
          case Tuples::BOND:
          reduction->item(REDUCTION_BOND_ENERGY) += energies_virials[ComputeBondedCUDAKernel::energyIndex_BOND];
          break;

          case Tuples::ANGLE:
          reduction->item(REDUCTION_ANGLE_ENERGY) += energies_virials[ComputeBondedCUDAKernel::energyIndex_ANGLE];
          break;

          case Tuples::DIHEDRAL:
          reduction->item(REDUCTION_DIHEDRAL_ENERGY) += energies_virials[ComputeBondedCUDAKernel::energyIndex_DIHEDRAL];
          break;

          case Tuples::IMPROPER:
          reduction->item(REDUCTION_IMPROPER_ENERGY) += energies_virials[ComputeBondedCUDAKernel::energyIndex_IMPROPER];
          break;

          case Tuples::EXCLUSION:
          reduction->item(REDUCTION_ELECT_ENERGY)      += energies_virials[ComputeBondedCUDAKernel::energyIndex_ELECT];
          reduction->item(REDUCTION_LJ_ENERGY)         += energies_virials[ComputeBondedCUDAKernel::energyIndex_LJ];
          reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += energies_virials[ComputeBondedCUDAKernel::energyIndex_ELECT_SLOW];
          break;

          case Tuples::CROSSTERM:
          reduction->item(REDUCTION_CROSSTERM_ENERGY) += energies_virials[ComputeBondedCUDAKernel::energyIndex_CROSSTERM];
          break;

          default:
          NAMD_bug("ComputeBondedCUDA::finishReductions, Unsupported tuple type");
          break;
        }
      }

      auto it = tupleList[tupletype].begin();
      (*it)->submitTupleCount(reduction, numTuplesPerType[tupletype]);
    }
  }

  if (doVirial) {
#ifdef WRITE_FULL_VIRIALS
#ifdef USE_FP_VIRIAL
    convertVirial(&energies_virials[ComputeBondedCUDAKernel::normalVirialIndex_XX]);
    convertVirial(&energies_virials[ComputeBondedCUDAKernel::nbondVirialIndex_XX]);
    convertVirial(&energies_virials[ComputeBondedCUDAKernel::slowVirialIndex_XX]);
    convertVirial(&energies_virials[ComputeBondedCUDAKernel::amdDiheVirialIndex_XX]);
#endif
#else
#error "non-WRITE_FULL_VIRIALS not implemented"
#endif
    ADD_TENSOR(reduction, REDUCTION_VIRIAL_NORMAL,energies_virials, ComputeBondedCUDAKernel::normalVirialIndex);
    ADD_TENSOR(reduction, REDUCTION_VIRIAL_NBOND, energies_virials, ComputeBondedCUDAKernel::nbondVirialIndex);
    ADD_TENSOR(reduction, REDUCTION_VIRIAL_SLOW,  energies_virials, ComputeBondedCUDAKernel::slowVirialIndex);
    ADD_TENSOR(reduction, REDUCTION_VIRIAL_AMD_DIHE, energies_virials, ComputeBondedCUDAKernel::amdDiheVirialIndex);
    // NOTE: AMD_DIHE virial is also added to NORMAL virial.
    // This is what happens in ComputeDihedrals.C and ComputeCrossterms.C
    ADD_TENSOR(reduction, REDUCTION_VIRIAL_NORMAL,   energies_virials, ComputeBondedCUDAKernel::amdDiheVirialIndex);
  }

  reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
  reduction->submit();
}

//
// Can only be called by master PE
//
void ComputeBondedCUDA::initialize() {

  if (CkMyPe() != masterPe)
    NAMD_bug("ComputeBondedCUDA::initialize() called on non master PE");

  // Build list of PEs
  for (int rank=0;rank < computes.size();rank++) {
    if (computes[rank].selfComputes.size() > 0 || computes[rank].homeCompute.patchIDs.size() > 0) {
      pes.push_back(CkNodeFirst(CkMyNode()) + rank);
    }
  }

  // Return if no work to do
  if (pes.size() == 0) return;

  initializeCalled = true;
  cudaCheck(cudaSetDevice(deviceID));

#if CUDA_VERSION >= 5050
  int leastPriority, greatestPriority;
  cudaCheck(cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority));
  cudaCheck(cudaStreamCreateWithPriority(&stream, cudaStreamDefault, greatestPriority));
#else
  cudaCheck(cudaStreamCreate(&stream));
#endif
  cudaCheck(cudaEventCreate(&forceDoneEvent));
  lock = CmiCreateLock();

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  PatchMap* patchMap = PatchMap::Object();

  // First, assign all patches in self computes.
  // NOTE: These never overlap between PEs. No proxies added.
  for (int rank=0;rank < computes.size();rank++) {
    std::vector< SelfCompute >& selfComputes = computes[rank].selfComputes;
    for (auto it=selfComputes.begin();it != selfComputes.end();it++) {
      for (auto jt=it->patchIDs.begin();jt != it->patchIDs.end();jt++) {
        if (!tuplePatchList.find( TuplePatchElem(*jt) ) ) {
          tuplePatchList.add( TuplePatchElem(*jt) );
          patchIDsPerRank[rank].push_back(*jt);
          allPatchIDs.push_back(*jt);
        }
      }
    }
  }

  // Second, assign all patches in home computes.
  // NOTE: The ranks always have these patches. No proxies added.
  for (int rank=0;rank < computes.size();rank++) {
    HomeCompute& homeCompute = computes[rank].homeCompute;
    std::vector<int>& patchIDs = homeCompute.patchIDs;
    for (int i=0;i < patchIDs.size();i++) {
      int patchID = patchIDs[i];
      if (!tuplePatchList.find( TuplePatchElem(patchID) ) ) {
        tuplePatchList.add( TuplePatchElem(patchID) );
        patchIDsPerRank[rank].push_back(patchID);
        allPatchIDs.push_back(patchID);
      }
    }
  }

  std::vector< std::vector<int> > patchIDsToAppend(CkMyNodeSize());
  // Find neighbors that are not added yet
  std::vector<int> neighborPids;
  for (int rank=0;rank < computes.size();rank++) {
    PatchID neighbors[PatchMap::MaxOneOrTwoAway];
    HomeCompute& homeCompute = computes[rank].homeCompute;
    std::vector<int>& patchIDs = homeCompute.patchIDs;
    for (int i=0;i < patchIDs.size();i++) {
      int patchID = patchIDs[i];
      int numNeighbors = patchMap->upstreamNeighbors(patchID, neighbors);
      for (int j=0;j < numNeighbors;j++) {
        if (!tuplePatchList.find( TuplePatchElem(neighbors[j]) ) ) {
          neighborPids.push_back(neighbors[j]);
        }
      }
    }
  }
  // Remove duplicates from neighborPids
  {
    std::sort(neighborPids.begin(), neighborPids.end());
    auto it_end = std::unique(neighborPids.begin(), neighborPids.end());
    neighborPids.resize(std::distance(neighborPids.begin(), it_end));
  }
  // Assign neighbors to the PEs on this node that have them
  for (int i=0;i < neighborPids.size();i++) {
    for (int rank=0;rank < computes.size();rank++) {
      int pid = neighborPids[i];
      int pe = rank + CkNodeFirst(CkMyNode());
      if (patchMap->node(pid) == pe) {
        // Patch pid found on PE "pe" on this node
        tuplePatchList.add( TuplePatchElem(pid) );
        patchIDsPerRank[rank].push_back(pid);
        allPatchIDs.push_back(pid);
        // Add to this rank's patches
        patchIDsToAppend[rank].push_back(pid);
        // Add to the list of PEs
        pes.push_back(CkNodeFirst(CkMyNode()) + rank);
        break;
      }
    }
  }
  // Remove duplicates from pes
  {
    std::sort(pes.begin(), pes.end());
    auto it_end = std::unique(pes.begin(), pes.end());
    pes.resize(std::distance(pes.begin(), it_end));
  }
  
  // Last, assign all patches in neighbors of home computes
  // NOTE: Will create proxies on multiple nodes
  for (int rank=0;rank < computes.size();rank++) {
    PatchID neighbors[PatchMap::MaxOneOrTwoAway];
    HomeCompute& homeCompute = computes[rank].homeCompute;
    std::vector<int>& patchIDs = homeCompute.patchIDs;
    std::vector<int> neighborPatchIDs;
    for (int i=0;i < patchIDs.size();i++) {
      int patchID = patchIDs[i];
      int numNeighbors = patchMap->upstreamNeighbors(patchID, neighbors);
      for (int j=0;j < numNeighbors;j++) {
        if (!tuplePatchList.find( TuplePatchElem(neighbors[j]) ) ) {
          // Patch not found => Add Proxy
          tuplePatchList.add( TuplePatchElem(neighbors[j]) );
          patchIDsPerRank[rank].push_back(neighbors[j]);
          allPatchIDs.push_back(neighbors[j]);
        }
        if ( std::count(patchIDs.begin(), patchIDs.end(), neighbors[j]) == 0 
          && std::count(neighborPatchIDs.begin(), neighborPatchIDs.end(), neighbors[j]) == 0 ) {
          neighborPatchIDs.push_back(neighbors[j]);
        }
      }
    }
    // Append neighboring patchIDs to homeCompute.patchIDs
    // int start = patchIDs.size();
    // patchIDs.resize(patchIDs.size() + neighborPatchIDs.size());
    // for (int i=0;i < neighborPatchIDs.size();i++) {
    //   patchIDs[start + i] = neighborPatchIDs[i];
    // }
    for (int i=0;i < neighborPatchIDs.size();i++) {
      patchIDsToAppend[rank].push_back(neighborPatchIDs[i]);
    }
  }

  for (int rank=0;rank < patchIDsToAppend.size();rank++) {
    for (int i=0;i < patchIDsToAppend[rank].size();i++) {
      computes[rank].homeCompute.patchIDs.push_back(patchIDsToAppend[rank][i]);
    }
  }

  // Remove duplicate patch IDs
  {
    std::sort(allPatchIDs.begin(), allPatchIDs.end());
    auto it_end = std::unique(allPatchIDs.begin(), allPatchIDs.end());
    allPatchIDs.resize(std::distance(allPatchIDs.begin(), it_end));
  }

  // Set number of (unique) patches
  setNumPatches(allPatchIDs.size());

  // Reset patchesCounter
  patchesCounter = getNumPatches();

  patches.resize(getNumPatches());

  // Setup tupleList
  for (int rank=0;rank < computes.size();rank++) {
    std::vector< SelfCompute >& selfComputes = computes[rank].selfComputes;
    for (auto it=selfComputes.begin();it != selfComputes.end();it++) {
      tupleList[it->tuples->getType()].push_back(it->tuples);
    }
    HomeCompute& homeCompute = computes[rank].homeCompute;
    for (int i=0;i < homeCompute.tuples.size();i++) {
      tupleList[homeCompute.tuples[i]->getType()].push_back(homeCompute.tuples[i]);
    }
  }  

  // Allocate host memory for energies and virials
  allocate_host<double>(&energies_virials, ComputeBondedCUDAKernel::energies_virials_SIZE);

  // Finally, do sanity checks
  std::vector<char> patchIDset(patchMap->numPatches(), 0);
  int numPatchIDset = 0;
  int numPatchIDs = 0;
  for (int rank=0;rank < computes.size();rank++) {
    numPatchIDs += patchIDsPerRank[rank].size();
    for (int i=0;i < patchIDsPerRank[rank].size();i++) {
      PatchID patchID = patchIDsPerRank[rank][i];
      if (patchIDset[patchID] == 0) numPatchIDset++;
      patchIDset[patchID] = 1;
      if ( !std::count(allPatchIDs.begin(), allPatchIDs.end(), patchID) ) {
        NAMD_bug("ComputeBondedCUDA::initialize(), inconsistent patch mapping");
      }
    }
  }
  if (numPatchIDs != getNumPatches() || numPatchIDset != getNumPatches()) {
    NAMD_bug("ComputeBondedCUDA::initialize(), inconsistent patch mapping");
  }

  // Warning: Direct indexing used, patchIndex could use up a lot of memory for large systems
  patchIndex.resize(patchMap->numPatches());
  atomMappers.resize(getNumPatches());
  for (int i=0;i < getNumPatches();i++) {
    atomMappers[i] = new AtomMapper(allPatchIDs[i], &atomMap);
    patchIndex[allPatchIDs[i]] = i;
  }

  // Copy coefficients to GPU
  Parameters* parameters = Node::Object()->parameters;
  for (int tupletype=0;tupletype < Tuples::NUM_TUPLE_TYPES;tupletype++) {
    if (tupleList[tupletype].size() > 0) {
      switch(tupletype) {

        case Tuples::BOND:
        {
          int NumBondParams = parameters->NumBondParams;
          BondValue* bond_array = parameters->bond_array;
          std::vector<CudaBondValue> bondValues(NumBondParams);
          for (int i=0;i < NumBondParams;i++) {
            bondValues[i].k  = bond_array[i].k;
            bondValues[i].x0 = bond_array[i].x0;
            bondValues[i].x1 = bond_array[i].x1;
          }
          bondedKernel.setupBondValues(NumBondParams, bondValues.data());
        }
        break;

        case Tuples::ANGLE:
        {
          int NumAngleParams = parameters->NumAngleParams;
          AngleValue* angle_array = parameters->angle_array;
          std::vector<CudaAngleValue> angleValues(NumAngleParams);
          bool normal_ub_error = false;
          for (int i=0;i < NumAngleParams;i++) {
            angleValues[i].k      = angle_array[i].k;
            if (angle_array[i].normal == 1) {
              angleValues[i].theta0 = angle_array[i].theta0;
            } else {
              angleValues[i].theta0 = cos(angle_array[i].theta0);
            }
            normal_ub_error |= (angle_array[i].normal == 0 && angle_array[i].k_ub);
            angleValues[i].k_ub   = angle_array[i].k_ub;
            angleValues[i].r_ub   = angle_array[i].r_ub;
            angleValues[i].normal = angle_array[i].normal;
          }
          if (normal_ub_error) NAMD_die("ERROR: Can't use cosAngles with Urey-Bradley angles");
          bondedKernel.setupAngleValues(NumAngleParams, angleValues.data());
        }
        break;

        case Tuples::DIHEDRAL:
        {
          int NumDihedralParams = parameters->NumDihedralParams;
          DihedralValue* dihedral_array = parameters->dihedral_array;
          int NumDihedralParamsMult = 0;
          for (int i=0;i < NumDihedralParams;i++) {
            NumDihedralParamsMult += std::max(0, dihedral_array[i].multiplicity);
          }
          std::vector<CudaDihedralValue> dihedralValues(NumDihedralParamsMult);
          dihedralMultMap.resize(NumDihedralParams);
          int k = 0;
          for (int i=0;i < NumDihedralParams;i++) {
            int multiplicity = dihedral_array[i].multiplicity;
            dihedralMultMap[i] = k;
            for (int j=0;j < multiplicity;j++) {
              dihedralValues[k].k     = dihedral_array[i].values[j].k;
              dihedralValues[k].n     = (dihedral_array[i].values[j].n << 1) | (j < (multiplicity - 1));
              dihedralValues[k].delta = dihedral_array[i].values[j].delta;
              k++;
            }
          }
          bondedKernel.setupDihedralValues(NumDihedralParamsMult, dihedralValues.data());
        }
        break;

        case Tuples::IMPROPER:
        {
          int NumImproperParams = parameters->NumImproperParams;
          ImproperValue* improper_array = parameters->improper_array;
          int NumImproperParamsMult = 0;
          for (int i=0;i < NumImproperParams;i++) {
            NumImproperParamsMult += std::max(0, improper_array[i].multiplicity);
          }
          std::vector<CudaDihedralValue> improperValues(NumImproperParamsMult);
          improperMultMap.resize(NumImproperParams);
          int k = 0;
          for (int i=0;i < NumImproperParams;i++) {
            int multiplicity = improper_array[i].multiplicity;
            improperMultMap[i] = k;
            for (int j=0;j < multiplicity;j++) {
              improperValues[k].k     = improper_array[i].values[j].k;
              improperValues[k].n     = (improper_array[i].values[j].n << 1) | (j < (multiplicity - 1));
              improperValues[k].delta = improper_array[i].values[j].delta;
              k++;
            }
          }
          bondedKernel.setupImproperValues(NumImproperParamsMult, improperValues.data());
        }
        break;

        case Tuples::CROSSTERM:
        {
          int NumCrosstermParams = parameters->NumCrosstermParams;
          CrosstermValue* crossterm_array = parameters->crossterm_array;
          std::vector<CudaCrosstermValue> crosstermValues(NumCrosstermParams);
          const int D = CrosstermValue::dim;
          const int N = CrosstermValue::dim - 1;
          for (int ipar=0;ipar < NumCrosstermParams;ipar++) {
            for (int i=0;i < N;i++) {
              for (int j=0;j < N;j++) {

                // Setups coefficients for bi-cubic interpolation.
                // See https://en.wikipedia.org/wiki/Bicubic_interpolation

                #define INDEX(ncols,i,j)  ((i)*ncols + (j))
                CrosstermData* table = &crossterm_array[ipar].c[0][0];

                const double Ainv[16][16] = {
                  { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
                  { 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
                  {-3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
                  { 2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
                  { 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
                  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
                  { 0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0},
                  { 0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0},
                  {-3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0},
                  { 0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0},
                  { 9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1},
                  {-6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1},
                  { 2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0},
                  { 0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0},
                  {-6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1},
                  { 4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1}
                };

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

                const double h = M_PI/12.0;

                const double x[16] = {
                  table[INDEX(D,i,j)].d00, table[INDEX(D,i+1,j)].d00, table[INDEX(D,i,j+1)].d00, table[INDEX(D,i+1,j+1)].d00,
                  table[INDEX(D,i,j)].d10*h, table[INDEX(D,i+1,j)].d10*h, table[INDEX(D,i,j+1)].d10*h, table[INDEX(D,i+1,j+1)].d10*h,
                  table[INDEX(D,i,j)].d01*h, table[INDEX(D,i+1,j)].d01*h, table[INDEX(D,i,j+1)].d01*h, table[INDEX(D,i+1,j+1)].d01*h,
                  table[INDEX(D,i,j)].d11*h*h, table[INDEX(D,i+1,j)].d11*h*h, table[INDEX(D,i,j+1)].d11*h*h, table[INDEX(D,i+1,j+1)].d11*h*h
                };

                // a = Ainv*x
                float* a = (float *)&crosstermValues[ipar].c[i][j][0];
                for (int k=0;k < 16;k++) {
                  double a_val = 0.0;
                  for (int l=0;l < 16;l++) {
                    a_val += Ainv[k][l]*x[l];
                  }
                  a[k] = (float)a_val;
                }

              }
            }
          }
          bondedKernel.setupCrosstermValues(NumCrosstermParams, crosstermValues.data());
        }
        break;

        case Tuples::EXCLUSION:
        // Nothing to do
        break;

        default:
        NAMD_bug("ComputeBondedCUDA::initialize, Undefined tuple type");
        break;
      }
    }
  }

  computeMgr->sendAssignPatchesOnPe(pes, this);
}

#endif // BONDED_CUDA
#endif // NAMD_CUDA
