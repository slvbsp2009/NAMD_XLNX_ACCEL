/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/LdbCoordinator.h,v $
 * $Author: jim $
 * $Date: 2013/09/06 19:11:37 $
 * $Revision: 1.47 $
 *****************************************************************************/

#ifndef LDBCOORDINATOR_H
#define LDBCOORDINATOR_H

#include <stdio.h>

#include <charm++.h>

#ifdef LB_MANAGER_VERSION
#include <LBManager.h>
typedef LBManager LdbInfra;
#else
#include <LBDatabase.h>
typedef LBDatabase LdbInfra;
#endif

#include "NamdTypes.h"
#include "BOCgroup.h"
#include "LdbCoordinator.decl.h"

/// In the new 64-bit id case defining LdbId as CmiUInt8,
/// the first 32 bits store the object's index and
/// the second 32 bits store the type.
///
/// In the old int[4] id case defining LdbId as LDObjid,
/// element 0 stores the object's index, element 1 stores the type,
/// and elements 2 and 3 are unused.
#if CMK_LBID_64BIT
typedef CmiUInt8 LdbId;
#else
typedef LDObjid LdbId;
#endif

inline const int& LdbIdField(const LdbId& id, const int index) {
#if CMK_LBID_64BIT
  return *(((int*)&id) + index);
#else
  return id.id[index];
#endif
}

inline int& LdbIdField(LdbId& id, const int index) {
#if CMK_LBID_64BIT
  return *(((int*)&id) + index);
#else
  return id.id[index];
#endif
}

/// Define the types encoded into the load balancing id.
/// Use negative numbers because the nonbonded/self types
/// are represented with the leading patch ID for that compute,
/// when available.
enum {
  NONBONDED_OR_SELF_TYPE = -1,  ///< represents nonbonded or self compute
  PATCH_TYPE = -2,              ///< represents a patch
  BONDED_TYPE = -3              ///< represents bonded compute
};

class PatchMap;
class ComputeMap;
class Controller;
class Sequencer;
class computeInfo;
class patchInfo;
class processorInfo;

enum {LDB_PATCHES = 4096};
enum {LDB_COMPUTES = 16384};
enum {COMPUTEMAX = 16384};
enum {PATCHMAX = 4096};
enum {PROCESSORMAX = 512};

void LdbCoordinator_initproc();

class LdbCoordinator : public CBase_LdbCoordinator
{
public:
  LdbCoordinator();
  ~LdbCoordinator(void);
  static LdbCoordinator *Object()  { 
    return CkpvAccess(LdbCoordinator_instance); 
  }

  void initialize(PatchMap *pmap, ComputeMap *cmap, int reinit=0);
  void createLoadBalancer();
  void patchLoad(PatchID id, int nAtoms, int timestep);

  void startWork(const LDObjHandle &handle) {  // start timer
    theLbdb->ObjectStart(handle);
  }
  void pauseWork(const LDObjHandle &handle) {  // stop timer only
    theLbdb->ObjectStop(handle);
  }
  void skipWork(const LDObjHandle &handle) {  // increment counter only
    nComputesReported++;
  }
  void endWork(const LDObjHandle &handle) {  // both
    theLbdb->ObjectStop(handle);
    nComputesReported++;
  }

  void rebalance(Sequencer *seq, PatchID id);
  void rebalance(Controller *seq);
  void nodeDone(CkReductionMsg *);
  void updateComputesReady();
  void barrier(void);
  void resume(void);
  void resumeReady(CkQdMsg *msg);
  void resume2(void);
  int getNumStepsToRun(void) { return numStepsToRun; }
  static void staticMigrateFn(LDObjHandle handle, int dest);
  static void staticStatsFn(LDOMHandle h, int state);
  static void staticQueryEstLoadFn(LDOMHandle h);
  static void staticReceiveAtSync(void* data);
  static void staticResumeFromSync(void* data);
  void AtSyncBarrierReached(void);
  void Migrate(LDObjHandle handle, int dest);
  void RecvMigrate(LdbMigrateMsg*);
  void ExpectMigrate(LdbMigrateMsg*);
  void ResumeFromSync(void);

public:
  void ExecuteMigrations(void);
  void awakenSequencers(void);
  int requiredProxies(PatchID id, int []);
  void printRequiredProxies(PatchID id, FILE *fp);
  void printLocalLdbReport(void);

  int stepsPerLdbCycle;
  int nLocalComputes;
  int nLocalPatches;
  int nPatchesReported;
  int nPatchesExpected;
  int nComputesReported;
  int nComputesExpected;
  int controllerReported;
  int controllerExpected;
  int nStatsMessagesReceived;
  int nStatsMessagesExpected;
  ComputeMap *computeMap;
  PatchMap *patchMap;
  int *patchNAtoms;
  int  nPatches;
  Controller *controllerThread;
  Sequencer **sequencerThreads;

  int ldbCycleNum;
  int numStepsToRun;	// tells Controller how many time steps to run 
			// before another load balancing
  int firstLdbStep;
  int totalStepsDone;	// keeps a count of the total number of
			// time steps to stop load balancing
  int takingLdbData;

  FILE *ldbStatsFP;
  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;

  LdbInfra *theLbdb;
  LDBarrierClient ldBarrierHandle;
  LDOMid myOMid;
  LDOMHandle myHandle;
  LdbMigrateMsg *migrateMsgs;
  int numComputes;
  int nRegisteredObjs;
  int reg_all_objs;
  LDObjHandle* patchHandles;

  void sendCollectLoads(CollectLoadsMsg*);
  void collectLoads(CollectLoadsMsg*);
private:
  int collPes;
  int reverted;
  int initTotalProxies;
  int finalTotalProxies;
  int initMaxPeProxies;
  int finalMaxPeProxies;
  int initMaxPatchProxies;
  int finalMaxPatchProxies;
  double initTime;
  double finalTime;
  double initMemory;
  double finalMemory;
  double initAvgPeLoad;
  double finalAvgPeLoad;
  double initMaxPeLoad;
  double finalMaxPeLoad;
};

class CollectLoadsMsg : public CMessage_CollectLoadsMsg {
public:
  int firstPe;
  int lastPe;
  int reverted;
  int initTotalProxies;
  int finalTotalProxies;
  int initMaxPeProxies;
  int finalMaxPeProxies;
  int initMaxPatchProxies;
  int finalMaxPatchProxies;
  double initTime;
  double finalTime;
  double initMemory;
  double finalMemory;
  double initAvgPeLoad;
  double finalAvgPeLoad;
  double initMaxPeLoad;
  double finalMaxPeLoad;
  char strategyName[16];
};

class LdbMigrateMsg : public CMessage_LdbMigrateMsg
{
public:
  LDObjHandle handle;
  int from;
  int to;
  LdbMigrateMsg *next;
};


#endif // LDBCOORDINATOR_H

