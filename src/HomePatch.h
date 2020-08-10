/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   HomePatch is the key distributed source/sink of Atom data
   including positions, velocities and forces applied
*/

#ifndef HOMEPATCH_H
#define HOMEPATCH_H

#include "charm++.h"

#include "NamdTypes.h"
#include "Patch.h"
#include "PatchMap.h"

#include "MigrateAtomsMsg.h"
#include "main.h"
#include "common.h"
#include "Migration.h"
#include "Settle.h"

#include <string>
#include <map>

class RegisterProxyMsg;
class UnregisterProxyMsg;
class ProxyResultVarsizeMsg;
class ProxyResultMsg;
class ProxyCombinedResultRawMsg;
class Sequencer;
class SubmitReduction;
class ProxyGBISP1ResultMsg;
class ProxyGBISP2ResultMsg;
class CheckpointAtomsMsg;
class ExchangeAtomsMsg;

class ProxyNodeAwareSpanningTreeMsg;

class ComputeQMMgr;


#ifdef TIMER_COLLECTION

#include <time.h>

struct TimerMicrosecond {
  struct timespec ts;
  inline void start() {
    clock_gettime(CLOCK_REALTIME, &ts);
  }
  inline double stop() {
    struct timespec tsend;
    clock_gettime(CLOCK_REALTIME, &tsend);
    return( (tsend.tv_sec - ts.tv_sec) * 1e6      // sec to microsec
        + (tsend.tv_nsec - ts.tv_nsec) * 1e-3 );  // nanosec to microsec
  }
};

#define TIMER_SLOTS  101
#define TIMER_SLOT_WIDTH  1

/// Timer entry
struct TimerEntry {
  TimerMicrosecond tmicro;
  double tcur;   ///< timer current entry
  double tavg;   ///< timer average
  double tvar;   ///< timer variance
  double tstd;   ///< timer standard deviation
  double tmin;   ///< timer minimum
  double tmax;   ///< timer maximum
  double tsum;   ///< timer summation
#if defined(DEBUG_TIMER_COLLECTION)
  double tcharm; ///< compare with timer from Charm++ for debugging
#endif
#if defined(TIMER_HISTOGRAM)
  double slotwidth;
  double inv_slotwidth;
  int hist[TIMER_SLOTS];
#endif
  int count;
  TimerEntry() { reset(); }
  /// Reset accumulators.
  inline void reset() {
    memset(this, 0, sizeof(TimerEntry));
  }
  inline void init(double t = TIMER_SLOT_WIDTH) {
    tmin = 1e20;
#if defined(TIMER_HISTOGRAM)
    slotwidth = t;
    inv_slotwidth = (slotwidth > 0 ? 1./slotwidth : 0);
#endif
  }
  /// Start timer, set tcur to seconds since program started.
  inline void start() {
#if defined(DEBUG_TIMER_COLLECTION)
    tcharm = CkWallTimer();
#endif
    tmicro.start();
  }
  /// Stop timer, set tcur to time difference between now and start.
  inline void stop() {
    tcur = tmicro.stop();
#if defined(DEBUG_TIMER_COLLECTION)
    tcharm = CkWallTimer() - tcharm;  // find ellapsed time
    tcharm *= 1e6;  // convert to microseconds
#endif
  }
  /// Update by including tcur into running statistics.
  /// tavg is running average
  /// tvar is unscaled variance
  /// tmin is minimum value so far
  /// tmax is maximum value so far
  inline void update() {
    count++;
    tsum += tcur;
    double delta = tcur - tavg;
    tavg = tavg + delta / count;
    double delta2 = tcur - tavg;
    tvar += delta * delta2;
    if (tcur > tmax) tmax = tcur;
    if (tcur < tmin) tmin = tcur;
#if defined(TIMER_HISTOGRAM)
    int index = int(floor(tcur * inv_slotwidth));
    if (index >= TIMER_SLOTS) index = TIMER_SLOTS - 1;
    hist[index]++;
#endif
  }
  /// Finalize the statistics.
  /// tvar becomes scaled by the count
  /// tstd calculated from tvar
  inline void finalize() {
    if (count > 0) tvar /= count;
    tstd = sqrt(tvar);
    if (tmin > tmax) tmin = tmax;
  }
};

struct TimerSet {
  enum {
    KICK,
    MAXMOVE,
    DRIFT,
    PISTON,
    SUBMITHALF,
    VELBBK1,
    VELBBK2,
    RATTLE1,
    SUBMITFULL,
    SUBMITCOLLECT,
    NUMTIMERS
  };
  TimerEntry t[NUMTIMERS];
  static const char *tlabel[NUMTIMERS];
};

#define TIMER_INIT(T,TYPE) \
  do { \
    (T).t[TimerSet::TYPE].init(); \
  } while(0)

#define TIMER_INIT_WIDTH(T,TYPE,WIDTH) \
  do { \
    (T).t[TimerSet::TYPE].init(WIDTH); \
  } while(0)

#define TIMER_START(T,TYPE) \
  do { \
    (T).t[TimerSet::TYPE].start(); \
  } while(0)

#if defined(DEBUG_TIMER_COLLECTION)

// For debugging, compare clock_gettime with CkWallTimer.
// The concern is regarding these routines that are on average less than
// 10 us (microseconds) but are sometimes an order of magnitude slower.
//
// Note:  After testing, everything with use of clock_gettime seems
// to be working correctly.
//
#define TIMER_STOP(T,TYPE)  \
  do { \
    (T).t[TimerSet::TYPE].stop(); \
    (T).t[TimerSet::TYPE].update(); \
    double tcur = (T).t[TimerSet::TYPE].tcur; \
    int count = (T).t[TimerSet::TYPE].count; \
    if (tcur >= 100 && patch->patchID == SPECIAL_PATCH_ID) { \
      printf("*** %s timing: %g   count: %d  line: %d  charm: %g\n", \
          (T).tlabel[TimerSet::TYPE], tcur, count, __LINE__, \
          (T).t[TimerSet::TYPE].tcharm); \
    } \
  } while(0)

#else   // no DEBUG

#define TIMER_STOP(T,TYPE)  \
  do { \
    (T).t[TimerSet::TYPE].stop(); \
    (T).t[TimerSet::TYPE].update(); \
  } while(0)

#endif  // DEBUG_TIMER_COLLECTION

#define TIMER_DONE(T) \
  do { \
    for (int i=0;  i < TimerSet::NUMTIMERS;  i++) { \
      (T).t[i].finalize(); \
    } \
  } while(0)

#if defined(TIMER_HISTOGRAM)

#define TIMER_REPORT(T) \
  do { \
    printf("%13s %11s %11s %8s %8s %11s %8s\n", \
        "name", "avg", "std", "min", "max", "sum", "calls"); \
    printf("---------------------------------------------------------------" \
        "-------------\n"); \
    for (int i=0;  i < TimerSet::NUMTIMERS;  i++) { \
      printf("%13s %11g %11g %8g %8g %11g %8d\n", \
          (T).tlabel[i], (T).t[i].tavg, (T).t[i].tstd, \
          (T).t[i].tmin, (T).t[i].tmax, (T).t[i].tsum, (T).t[i].count); \
    } \
    printf("---------------------------------------------------------------" \
        "-------------\n"); \
    for (int i=0;  i < TimerSet::NUMTIMERS;  i++) { \
      printf("%13s %8s %8s %8s\n", \
          (T).tlabel[i], "slot", "time", "count"); \
      for (int j=0;  j < TIMER_SLOTS;  j++) { \
        printf("%13s %8d %8g %8d\n", \
            " ", j, (j+1)*(T).t[i].slotwidth, (T).t[i].hist[j]); \
      } \
      printf("---------------------------------------------------------------" \
          "-------------\n"); \
    } \
  } while(0)

#else  // no HISTOGRAM

#define TIMER_REPORT(T) \
  do { \
    printf("%13s %11s %11s %8s %8s %11s %8s\n", \
        "name", "avg", "std", "min", "max", "sum", "calls"); \
    printf("---------------------------------------------------------------" \
        "-------------\n"); \
    for (int i=0;  i < TimerSet::NUMTIMERS;  i++) { \
      printf("%13s %11g %11g %8g %8g %11g %8d\n", \
          (T).tlabel[i], (T).t[i].tavg, (T).t[i].tstd, \
          (T).t[i].tmin, (T).t[i].tmax, (T).t[i].tsum, (T).t[i].count); \
    } \
    printf("---------------------------------------------------------------" \
        "-------------\n"); \
  } while(0)

#endif  // TIMER_HISTOGRAM

#else   // no TIMER

#define TIMER_INIT(T,TYPE)   do { } while(0)
#define TIMER_INIT_WIDTH(T,TYPE,WIDTH)  do{ } while(0)
#define TIMER_START(T,TYPE)  do { } while(0)
#define TIMER_STOP(T,TYPE)   do { } while(0)
#define TIMER_DONE(T)        do { } while(0)
#define TIMER_REPORT(T)      do { } while(0)

#endif  // TIMER_COLLECTION


class HomePatch : public Patch {
  friend class PatchMgr;
  friend class Sequencer;
  friend class ComputeGlobal;

private: 
  // for PatchMgr to use only
  HomePatch(PatchID, FullAtomList&);

  void reinitAtoms(FullAtomList&);
  ScaledPosition min, max, center;
  BigReal aAwayDist, bAwayDist, cAwayDist;

  Bool doAtomUpdate;  // atom changes other than migration

  //Note: If new proxies are added to this HomePatch
  // after load balancing, and it is not the immediate step
  // after atom migration (where ProxyAllMsg will be sent), 
  // then the CompAtomExt list has to be resent with the 
  // ProxyDataMsg (the normal proxy msg when atoms don't 
  // migrate), otherwise, program will crash without such 
  // information when doing force calculations --Chao Mei
  Bool isNewProxyAdded;
  int numGBISP1Arrived, numGBISP2Arrived, numGBISP3Arrived;
  bool phase1BoxClosedCalled;
  bool phase2BoxClosedCalled;
  bool phase3BoxClosedCalled;

public:
  ~HomePatch();

  // Message from ProxyPatch (via ProxyMgr) which registers its existence
  void registerProxy(RegisterProxyMsg *);
  // opposite of above
  void unregisterProxy(UnregisterProxyMsg *);

  // ProxyPatch sends Forces back to here (via ProxyMgr)  
  void receiveResults(ProxyResultVarsizeMsg *msg);
  void receiveResults(ProxyResultMsg *msg);     
  //gbis receiving results from intermediate phases
  void receiveResult(ProxyGBISP1ResultMsg *msg);//after P1
  void receiveResult(ProxyGBISP2ResultMsg *msg);//after P2
  
  //direct function calls, not as entry methods
  void receiveResults(ProxyCombinedResultRawMsg *msg);

  // AtomMigration messages passes from neighbor HomePatches to here.
  void depositMigration(MigrateAtomsMsg *);

  // Bind a Sequencer to this HomePatch
  void useSequencer(Sequencer *sequencerPtr);
  // start simulation over this Patch of atoms
  void runSequencer(void);
  
  //--------------------------------------------------------------------
  // methods for Sequencer to use
  //

  // Signal HomePatch that positions stored are to be now to be used
  void positionsReady(int doMigration=0);
  int marginViolations;

  // methods to implement integration
  void saveForce(const int ftag = Results::normal);
  void addForceToMomentum(
      FullAtom       * __restrict atom_arr,
      const Force    * __restrict force_arr,
      const BigReal    dt,
      int              num_atoms
      )
#if !defined(WIN32) && !defined(WIN64)
    __attribute__((__noinline__))
#endif
    ;
  void addForceToMomentum3(
      FullAtom       * __restrict atom_arr,
      const Force    * __restrict force_arr1,
      const Force    * __restrict force_arr2,
      const Force    * __restrict force_arr3,
      const BigReal    dt1,
      const BigReal    dt2,
      const BigReal    dt3,
      int              num_atoms
      ) 
#if !defined(WIN32) && !defined(WIN64)
    __attribute__((__noinline__))
#endif
    ;
  void addVelocityToPosition(
      FullAtom       * __restrict atom_arr,
      const BigReal    dt,
      int              num_atoms
      ) 
#if !defined(WIN32) && !defined(WIN64)
    __attribute__((__noinline__))
#endif
    ;

  // impose hard wall constraint on Drude bond length
  int hardWallDrude(const BigReal, Tensor *virial, SubmitReduction *);

  // methods for rigidBonds
  struct RattleList {
    int ig;
    int icnt;
  };

  std::vector<int> settleList;
  std::vector<RattleList> rattleList;
  std::vector<RattleParam> rattleParam;
  std::vector<int> noconstList;

  bool rattleListValid;

  // Array to store new positions and velocities. Allocated in "buildRattleList" to size numAtoms
  std::vector<Vector> velNew;
  std::vector<Vector> posNew;

  void addRattleForce(const BigReal invdt, Tensor& wc);

  void buildRattleList();
  int rattle1old(const BigReal, Tensor *virial, SubmitReduction *);
  int rattle1(const BigReal, Tensor *virial, SubmitReduction *);
  void rattle2(const BigReal, Tensor *virial);
  void minimize_rattle2(const BigReal, Tensor *virial, bool forces=false);

  // methods for mollified impluse (MOLLY)
  void mollyAverage();
  void mollyMollify(Tensor *virial);
//  Bool average(Vector qtilde[],const Vector q[],BigReal lambda[],const int n,const int m, const BigReal imass[], const BigReal length2[], const int ial[], const int ilb[], const Vector qji[], const BigReal tolf, const int ntrial);
//  void mollify(Vector qtilde[],const Vector q0[],const BigReal lambda[], Vector force[],const int n, const int m, const BigReal imass[],const int ial[],const int ibl[],const Vector refab[]); 
  
  // BEGIN LA
  void loweAndersenVelocities();
  void loweAndersenFinish();
  // END LA

  void setGBISIntrinsicRadii();
  void gbisComputeAfterP1();//calculate bornRad
  void gbisComputeAfterP2();//calculate dHdrPrefix or self energies
  void gbisP2Ready();
  void gbisP3Ready();

  //LCPO
  void setLcpoType();

  // methods for CONTRA, etc
  void checkpoint(void);
  void revert(void);

  void exchangeCheckpoint(int scriptTask, int &bpc);
  void recvCheckpointReq(int task, const char *key, int replica, int pe);
  void recvCheckpointLoad(CheckpointAtomsMsg *msg);
  void recvCheckpointStore(CheckpointAtomsMsg *msg);
  void recvCheckpointAck();
  int checkpoint_task;
  struct checkpoint_t {
    Lattice lattice;
    int berendsenPressure_count;
    int numAtoms;
    ResizeArray<FullAtom> atoms;
  };
  std::map<std::string,checkpoint_t*> checkpoints;

  // replica exchange
  void exchangeAtoms(int scriptTask);
  void recvExchangeReq(int req);
  void recvExchangeMsg(ExchangeAtomsMsg *msg);
  int exchange_dst;
  int exchange_src;
  int exchange_req;
  ExchangeAtomsMsg *exchange_msg;

  // methods for QM (ExtForces replacement)
  void replaceForces(ExtForce *f);

  void qmSwapAtoms();
  
  // load-balancing trigger
  void submitLoadStats(int timestep);

  // for ComputeHomePatches
  FullAtomList &getAtomList() { return (atom); }

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  // build spanning tree for proxy nodes
  void buildNodeAwareSpanningTree(void);
  void setupChildrenFromProxySpanningTree();
#else
    // build spanning tree for proxy nodes
  void buildSpanningTree(void);
#endif

  void sendNodeAwareSpanningTree();
  void recvNodeAwareSpanningTree(ProxyNodeAwareSpanningTreeMsg *msg);

  void sendSpanningTree();
  void recvSpanningTree(int *t, int n);


  void sendProxies();

#if USE_TOPOMAP 
  int findSubroots(int dim, int* subroots, int psize, int* pidscopy);
#endif

  LDObjHandle ldObjHandle;

#ifdef TIMER_COLLECTION
  TimerSet timerSet;
#endif
protected:
  virtual void boxClosed(int);

  // Internal Atom Migration methods and data
  void doPairlistCheck();
  void doGroupSizeCheck();
  void doMarginCheck();
  void doAtomMigration();
  int inMigration;
  int numMlBuf;
  MigrateAtomsMsg *msgbuf[PatchMap::MaxOneAway];
  
private:
  // Store of Atom-wise variables
  FullAtomList  atom;
  ForceList f_saved[Results::maxNumForces];
  ExtForce *replacementForces;

  CudaAtomList cudaAtomList;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    FullAtomList tempAtom;  // A temporary array used to sort waters
                            //   from non-waters in the atom array
    void separateAtoms();   // Function to separate the atoms currently in atoms.
    void mergeAtomList(FullAtomList &al);  // Function to combine and separate
                                           //   the atoms in al with atoms.
  #endif


  // checkpointed state
  FullAtomList  checkpoint_atom;
  Lattice  checkpoint_lattice;

  // DMK - Atom Separation (water vs. non-water)
  #if NAMD_SeparateWaters != 0
    int checkpoint_numWaterAtoms;
  #endif


  // checkPairlist data
  CompAtomList doPairlistCheck_positions;
  Lattice doPairlistCheck_lattice;
  BigReal doPairlistCheck_newTolerance;

  // MOLLY data
  ResizeArray<BigReal> molly_lambda;
  
  // List of Proxies
  NodeIDList proxy;
  
  Sequencer  *sequencer;

  // Needed for initialization
  int patchMapRead;
  void readPatchMap();

  // Atom Migration internals
  int allMigrationIn;
  int migrationSuspended;
  int patchMigrationCounter;
  int numNeighbors;
  MigrationInfo realInfo[PatchMap::MaxOneAway];
  MigrationInfo *mInfo[3][3][3];

#ifdef NODEAWARE_PROXY_SPANNINGTREE
  //the whole spanning tree for all the proxies this home patch has
  proxyTreeNodeList ptnTree;
  //the immediate children (recording pe ids) containing two parts: 
  //one part of them all belong to the physical node this home patch
  // resides on; the other part of pes belong to all external nodes.
  /* Moved to Patch.h */ 
  //int *children;
  //int numChild;
#else
  NodeIDList tree;              // the whole tree
  int *child;	// spanning tree of proxies - immediate children
  int nChild;
#endif

  // Cached settle1 parameters
  int settle_initialized;
  BigReal settle_mOrmT; BigReal settle_mHrmT; BigReal settle_ra;
  BigReal settle_rb; BigReal settle_rc; BigReal settle_rra;

  /**
   * Redistribute all lonepair forces (of any kind). This may include a direct
   * correction to the virial.
   */
  void redistrib_lonepair_forces(const int, Tensor *);
  
  // single topology force redistribution
  void redistrib_alchpair_forces(const int);

  // PLF -- for TIP4P
  //void redistrib_tip4p_force(Vector&, Vector&, Vector&, Vector&, int, Tensor*);
  void redistrib_tip4p_forces(const int, Tensor*);
  void tip4_omrepos(Vector*, Vector*, Vector*, BigReal);
  void init_tip4();

  // Drude SWM4
  void redistrib_swm4_forces(const int, Tensor*);
  void swm4_omrepos(Vector*, Vector*, Vector*, BigReal);
  void init_swm4();

  /**
   * Reposition lonepair i in a colinear fashion relative to its hosts j and k
   * and according to a fixed distance and scaled vector magnitude between the
   * two hosts.
   */
  void reposition_colinear_lonepair(
      Vector& ri, const Vector& rj, const Vector& rk, Real distance,
      Real scale);

  /**
   * Reposition a lonepair i relative to its hosts j, k, and l according to a
   * given distance, angle, and dihedral formed with the three hosts.
   */
  void reposition_relative_lonepair(
      Vector& ri, const Vector& rj, const Vector& rk, const Vector& rl,
      Real distance, Real angle, Real dihedral);

  /**
   * Reposition all lonepairs (of any kind).
   */
  void reposition_all_lonepairs(void);

  //single topology end state reposition
  void reposition_alchpair(Vector& ri, Vector& rj, Mass& Mi, Mass& Mj);
  void reposition_all_alchpairs(void); //single topolofy alch 

  /**
   * Redistribute the force on a colinear lonepair onto its hosts.
   */
  void redistrib_colinear_lp_force(
       Vector& fi, Vector& fj, Vector& fk,
       const Vector& ri, const Vector& rj, const Vector& rk,
       Real distance, Real scale);

  /**
   * Redistribute the force on a relative lonepair onto its hosts.
   */
  void redistrib_relative_lp_force(
      Vector& fi, Vector& fj, Vector& fk, Vector& fl,
      const Vector& ri, const Vector& rj, const Vector& rk, const Vector& rl,
      Tensor *virial, int midpt);

  //single topology force transfer  
  void redistrib_ap_force(Vector& fi, Vector& fj);
  /**
   * Redistribute the force on a water (TIP4P, SWM4) lonepair onto its hosts.
   * This is similar to redistrib_relative_lp_force but specialized for the
   * bisector case.
   */
  void redistrib_lp_water_force(
      Vector& f_ox, Vector& f_h1, Vector& f_h2, Vector& f_lp,
      const Vector& p_ox, const Vector& p_h1, const Vector& p_h2,
      const Vector& p_lp, Tensor *virial);

  BigReal r_om, r_ohc;
  void write_tip4_props(void);

  int isProxyChanged;

#if CMK_PERSISTENT_COMM
  PersistentHandle *localphs;
  int nphs;
#endif
};

#endif

