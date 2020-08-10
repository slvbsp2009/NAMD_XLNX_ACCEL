/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef NAMDTYPES_H

#define NAMDTYPES_H

#include "common.h"
#include "Vector.h"
#ifndef __CUDACC__
#include "ResizeArray.h"
#endif

class Patch;
class Compute;

typedef Vector Position;
typedef Vector Velocity;

//#ifdef ARCH_POWERPC
//typedef AlignVector Force;
//#else
typedef Vector Force;
//#endif

typedef int AtomID;
typedef int AtomType;
typedef float Mass;
typedef float Charge;

typedef double Coordinate;

struct Transform
{
  signed char i,j,k;
  Transform(void) { i=0; j=0; k=0; }
};

/*
 * 1. "position" field in this structure is very important since it
 * needs to be sent to every patch after every timestep.
 * 2. Anything that is static (value is decided before computation)
 * or only changes after atom migration should be put into the CompAtomExt structure
 * 3. This data structure is 32-byte long which is particularly optimized for some machines
 * (including BG/L) for better cache and message performance. Therefore, changes
 * to this structure should be cautioned for the sake of performance.
 */

struct CompAtom {
  Position position;
  Charge charge;
  short vdwType;
  unsigned char partition;
  unsigned int nonbondedGroupSize : 3;
  unsigned int hydrogenGroupSize : 4;  // could be 3 if unsigned
  unsigned int isWater : 1;  // 0 = particle is not in water, 1 = is in water
};

#ifdef NAMD_KNL
struct CompAtomFlt {
  FloatVector position;
  int32 vdwType;
};
#endif

//CompAtomExt is now needed even in normal case
//for changing the packed msg type related to
//ProxyPatch into varsize msg type where
// two types of proxy msgs (originally, the msg 
// for the step where atoms migrate (ProxyAllMsg), 
// and  the msg for normal steps (ProxyDataMsg))
// are declared as the same class (ProxyDataMsg).
// Note that in normal case, the class is needed
// just for passing the compilation, but not involved
// in the actual force calculation.
// --Chao Mei

typedef int SigIndex;
typedef int AtomSigID;
typedef int ExclSigID;

struct CompAtomExt {
  #if defined(NAMD_CUDA) || defined(NAMD_MIC)
  int sortOrder;  // used to reorder atoms for CUDA
  #endif
  #ifdef MEM_OPT_VERSION
  int id;
  ExclSigID exclId;
  int sigId : 30;  // AtomSigID sigId;
  #else
  int id : 30;  // minimum for 100M atoms is 28 signed, 27 unsigned
  #endif
  unsigned int atomFixed : 1;
  unsigned int groupFixed : 1;
};

struct FullAtom : CompAtom, CompAtomExt{
  Velocity velocity;
  Position fixedPosition;
  double recipMass;
  /**< The reciprocal mass is set to 1/mass or to 0 for massless particles.
   * Calculating this apriori allows us to remove the divide instruction 
   * from the integration loops and the Langevin velocity updates. 
   */
  Mass mass;
  union{
      Real langevinParam;
#ifdef MEM_OPT_VERSION
      int hydVal;
#endif      
  };  
  int32 status;
  Transform transform;
  int migrationGroupSize;
  Real rigidBondLength;

#ifdef MEM_OPT_VERSION
  int outputRank;
#endif

#ifdef MEM_OPT_VERSION
  //a HACK to re-sort FullAtom list used in Parallel IO
  //When every home patch processor receives its atoms list for a patch,
  //the atoms inside this patch may not sorted according to hydList value
  //To save space, use anonymous union data structure to share the space
  //of "langevinParam" to store "hydList" from an InputAtom and then sort the 
  //atom list. The "langevinParam" value is not initialized until home 
  //patch creation -Chao Mei
  int operator < (const FullAtom &a) const {
      return hydVal < a.hydVal;
  }
#endif
};

//InputAtom is used to contain the info of the atoms
//loaded into input processors.
struct InputAtom: FullAtom{
	bool isValid;
	short isGP;
	short isMP;
	int hydList;
	int GPID;
	int MPID;
    	
	int operator < (const InputAtom &a) const{
		return hydList < a.hydList;
	}
};

struct CudaAtom {
  float x,y,z,q;
};

struct CudaForce {
  float x, y, z;
};

#ifndef __CUDACC__
typedef ResizeArray<CudaAtom> CudaAtomList;
typedef ResizeArray<CompAtom> CompAtomList;
typedef ResizeArray<CompAtomExt> CompAtomExtList;
#ifdef NAMD_KNL
typedef ResizeArray<CompAtomFlt> CompAtomFltList;
#endif
typedef ResizeArray<FullAtom> FullAtomList;
typedef ResizeArray<InputAtom> InputAtomList;
typedef ResizeArray<Position> PositionList;
typedef ResizeArray<Velocity> VelocityList;
typedef ResizeArray<Force> ForceList;
typedef ResizeArray<Transform> TransformList;

typedef ResizeArray<AtomID> AtomIDList;
typedef ResizeArray<BigReal> BigRealList;
typedef ResizeArray<Real> RealList;
typedef float GBReal;
typedef ResizeArray<GBReal> GBRealList;
typedef ResizeArray<int> IntList;

typedef int PatchID;
typedef int ComputeID;
typedef int NodeID;

typedef ResizeArray<PatchID> PatchIDList;
typedef ResizeArray<Patch *> PatchList;

typedef ResizeArray<Compute *> ComputeList;

// See AtomMap
struct LocalID
{
  PatchID pid;
  int index;
};

typedef ResizeArray<NodeID> NodeIDList;

struct ExtForce {
  int replace;
  Force force;
  ExtForce() : replace(0) {;}
};


// DMK - Atom Sort
#if NAMD_ComputeNonbonded_SortAtoms != 0

  typedef struct __sort_entry {
    int index;  // Index of atom in CompAtom array
    BigReal sortValue;   // Distance of PAp from P0 (see calculation code)
  } SortEntry;

#endif

//This class represents a tree node of proxy spanning tree
//All pes in this array have the same "nodeID". In other words,
//all those pes are in the same physical node.
//This is a structure for adapting NAMD to multicore processors
struct proxyTreeNode{
    int nodeID;
    int *peIDs;
    int numPes;

    proxyTreeNode(){
        nodeID = -1;
        peIDs = NULL;
        numPes = 0;
    }
    proxyTreeNode(int nid, int numPes_, int *pes){
        nodeID = nid;
        numPes = numPes_;
        peIDs = new int[numPes];
        memcpy(peIDs, pes, sizeof(int)*numPes);
    }

    inline proxyTreeNode(const proxyTreeNode &n){
        nodeID = n.nodeID;
        numPes = n.numPes;
        if(numPes==0) {
            peIDs = NULL;
        }else{
            peIDs = new int[n.numPes];
            memcpy(peIDs, n.peIDs, sizeof(int)*numPes);
        }
    }
    inline proxyTreeNode &operator=(const proxyTreeNode &n){
        nodeID = n.nodeID;
        numPes = n.numPes;
        delete [] peIDs;
        if(numPes==0) {
            peIDs = NULL;
            return (*this);
        }
        peIDs = new int[n.numPes];
        memcpy(peIDs, n.peIDs, sizeof(int)*numPes);
        return (*this);
    }
    ~proxyTreeNode(){
        delete [] peIDs;
    }
};

typedef ResizeArray<proxyTreeNode> proxyTreeNodeList;
#endif // __CUDACC__


//
// When defined, use NVTX to record CPU activity into CUDA profiling run.
//
#undef PUSH_RANGE
#undef POP_RANGE
#undef RANGE

#if defined(NAMD_CUDA) && defined(NAMD_USE_NVTX)

#include <nvToolsExt.h>

// C++ note: declaring const variables implies static (internal) linkage,
// and you have to explicitly specify "extern" to get external linkage.
const uint32_t NAMD_nvtx_colors[] = {
  0x0000ff00,
  0x000000ff,
  0x00ffff00,
  0x00ff00ff,
  0x0000ffff,
  0x00ff0000,
  0x00ffffff,
};
const int NAMD_nvtx_colors_len = sizeof(NAMD_nvtx_colors)/sizeof(uint32_t);

// start recording an event
#define PUSH_RANGE(name,cid) \
  do { \
    int color_id = cid; \
    color_id = color_id % NAMD_nvtx_colors_len; \
    nvtxEventAttributes_t eventAttrib = {0}; \
    eventAttrib.version = NVTX_VERSION; \
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
    eventAttrib.colorType = NVTX_COLOR_ARGB; \
    eventAttrib.color = NAMD_nvtx_colors[color_id]; \
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
    eventAttrib.message.ascii = name; \
    nvtxRangePushEx(&eventAttrib); \
  } while(0)  // must terminate with semi-colon

// stop recording an event
#define POP_RANGE \
  nvtxRangePop()
  // must terminate with semi-colon

// embed event recording in class to automatically pop when destroyed
class NAMD_NVTX_Tracer {
  public:
    NAMD_NVTX_Tracer(const char *name, int cid = 0) { PUSH_RANGE(name, cid); }
    ~NAMD_NVTX_Tracer() { POP_RANGE; }
};

// include cid as part of the name
// call RANGE at beginning of function to push event recording
// destructor is automatically called on return to pop event recording
#define RANGE(name,cid) \
  NAMD_NVTX_Tracer namd_nvtx_tracer##cid(name,cid)
  // must terminate with semi-colon

#else

//
// Otherwise the NVTX profiling macros become no-ops.
//
#define PUSH_RANGE(name,cid) do { } while(0)  // must terminate with semi-colon
#define POP_RANGE            do { } while(0)  // must terminate with semi-colon
#define RANGE(namd,cid)      do { } while(0)  // must terminate with semi-colon

#endif // NAMD_CUDA && NAMD_USE_NVTX


#endif /* NAMDTYPES_H */

