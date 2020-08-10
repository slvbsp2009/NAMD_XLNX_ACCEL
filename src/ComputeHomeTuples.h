/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTEHOMETUPLES_H
#define COMPUTEHOMETUPLES_H

#ifdef USE_HOMETUPLES
#include <vector>
#endif
#include "NamdTypes.h"
#include "common.h"
#include "structures.h"
#include "Compute.h"
#include "HomePatch.h"

#include "Box.h"
#include "OwnerBox.h"
#include "UniqueSet.h"

#include "Node.h"
#include "SimParameters.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "PatchMgr.h"
#include "ProxyMgr.h"
#include "HomePatchList.h"
#include "Molecule.h"
#include "Parameters.h"
#include "ReductionMgr.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "Priorities.h"
#include "LdbCoordinator.h"

class TuplePatchElem {
  public:
    PatchID patchID;
    Patch *p;
    Box<Patch,CompAtom> *positionBox;
    Box<Patch,CompAtom> *avgPositionBox;
    Box<Patch,Results> *forceBox;
    CompAtom *x;
    CompAtomExt *xExt;
    CompAtom *x_avg;
    Results *r;
    Force *f;
    Force *af;

    int hash() const { return patchID; }

  TuplePatchElem(PatchID pid = -1) {
    patchID = pid;
    p = NULL;
    positionBox = NULL;
    avgPositionBox = NULL;
    forceBox = NULL;
    x = NULL;
    xExt = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
    af = NULL;
  }

  TuplePatchElem(Patch *p_param, Compute *cid) {
    patchID = p_param->getPatchID();
    p = p_param;
    positionBox = p_param->registerPositionPickup(cid);
    avgPositionBox = p_param->registerAvgPositionPickup(cid);
    forceBox = p_param->registerForceDeposit(cid);
    x = NULL;
    xExt = NULL;
    x_avg = NULL;
    r = NULL;
    f = NULL;
    af = NULL;
  }
    
  ~TuplePatchElem() {};

  int operator==(const TuplePatchElem &elem) const {
    return (elem.patchID == patchID);
  }

  int operator<(const TuplePatchElem &elem) const {
    return (patchID < elem.patchID);
  }
};

typedef UniqueSet<TuplePatchElem> TuplePatchList;
typedef UniqueSetIter<TuplePatchElem> TuplePatchListIter;

class AtomMap;
class ReductionMgr;

#ifdef MEM_OPT_VERSION
template <class T> struct ElemTraits {
  typedef AtomSignature signature;
  static signature* get_sig_pointer(Molecule *mol) { return mol->atomSigPool; }
  static int get_sig_id(const CompAtomExt &a) { return a.sigId; }
};

template <> struct ElemTraits <ExclElem> {
  typedef ExclusionSignature signature;
  static signature* get_sig_pointer(Molecule *mol) { return mol->exclSigPool; }
  static int get_sig_id(const CompAtomExt &a) { return a.exclId; }
};
#endif

#ifdef USE_HOMETUPLES
//
// Simple base class for HomeTuples and SelfTuples that stores the type of the tuple
//
class Tuples {
private:
  int type;
protected:
  Tuples(int type) : type(type) {}
public:
  // Tuple types
  enum {BOND=0, ANGLE, DIHEDRAL, IMPROPER, EXCLUSION, CROSSTERM, NUM_TUPLE_TYPES};

  int getType() {return type;}
  virtual void submitTupleCount(SubmitReduction *reduction, int tupleCount)=0;
  // virtual void copyTupleData(void* tupleData)=0;
  virtual int getNumTuples()=0;
  virtual void* getTupleList()=0;
  virtual void loadTuples(TuplePatchList& tuplePatchList, const char* isBasePatch, AtomMap *atomMap,
      const std::vector<int>& pids = std::vector<int>())=0;
};

//
// HomeTuples class. These are created and stored in ComputeBondedCUDA::registerCompute()
// e.g.: new HomeTuples<BondElem, Bond, BondValue>(BOND)
//
template <class T, class S, class P> class HomeTuples : public Tuples {
  protected:
    std::vector<T> tupleList;

  public:

    HomeTuples(int type=-1) : Tuples(type) {}

#if __cplusplus < 201103L
#define final
#endif

    virtual void* getTupleList() final {
      return (void*)tupleList.data();
    }

    virtual void submitTupleCount(SubmitReduction *reduction, int tupleCount) final {
      reduction->item(T::reductionChecksumLabel) += (BigReal)tupleCount;
    }

    // virtual void copyTupleData(void* tupleData) final {
      // for (int i=0;i < tupleList.size();i++) {
      //   tupleData[i] = 
      // }
      // T::loadTupleData(tupleData);
    // }

    virtual int getNumTuples() final {
      return tupleList.size();
    }

    virtual void loadTuples(TuplePatchList& tuplePatchList, const char* isBasePatch, AtomMap *atomMap,
      const std::vector<int>& pids = std::vector<int>()) {

      if (isBasePatch == NULL) {
        NAMD_bug("NULL isBasePatch detected in HomeTuples::loadTuples()");
      }

      int numTuples;

#ifdef MEM_OPT_VERSION
      typename ElemTraits<T>::signature *allSigs;      
#else
      int32 **tuplesByAtom;
      /* const (need to propagate const) */ S *tupleStructs;
#endif
      
      const P *tupleValues;
      Node *node = Node::Object();
      PatchMap *patchMap = PatchMap::Object();
      // AtomMap *atomMap = AtomMap::Object();

#ifdef MEM_OPT_VERSION
      allSigs = ElemTraits<T>::get_sig_pointer(node->molecule);
#else      
      T::getMoleculePointers(node->molecule,
        &numTuples, &tuplesByAtom, &tupleStructs);      
#endif
      
      T::getParameterPointers(node->parameters, &tupleValues);

      tupleList.clear();

      LocalID aid[T::size];
      int partition[T::size];

      const int lesOn = node->simParameters->lesOn;
      const int soluteScalingOn = node->simParameters->soluteScalingOn;
      const int fepOn = node->simParameters->singleTopology;
      const int sdScaling = node->simParameters->sdScaling;
      Real invLesFactor = lesOn ? 
                          1.0/node->simParameters->lesFactor :
                          1.0;
      const Real soluteScalingFactor = node->simParameters->soluteScalingFactor;
      const Bool soluteScalingAll = node->simParameters->soluteScalingAll;
      BigReal OneMinusLambda = 1.0 - node->simParameters->alchLambda;
      BigReal Lambda = node->simParameters->alchLambda;
      const int num_unpert_bonds = node->molecule->num_alch_unpert_Bonds;
      const int num_unpert_angles = node->molecule->num_alch_unpert_Angles;
      const int num_unpert_dihedrals = node->molecule->num_alch_unpert_Dihedrals;
      Bond *unpert_bonds = node->molecule->alch_unpert_bonds;
      Angle *unpert_angles = node->molecule->alch_unpert_angles;
      Dihedral *unpert_dihedrals = node->molecule->alch_unpert_dihedrals;

      // cycle through each patch and gather all tuples
      TuplePatchListIter ai(tuplePatchList);
      if (pids.size() == 0) ai = ai.begin();

      int numPid = (pids.size() == 0) ? tuplePatchList.size() : pids.size();

      for (int ipid=0;ipid < numPid;ipid++) {
        // Patch *patch;
        int numAtoms;
        CompAtomExt *atomExt;
        // Take next patch
        if (pids.size() == 0) {
          Patch* patch = (*ai).p;
          numAtoms = patch->getNumAtoms();
          atomExt = (*ai).xExt;
          ai++;
        } else {
          TuplePatchElem *tpe = tuplePatchList.find(TuplePatchElem(pids[ipid]));
          Patch* patch = tpe->p;
          numAtoms = patch->getNumAtoms();
          atomExt = tpe->xExt;          
        }
   
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {
          /* cycle through each tuple */
#ifdef MEM_OPT_VERSION
          typename ElemTraits<T>::signature *thisAtomSig =
                   &allSigs[ElemTraits<T>::get_sig_id(atomExt[j])];
          TupleSignature *allTuples;
          T::getTupleInfo(thisAtomSig, &numTuples, &allTuples);
          for(int k=0; k<numTuples; k++) {
            T t(atomExt[j].id, &allTuples[k], tupleValues);
#else
          /* get list of all tuples for the atom */
          int32 *curTuple = tuplesByAtom[atomExt[j].id];
          for( ; *curTuple != -1; ++curTuple) {
            T t(&tupleStructs[*curTuple],tupleValues);
#endif            
            register int i;
            aid[0] = atomMap->localID(t.atomID[0]);
            int homepatch = aid[0].pid;
            int samepatch = 1;
            partition[0] = fepOn ? node->molecule->get_fep_type(t.atomID[0]) : 0;  //using atom partition to determine if a bonded term to be scaled by lambda or 1-lambda in single topology relative FEP. 
            int has_les = lesOn ? node->molecule->get_fep_type(t.atomID[0]) : 0;
            int has_ss = soluteScalingOn ? node->molecule->get_ss_type(t.atomID[0]) : 0;
            int is_fep_ss = partition[0] > 2;
            int is_fep_sd = 0;
            int fep_tuple_type = 0;
            for (i=1; i < T::size; i++) {
              aid[i] = atomMap->localID(t.atomID[i]);
              samepatch = samepatch && ( homepatch == aid[i].pid );
              partition[i] = fepOn ? node->molecule->get_fep_type(t.atomID[i]) : 0;
              has_les |= lesOn ? node->molecule->get_fep_type(t.atomID[i]) : 0;
              has_ss |= soluteScalingOn ? node->molecule->get_ss_type(t.atomID[i]) : 0;
              if (fepOn) {
              is_fep_ss &= partition[i] > 2;
              is_fep_sd |= (abs(partition[i] - partition[0]) == 2);
              fep_tuple_type = partition[i]; }
            }
            if (sdScaling && is_fep_sd) {
              // check if this bonded term is one of Shobana term.
              // This segment looks ugly and not GPU friendly,
              // and might not appear in GPU code.
              //
              // XXX Could optimize in a number of ways:
              // - could switch on T::size, then loop for that sized tuple
              // - could use hash table to look up unpert_*[] elements
              // - could add flag field to BondElem, et al., classes
              for (i=0; i < num_unpert_bonds; i++) {
                if (T::size == 2
                    && t.atomID[0]==unpert_bonds[i].atom1
                    && t.atomID[1]==unpert_bonds[i].atom2) is_fep_sd = 0;
              }
              for (i=0; i < num_unpert_angles; i++) {
                if (T::size == 3
                    && t.atomID[0]==unpert_angles[i].atom1
                    && t.atomID[1]==unpert_angles[i].atom2
                    && t.atomID[2]==unpert_angles[i].atom3) is_fep_sd = 0;
              }
              for (i=0; i < num_unpert_dihedrals; i++) {
                if (T::size == 4
                    && t.atomID[0]==unpert_dihedrals[i].atom1
                    && t.atomID[1]==unpert_dihedrals[i].atom2
                    && t.atomID[2]==unpert_dihedrals[i].atom3
                    && t.atomID[3]==unpert_dihedrals[i].atom4) is_fep_sd = 0;
              }
            }
            if (T::size < 4 && !soluteScalingAll) has_ss = false;
            if ( samepatch ) continue;
            t.scale = (!has_les && !has_ss) ? 1.0 : ( has_les ? invLesFactor : soluteScalingFactor );
            if (is_fep_ss) t.scale = (fep_tuple_type == 4) ? OneMinusLambda : Lambda;
            if (is_fep_sd && sdScaling) t.scale = (fep_tuple_type == 4 || fep_tuple_type == 2) ? OneMinusLambda : Lambda;

            for (i=1; i < T::size; i++) {
              homepatch = patchMap->downstream(homepatch,aid[i].pid);
            }
            if ( homepatch != notUsed && isBasePatch[homepatch] ) {
              TuplePatchElem *p;
              for (i=0; i < T::size; i++) {
                t.p[i] = p = tuplePatchList.find(TuplePatchElem(aid[i].pid));
                if ( ! p ) {
#ifdef MEM_OPT_VERSION
                  iout << iWARN << "Tuple with atoms ";
#else
                  iout << iWARN << "Tuple " << *curTuple << " with atoms ";
#endif
                  int erri;
                  for( erri = 0; erri < T::size; erri++ ) {
                    iout << t.atomID[erri] << "(" <<  aid[erri].pid << ") ";
                  }
                  iout << "missing patch " << aid[i].pid << "\n" << endi;
                  break;
                }
                t.localIndex[i] = aid[i].index;
              }
              if ( ! p ) continue;
#ifdef MEM_OPT_VERSION
              //avoid adding Tuples whose atoms are all fixed
              if(node->simParameters->fixedAtomsOn && !node->simParameters->fixedAtomsForces) {
                int allfixed = 1;
                for(i=0; i<T::size; i++){
                  CompAtomExt *one = &(t.p[i]->xExt[aid[i].index]);
                  allfixed = allfixed & one->atomFixed;
                }
                if(!allfixed) tupleList.push_back(t);
              }else{
                tupleList.push_back(t);
              }
#else
              tupleList.push_back(t);
#endif               
            }
          }
        }
      }
    }

};
#endif

template <class T, class S, class P> class ComputeHomeTuples : public Compute {

  protected:
  
#ifndef USE_HOMETUPLES
    virtual void loadTuples(void) {
      int numTuples;

      #ifdef MEM_OPT_VERSION
      typename ElemTraits<T>::signature *allSigs;      
      #else
      int32 **tuplesByAtom;
      /* const (need to propagate const) */ S *tupleStructs;
      #endif
      
      const P *tupleValues;
      Node *node = Node::Object();

      #ifdef MEM_OPT_VERSION
      allSigs = ElemTraits<T>::get_sig_pointer(node->molecule);
      #else      
      T::getMoleculePointers(node->molecule,
		    &numTuples, &tuplesByAtom, &tupleStructs);      
      #endif
      
      T::getParameterPointers(node->parameters, &tupleValues);

      tupleList.resize(0);

      LocalID aid[T::size];
      int partition[T::size];

      const int lesOn = node->simParameters->lesOn;
      const int soluteScalingOn = node->simParameters->soluteScalingOn;
      const int fepOn = node->simParameters->singleTopology;
      const int sdScaling = node->simParameters->sdScaling;
      Real invLesFactor = lesOn ? 
                          1.0/node->simParameters->lesFactor :
                          1.0;
      const Real soluteScalingFactor = node->simParameters->soluteScalingFactor;
      const Bool soluteScalingAll = node->simParameters->soluteScalingAll;
      BigReal OneMinusLambda = 1.0 - node->simParameters->alchLambda;
      BigReal Lambda = node->simParameters->alchLambda;
      const int num_unpert_bonds = node->molecule->num_alch_unpert_Bonds;
      const int num_unpert_angles = node->molecule->num_alch_unpert_Angles;
      const int num_unpert_dihedrals = node->molecule->num_alch_unpert_Dihedrals;
      Bond *unpert_bonds = node->molecule->alch_unpert_bonds;
      Angle *unpert_angles = node->molecule->alch_unpert_angles;
      Dihedral *unpert_dihedrals = node->molecule->alch_unpert_dihedrals;

      // cycle through each patch and gather all tuples
      TuplePatchListIter ai(tuplePatchList);
    
      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
        // CompAtom *atom = (*ai).x;
        Patch *patch = (*ai).p;
        int numAtoms = patch->getNumAtoms();
	CompAtomExt *atomExt = (*ai).xExt; //patch->getCompAtomExtInfo();
    
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {              
           /* cycle through each tuple */
           #ifdef MEM_OPT_VERSION
           typename ElemTraits<T>::signature *thisAtomSig =
                   &allSigs[ElemTraits<T>::get_sig_id(atomExt[j])];
           TupleSignature *allTuples;
           T::getTupleInfo(thisAtomSig, &numTuples, &allTuples);
           for(int k=0; k<numTuples; k++) {
               T t(atomExt[j].id, &allTuples[k], tupleValues);
           #else
           /* get list of all tuples for the atom */
           int32 *curTuple = tuplesByAtom[atomExt[j].id];
           for( ; *curTuple != -1; ++curTuple) {             
             T t(&tupleStructs[*curTuple],tupleValues);
           #endif            
             register int i;
             aid[0] = atomMap->localID(t.atomID[0]);
             int homepatch = aid[0].pid;
             int samepatch = 1;
             partition[0] = fepOn ? node->molecule->get_fep_type(t.atomID[0]) : 0;
             int has_les = lesOn ? node->molecule->get_fep_type(t.atomID[0]) : 0;
             int has_ss = soluteScalingOn ? node->molecule->get_ss_type(t.atomID[0]) : 0;
             int is_fep_ss = partition[0] > 2;
             int is_fep_sd = 0;
             int fep_tuple_type = 0;
             for (i=1; i < T::size; i++) {
               aid[i] = atomMap->localID(t.atomID[i]);
               samepatch = samepatch && ( homepatch == aid[i].pid );
               partition[i] = fepOn ? node->molecule->get_fep_type(t.atomID[i]) : 0;
               has_les |= lesOn ? node->molecule->get_fep_type(t.atomID[i]) : 0;
               has_ss |= soluteScalingOn ? node->molecule->get_ss_type(t.atomID[i]) : 0;  
               if (fepOn) {
                 is_fep_ss &= partition[i] > 2;
                 is_fep_sd |= (abs(partition[i] - partition[0]) == 2);
                 fep_tuple_type = partition[i];
               }
             }
             if (sdScaling && is_fep_sd) {
               for (i=0; i < num_unpert_bonds; i++) {
                 if (T::size == 2
                     && t.atomID[0]==unpert_bonds[i].atom1
                     && t.atomID[1]==unpert_bonds[i].atom2) is_fep_sd = 0;
               }
               for (i=0; i < num_unpert_angles; i++) {
                 if (T::size == 3
                     && t.atomID[0]==unpert_angles[i].atom1
                     && t.atomID[1]==unpert_angles[i].atom2
                     && t.atomID[2]==unpert_angles[i].atom3) is_fep_sd = 0;
               }
               for (i=0; i < num_unpert_dihedrals; i++) {
                 if (T::size == 4
                     && t.atomID[0]==unpert_dihedrals[i].atom1
                     && t.atomID[1]==unpert_dihedrals[i].atom2
                     && t.atomID[2]==unpert_dihedrals[i].atom3
                     && t.atomID[3]==unpert_dihedrals[i].atom4) is_fep_sd = 0;
               }
             }
             if (T::size < 4 && !soluteScalingAll) has_ss = false;
             if ( samepatch ) continue;
             t.scale = (!has_les && !has_ss) ? 1.0 : ( has_les ? invLesFactor : soluteScalingFactor );
             if (is_fep_ss) t.scale = (fep_tuple_type == 4) ? OneMinusLambda : Lambda;
             if (is_fep_sd && sdScaling) t.scale = (fep_tuple_type == 4 || fep_tuple_type == 2) ? OneMinusLambda : Lambda;

             for (i=1; i < T::size; i++) {
	         homepatch = patchMap->downstream(homepatch,aid[i].pid);
             }
             if ( homepatch != notUsed && isBasePatch[homepatch] ) {
      	       TuplePatchElem *p;
               for (i=0; i < T::size; i++) {
      	         t.p[i] = p = tuplePatchList.find(TuplePatchElem(aid[i].pid));
      	         if ( ! p ) {
                     #ifdef MEM_OPT_VERSION
                     iout << iWARN << "Tuple with atoms ";
                     #else
      	           iout << iWARN << "Tuple " << *curTuple << " with atoms ";
                     #endif
      	           int erri;
      	           for( erri = 0; erri < T::size; erri++ ) {
      	             iout << t.atomID[erri] << "(" <<  aid[erri].pid << ") ";
      	           }
      	           iout << "missing patch " << aid[i].pid << "\n" << endi;
      	           break;
      	         }
      	         t.localIndex[i] = aid[i].index;
               }
      	       if ( ! p ) continue;
             #ifdef MEM_OPT_VERSION
               //avoid adding Tuples whose atoms are all fixed
               if(node->simParameters->fixedAtomsOn &&
		  !node->simParameters->fixedAtomsForces) {
                 int allfixed = 1;
                 for(i=0; i<T::size; i++){
                   CompAtomExt *one = &(t.p[i]->xExt[aid[i].index]);
                   allfixed = allfixed & one->atomFixed;
                 }
                 if(!allfixed) tupleList.add(t);
               }else{
                 tupleList.add(t);
               }
             #else
               tupleList.add(t);
             #endif               
             }
           }
        }
      }
    }
#endif

    int doLoadTuples;
  
  protected:
  
#ifdef USE_HOMETUPLES
    HomeTuples<T, S, P>* tuples;
    TuplePatchList tuplePatchList;
#else
    ResizeArray<T> tupleList;
    TuplePatchList tuplePatchList;
#endif

    PatchMap *patchMap;
    AtomMap *atomMap;
    SubmitReduction *reduction;
    int accelMDdoDihe;
    SubmitReduction *pressureProfileReduction;
    BigReal *pressureProfileData;
    int pressureProfileSlabs;
    char *isBasePatch;
  
    ComputeHomeTuples(ComputeID c) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      
      SimParameters *params = Node::Object()->simParameters;
      accelMDdoDihe=false;
      if (params->accelMDOn) {
         if (params->accelMDdihe || params->accelMDdual) accelMDdoDihe=true;
      }
      if (params->pressureProfileOn) {
        pressureProfileSlabs = T::pressureProfileSlabs = 
          params->pressureProfileSlabs;
        int n = T::pressureProfileAtomTypes = params->pressureProfileAtomTypes;
        pressureProfileReduction = ReductionMgr::Object()->willSubmit(
          REDUCTIONS_PPROF_BONDED, 3*pressureProfileSlabs*((n*(n+1))/2));
        int numAtomTypePairs = n*n;
        pressureProfileData = new BigReal[3*pressureProfileSlabs*numAtomTypePairs];
      } else {
        pressureProfileReduction = NULL;
        pressureProfileData = NULL;
      }
      doLoadTuples = false;
      isBasePatch = 0;
#ifdef USE_HOMETUPLES
      tuples = NULL;
#endif
    }

    ComputeHomeTuples(ComputeID c, PatchIDList &pids) : Compute(c) {
      patchMap = PatchMap::Object();
      atomMap = AtomMap::Object();
      reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
      SimParameters *params = Node::Object()->simParameters;
      accelMDdoDihe=false;
      if (params->accelMDOn) {
         if (params->accelMDdihe || params->accelMDdual) accelMDdoDihe=true;
      }
      if (params->pressureProfileOn) {
        pressureProfileSlabs = T::pressureProfileSlabs = 
          params->pressureProfileSlabs;
        int n = T::pressureProfileAtomTypes = params->pressureProfileAtomTypes;
        pressureProfileReduction = ReductionMgr::Object()->willSubmit(
          REDUCTIONS_PPROF_BONDED, 3*pressureProfileSlabs*((n*(n+1))/2));
        int numAtomTypePairs = n*n;
        pressureProfileData = new BigReal[3*pressureProfileSlabs*numAtomTypePairs];
      } else {
        pressureProfileReduction = NULL;
        pressureProfileData = NULL;
      }
      doLoadTuples = false;
      int nPatches = patchMap->numPatches();
      isBasePatch = new char[nPatches];
      int i;
      for (i=0; i<nPatches; ++i) { isBasePatch[i] = 0; }
      for (i=0; i<pids.size(); ++i) { isBasePatch[pids[i]] = 1; }
#ifdef USE_HOMETUPLES
      tuples = NULL;
#endif
    }

  public:
  
    virtual ~ComputeHomeTuples() {
      delete reduction;
      delete [] isBasePatch;
      delete pressureProfileReduction;
      delete pressureProfileData;
#ifdef USE_HOMETUPLES
      if (tuples != NULL) delete tuples;
#endif
    }

    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    virtual void initialize(void) {

#ifdef NAMD_CUDA
      ProxyMgr *proxyMgr = ProxyMgr::Object();
#endif

#ifdef USE_HOMETUPLES
      tuples = new HomeTuples<T, S, P>();
#endif

      // Start with empty list
      tuplePatchList.clear();
    
      int nPatches = patchMap->numPatches();
      int pid;
      for (pid=0; pid<nPatches; ++pid) {
        if ( isBasePatch[pid] ) {
#ifdef NAMD_CUDA
          proxyMgr->createProxy(pid);
#endif
          Patch *patch = patchMap->patch(pid);
	  tuplePatchList.add(TuplePatchElem(patch, this));
        }
      }
    
      // Gather all proxy patches (neighbors, that is)
      PatchID neighbors[PatchMap::MaxOneOrTwoAway];
    
      for (pid=0; pid<nPatches; ++pid) if ( isBasePatch[pid] ) {
        int numNeighbors = patchMap->upstreamNeighbors(pid,neighbors);
        for ( int i = 0; i < numNeighbors; ++i ) {
          if ( ! tuplePatchList.find(TuplePatchElem(neighbors[i])) ) {
#ifdef NAMD_CUDA
            proxyMgr->createProxy(neighbors[i]);
#endif
            Patch *patch = patchMap->patch(neighbors[i]);
	    tuplePatchList.add(TuplePatchElem(patch, this));
          }
        }
      }
      setNumPatches(tuplePatchList.size());
      doLoadTuples = true;

      basePriority = COMPUTE_PROXY_PRIORITY;  // no patch dependence
    }

    //======================================================================
    // atomUpdate() - Method is invoked after anytime that atoms have been
    // changed in patches used by this Compute object.
    //======================================================================
    void atomUpdate(void) {
      doLoadTuples = true;
    }

//-------------------------------------------------------------------
// Routine which is called by enqueued work msg.  It wraps
// actualy Force computation with the apparatus needed
// to get access to atom positions, return forces etc.
//-------------------------------------------------------------------
    virtual void doWork(void) {

      LdbCoordinator::Object()->startWork(ldObjHandle);

      // Open Boxes - register that we are using Positions
      // and will be depositing Forces.
      UniqueSetIter<TuplePatchElem> ap(tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->x = ap->positionBox->open();
	ap->xExt = ap->p->getCompAtomExtInfo();
        if ( ap->p->flags.doMolly ) ap->x_avg = ap->avgPositionBox->open();
        ap->r = ap->forceBox->open();
        ap->f = ap->r->f[Results::normal];
        if (accelMDdoDihe) ap->af = ap->r->f[Results::amdf]; // for dihedral-only or dual-boost accelMD
      } 
    
      BigReal reductionData[T::reductionDataSize];
      int tupleCount = 0;
      int numAtomTypes = T::pressureProfileAtomTypes;
      int numAtomTypePairs = numAtomTypes*numAtomTypes;
    
      for ( int i = 0; i < T::reductionDataSize; ++i ) reductionData[i] = 0;
      if (pressureProfileData) {
        memset(pressureProfileData, 0, 3*pressureProfileSlabs*numAtomTypePairs*sizeof(BigReal));
        // Silly variable hiding of the previous iterator
        UniqueSetIter<TuplePatchElem> newap(tuplePatchList);
        newap = newap.begin();
        const Lattice &lattice = newap->p->lattice;
        T::pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
        T::pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;
      }

      if ( ! Node::Object()->simParameters->commOnly ) {
      if ( doLoadTuples ) {
#ifdef USE_HOMETUPLES
        tuples->loadTuples(tuplePatchList, isBasePatch, AtomMap::Object());
#else
        loadTuples();
#endif
        doLoadTuples = false;
      }
      // take triplet and pass with tuple info to force eval
#ifdef USE_HOMETUPLES
      T *al = (T *)tuples->getTupleList();
      const int ntuple = tuples->getNumTuples();
#else
      T *al = tupleList.begin();
      const int ntuple = tupleList.size();
#endif
      if ( ntuple ) T::computeForce(al, ntuple, reductionData, pressureProfileData);
      tupleCount += ntuple;
      }
 
    LdbCoordinator::Object()->endWork(ldObjHandle);

      T::submitReductionData(reductionData,reduction);
      reduction->item(T::reductionChecksumLabel) += (BigReal)tupleCount;
      reduction->submit();

      if (pressureProfileReduction) {
        // For ease of calculation we stored interactions between types
        // i and j in (ni+j).  For efficiency now we coalesce the
        // cross interactions so that just i<=j are stored.
        const int arraysize = 3*pressureProfileSlabs;
        const BigReal *data = pressureProfileData;
        for (int i=0; i<numAtomTypes; i++) {
          for (int j=0; j<numAtomTypes; j++) {
            int ii=i;
            int jj=j;
            if (ii > jj) { int tmp=ii; ii=jj; jj=tmp; }
            const int reductionOffset = 
              (ii*numAtomTypes - (ii*(ii+1))/2 + jj)*arraysize;
            for (int k=0; k<arraysize; k++) {
              pressureProfileReduction->item(reductionOffset+k) += data[k];
            }
            data += arraysize;
          }
        }
        pressureProfileReduction->submit();
      }
    
      // Close boxes - i.e. signal we are done with Positions and
      // AtomProperties and that we are depositing Forces
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->positionBox->close(&(ap->x));
        if ( ap->p->flags.doMolly ) ap->avgPositionBox->close(&(ap->x_avg));
        ap->forceBox->close(&(ap->r));
      }
    }
};


#endif

