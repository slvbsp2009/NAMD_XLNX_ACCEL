/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTESELFTUPLES_H
#define COMPUTESELFTUPLES_H

#include "ComputeHomeTuples.h"
#include "LdbCoordinator.h"

#ifdef USE_HOMETUPLES
template <class T, class S, class P> class SelfTuples : public HomeTuples<T, S, P> {

  public:
    SelfTuples(int type=-1) : HomeTuples<T,S,P>(type) {}

  private:

    virtual void loadTuples(TuplePatchList& tuplePatchList, const char* isBasePatch, AtomMap *atomMap,
      const std::vector<int>& pids = std::vector<int>()) {

      if (isBasePatch != NULL) {
        iout << iWARN << "Non-NULL isBasePatch detected in SelfTuples::loadTuples()" << endi;
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

      this->tupleList.clear();

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
      // There should be only one!
      TuplePatchListIter ai(tuplePatchList);
      if (pids.size() == 0) ai = ai.begin();

      int numPid = (pids.size() == 0) ? tuplePatchList.size() : pids.size();

      for (int ipid=0;ipid < numPid;ipid++) {
        Patch *patch;
        int numAtoms;
        CompAtomExt *atomExt;
        // Take next patch
        if (pids.size() == 0) {
          patch = (*ai).p;
          numAtoms = patch->getNumAtoms();
          atomExt = (*ai).xExt;
          ai++;
        } else {
          TuplePatchElem *tpe = tuplePatchList.find(TuplePatchElem(pids[ipid]));
          patch = tpe->p;
          numAtoms = patch->getNumAtoms();
          atomExt = tpe->xExt;          
        }
  
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {
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
          /* cycle through each tuple */
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
            if (sdScaling && is_fep_sd) {   //check if this bonded term is one of Shobana term. This segment looks ugly and not GPU friendly, and might not appear in GPU code.
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
            if ( samepatch ) {
              t.scale = (!has_les && !has_ss) ? 1.0 : ( has_les ? invLesFactor : soluteScalingFactor );
              if (is_fep_ss) t.scale = (fep_tuple_type == 4) ? OneMinusLambda : Lambda;
              if (is_fep_sd && sdScaling) t.scale = (fep_tuple_type == 4 || fep_tuple_type == 2) ? OneMinusLambda : Lambda;
              TuplePatchElem *p;
              p = tuplePatchList.find(TuplePatchElem(homepatch));
              for(i=0; i < T::size; i++) {
                t.p[i] = p;
                t.localIndex[i] = aid[i].index;
              }
#ifdef MEM_OPT_VERSION
              //avoid adding Tuples whose atoms are all fixed
              if(node->simParameters->fixedAtomsOn &&
                 !node->simParameters->fixedAtomsForces) {
                int allfixed = 1;
                for(i=0; i<T::size; i++){
                  CompAtomExt *one = &(p->xExt[aid[i].index]);
                  allfixed = allfixed & one->atomFixed;
                }
                if(!allfixed) this->tupleList.push_back(t);
              } else {
                this->tupleList.push_back(t);
              }
#else
              this->tupleList.push_back(t);
#endif
            }
          }
        }
      }
    }

};
#endif

template <class T, class S, class P> class ComputeSelfTuples :
	public ComputeHomeTuples<T,S,P> {

#ifndef USE_HOMETUPLES
  private:
  
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

      this->tupleList.resize(0);

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
      // There should be only one!
      TuplePatchListIter ai(this->tuplePatchList);

      for ( ai = ai.begin(); ai != ai.end(); ai++ )
      {
    
        // CompAtomExt *atom = (*ai).x;
        Patch *patch = (*ai).p;
        int numAtoms = patch->getNumAtoms();
	CompAtomExt *atomExt = (*ai).xExt; //patch->getCompAtomExtInfo();
    
        // cycle through each atom in the patch and load up tuples
        for (int j=0; j < numAtoms; j++)
        {
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
           /* cycle through each tuple */
           for( ; *curTuple != -1; ++curTuple) {
             T t(&tupleStructs[*curTuple],tupleValues);
           #endif
             register int i;
             aid[0] = this->atomMap->localID(t.atomID[0]);
             int homepatch = aid[0].pid;
             int samepatch = 1;
             partition[0] = fepOn ? node->molecule->get_fep_type(t.atomID[0]) : 0; 
             int has_les = lesOn ? node->molecule->get_fep_type(t.atomID[0]) : 0;
             int has_ss = soluteScalingOn ? node->molecule->get_ss_type(t.atomID[0]) : 0;
             int is_fep_ss = partition[0] > 2;
             int is_fep_sd = 0;
             int fep_tuple_type = 0;
             for (i=1; i < T::size; i++) {
	         aid[i] = this->atomMap->localID(t.atomID[i]);
	         samepatch = samepatch && ( homepatch == aid[i].pid );
                 partition[i] = fepOn ? node->molecule->get_fep_type(t.atomID[i]) : 0;
                 has_les |= lesOn ? node->molecule->get_fep_type(t.atomID[i]) : 0;
                 has_ss |= soluteScalingOn ? node->molecule->get_ss_type(t.atomID[i]) : 0;
                 if (fepOn) {
                 is_fep_ss &= partition[i] > 2;
                 is_fep_sd |= (abs(partition[i] - partition[0]) == 2);
                 fep_tuple_type = partition[i]; }
             }
             if (T::size < 4 && !soluteScalingAll) has_ss = false;
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
             if ( samepatch ) {
               t.scale = (!has_les && !has_ss) ? 1.0 : ( has_les ? invLesFactor : soluteScalingFactor );
               if (is_fep_ss) t.scale = (fep_tuple_type == 4) ? OneMinusLambda : Lambda;
               if (is_fep_sd && sdScaling) t.scale = (fep_tuple_type == 4 || fep_tuple_type == 2) ? OneMinusLambda : Lambda;
               TuplePatchElem *p;
               p = this->tuplePatchList.find(TuplePatchElem(homepatch));
               for(i=0; i < T::size; i++) {
                   t.p[i] = p;
                   t.localIndex[i] = aid[i].index;
               }
             #ifdef MEM_OPT_VERSION
               //avoid adding Tuples whose atoms are all fixed
	       if(node->simParameters->fixedAtomsOn &&
                  !node->simParameters->fixedAtomsForces) {
   		 int allfixed = 1;
                 for(i=0; i<T::size; i++){
                   CompAtomExt *one = &(p->xExt[aid[i].index]);
                   allfixed = allfixed & one->atomFixed;
                 }
                 if(!allfixed) this->tupleList.add(t);
               }else{
                 this->tupleList.add(t);
               }
             #else
               this->tupleList.add(t);
             #endif               
             }
           }
        }
      }
    }
#endif

    PatchID patchID;

  public:

    ComputeSelfTuples(ComputeID c, PatchID p) : ComputeHomeTuples<T,S,P>(c) {
      patchID = p;
    }

    virtual ~ComputeSelfTuples() {
      UniqueSetIter<TuplePatchElem> ap(this->tuplePatchList);
      for (ap = ap.begin(); ap != ap.end(); ap++) {
        ap->p->unregisterPositionPickup(this,&(ap->positionBox));
        ap->p->unregisterAvgPositionPickup(this,&(ap->avgPositionBox));
        ap->p->unregisterForceDeposit(this,&(ap->forceBox));
      }
    }


    //======================================================================
    // initialize() - Method is invoked only the first time
    // atom maps, patchmaps etc are ready and we are about to start computations
    //======================================================================
    virtual void initialize(void) {
#ifdef USE_HOMETUPLES
      this->tuples = new SelfTuples<T, S, P>();
#endif    
      // Start with empty list
      this->tuplePatchList.clear();
    
      this->tuplePatchList.add(TuplePatchElem(ComputeHomeTuples<T,S,P>::patchMap->patch(patchID), this));
    
      this->setNumPatches(this->tuplePatchList.size());

      this->doLoadTuples = true;

      int myNode = CkMyPe();
      if ( PatchMap::Object()->node(patchID) != myNode )
      {
        this->basePriority = COMPUTE_PROXY_PRIORITY + PATCH_PRIORITY(patchID);
      }
      else
      {
        this->basePriority = COMPUTE_HOME_PRIORITY + PATCH_PRIORITY(patchID);
      }
    }

    void doWork(void) {
//      LdbCoordinator::Object()->startWork(this->ldObjHandle);

#ifdef TRACE_COMPUTE_OBJECTS
    double traceObjStartTime = CmiWallTimer();
#endif

      ComputeHomeTuples<T,S,P>::doWork();

#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+this->cid, traceObjStartTime, CmiWallTimer());
#endif

//      LdbCoordinator::Object()->endWork(this->ldObjHandle);
    }

};


#endif

