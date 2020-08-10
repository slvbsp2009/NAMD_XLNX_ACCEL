/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeGlobal.h"
#include "ComputeGlobalMsgs.h"
#include "GridForceGrid.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "SimParameters.h"
#include <stdio.h>
#include <algorithm>

#include "Debug.h"

#include "GridForceGrid.inl"
#include "MGridforceParams.h"

// CLIENTS

ComputeGlobal::ComputeGlobal(ComputeID c, ComputeMgr *m)
  : ComputeHomePatches(c)
{
  DebugM(3,"Constructing client\n");
  aid.resize(0);
  gdef.resize(0);
  comm = m;
  firsttime = 1;
  isRequested = 0;
  isRequestedAllocSize = 0;
  endRequested = 0;
  numGroupsRequested = 0;
  SimParameters *sp = Node::Object()->simParameters;
  dofull = (sp->GBISserOn || sp->GBISOn || sp->fullDirectOn || sp->FMAOn || sp->PMEOn);
  forceSendEnabled = 0;
  if ( sp->tclForcesOn ) forceSendEnabled = 1;
  if ( sp->colvarsOn ) forceSendEnabled = 1;
  forceSendActive = 0;
  fid.resize(0);
  totalForce.resize(0);
  gfcount = 0;
  groupTotalForce.resize(0);
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  int numPatches = PatchMap::Object()->numPatches();
  forcePtrs = new Force*[numPatches];
  atomPtrs = new FullAtom*[numPatches];
  gridForcesPtrs = new ForceList **[numPatches];
  numGridObjects = numActiveGridObjects = 0;
  for ( int i = 0; i < numPatches; ++i ) {
    forcePtrs[i] = NULL; atomPtrs[i] = NULL;
    gridForcesPtrs[i] = NULL;
  }
}

ComputeGlobal::~ComputeGlobal()
{
  delete[] isRequested;
  delete[] forcePtrs;
  deleteGridObjects();
  delete[] gridForcesPtrs;
  delete[] atomPtrs;
  delete reduction;
}

void ComputeGlobal::configure(AtomIDList &newaid, AtomIDList &newgdef, IntList &newgridobjid) {
  DebugM(4,"Receiving configuration (" << newaid.size() <<
         " atoms, " << newgdef.size() << " atoms/groups and " <<
         newgridobjid.size() << " grid objects) on client\n" << endi);

  AtomIDList::iterator a, a_e;

  if ( forceSendEnabled ) {
    // clear previous data
    int max = -1;
    for (a=newaid.begin(),a_e=newaid.end(); a!=a_e; ++a) {
      if ( *a > max ) max = *a;
    }
    for (a=newgdef.begin(),a_e=newgdef.end(); a!=a_e; ++a) {
      if ( *a > max ) max = *a;
    }
    endRequested = max+1;
    if ( endRequested > isRequestedAllocSize ) {
      delete [] isRequested;
      isRequestedAllocSize = endRequested+10;
      isRequested = new char[isRequestedAllocSize];
      memset(isRequested, 0, isRequestedAllocSize);
    } else {
      for (a=aid.begin(),a_e=aid.end(); a!=a_e; ++a) {
        isRequested[*a] = 0;
      }
      for (a=gdef.begin(),a_e=gdef.end(); a!=a_e; ++a) {
        if ( *a != -1 ) isRequested[*a] = 0;
      }
    }
    // reserve space
    gpair.resize(0);
    gpair.resize(newgdef.size());
    gpair.resize(0);
  }

  // store data
  aid.swap(newaid);
  gdef.swap(newgdef);

  if (newgridobjid.size()) configureGridObjects(newgridobjid);

  if ( forceSendEnabled ) {
    int newgcount = 0;
    for (a=aid.begin(),a_e=aid.end(); a!=a_e; ++a) {
      isRequested[*a] = 1;
    }
    for (a=gdef.begin(),a_e=gdef.end(); a!=a_e; ++a) {
      if ( *a == -1 ) ++newgcount;
      else {
        isRequested[*a] |= 2;
        gpair.add(intpair(*a,newgcount));
      }
    }
    std::sort(gpair.begin(),gpair.end());
    numGroupsRequested = newgcount;
  }
}

void ComputeGlobal::deleteGridObjects()
{
  if (numGridObjects == 0) return;
  ResizeArrayIter<PatchElem> ap(patchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    ForceList **gridForces = gridForcesPtrs[ap->p->getPatchID()];
    if (gridForces != NULL) {
      for (size_t ig = 0; ig < numGridObjects; ig++) {
        if (gridForces[ig] != NULL) {
          delete gridForces[ig];
          gridForces[ig] = NULL;
        }
      }
      delete [] gridForces;
      gridForces = NULL;
    }
  }
  numGridObjects = numActiveGridObjects = 0;
}

void ComputeGlobal::configureGridObjects(IntList &newgridobjid)
{
  Molecule *mol = Node::Object()->molecule;

  deleteGridObjects();

  numGridObjects = mol->numGridforceGrids;
  numActiveGridObjects = 0;

  gridObjActive.resize(numGridObjects);
  gridObjActive.setall(0);

  IntList::const_iterator goid_i = newgridobjid.begin();
  IntList::const_iterator goid_e = newgridobjid.end();
  for ( ; goid_i != goid_e; goid_i++) {
    if ((*goid_i < 0) || (*goid_i >= numGridObjects)) {
      NAMD_bug("Requested illegal gridForceGrid index.");
    } else {
      DebugM(3,"Adding grid with index " << *goid_i << " to ComputeGlobal\n");
      gridObjActive[*goid_i] = 1;
      numActiveGridObjects++;
    }
  }

  for (size_t ig = 0; ig < numGridObjects; ig++) {
    DebugM(3,"Grid index " << ig << " is active or inactive? "
           << gridObjActive[ig] << "\n" << endi);
  }

  ResizeArrayIter<PatchElem> ap(patchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    gridForcesPtrs[ap->p->getPatchID()] = new ForceList *[numGridObjects];
    ForceList **gridForces = gridForcesPtrs[ap->p->getPatchID()];
    for (size_t ig = 0; ig < numGridObjects; ig++) {
      if (gridObjActive[ig]) {
        gridForces[ig] = new ForceList;
      } else {
        gridForces[ig] = NULL;
      }
    }
  }
}

#if 0
void ComputeGlobal::recvConfig(ComputeGlobalConfigMsg *msg) {
  DebugM(3,"Receiving configure on client\n");
  configure(msg->aid,msg->gdef);
  delete msg;
  sendData();
}
#endif

void ComputeGlobal::recvResults(ComputeGlobalResultsMsg *msg) {
  DebugM(3,"Receiving results (" << msg->aid.size() << " forces, "
	 << msg->newgdef.size() << " new group atoms) on client\n");

  forceSendActive = msg->totalforces;
  if ( forceSendActive && ! forceSendEnabled ) NAMD_bug("ComputeGlobal::recvResults forceSendActive without forceSendEnabled");

  // set the forces only if we aren't going to resend the data
  int setForces = !msg->resendCoordinates;

  if(setForces) { // we are requested to 
    // Store forces to patches
    AtomMap *atomMap = AtomMap::Object();
    const Lattice & lattice = patchList[0].p->lattice;
    ResizeArrayIter<PatchElem> ap(patchList);
    Force **f = forcePtrs;
    FullAtom **t = atomPtrs;
    Force extForce = 0.;
    Tensor extVirial;

    for (ap = ap.begin(); ap != ap.end(); ap++) {
      (*ap).r = (*ap).forceBox->open();
      f[(*ap).patchID] = (*ap).r->f[Results::normal];
      t[(*ap).patchID] = (*ap).p->getAtomList().begin();
    }

    AtomIDList::iterator a = msg->aid.begin();
    AtomIDList::iterator a_e = msg->aid.end();
    ForceList::iterator f2 = msg->f.begin();
    for ( ; a != a_e; ++a, ++f2 ) {
      DebugM(1,"processing atom "<<(*a)<<", F="<<(*f2)<<"...\n");
      /* XXX if (*a) is out of bounds here we get a segfault */
      LocalID localID = atomMap->localID(*a);
      if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
      Force f_atom = (*f2);
      f[localID.pid][localID.index] += f_atom;
      FullAtom &atom = t[localID.pid][localID.index];
      Position x_orig = atom.position;
      Transform trans = atom.transform;
      Position x_atom = lattice.reverse_transform(x_orig,trans);
      extForce += f_atom;
      extVirial += outer(f_atom,x_atom);
    }
    DebugM(1,"done with the loop\n");

    // calculate forces for atoms in groups
    AtomIDList::iterator g_i, g_e;
    g_i = gdef.begin(); g_e = gdef.end();
    ForceList::iterator gf_i = msg->gforce.begin();
    //iout << iDEBUG << "recvResults\n" << endi;
    for ( ; g_i != g_e; ++g_i, ++gf_i ) {
      //iout << iDEBUG << *gf_i << '\n' << endi;
      Vector accel = (*gf_i);
      for ( ; *g_i != -1; ++g_i ) {
	//iout << iDEBUG << *g_i << '\n' << endi;
	LocalID localID = atomMap->localID(*g_i);
	if ( localID.pid == notUsed || ! f[localID.pid] ) continue;
        FullAtom &atom = t[localID.pid][localID.index];
	Force f_atom = accel * atom.mass;
	f[localID.pid][localID.index] += f_atom;
        Position x_orig = atom.position;
        Transform trans = atom.transform;
        Position x_atom = lattice.reverse_transform(x_orig,trans);
        extForce += f_atom;
        extVirial += outer(f_atom,x_atom);
      }
    }
    DebugM(1,"done with the groups\n");

    if (numActiveGridObjects > 0) {
      applyGridObjectForces(msg, &extForce, &extVirial);
    }

    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,extForce);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,extVirial);
    reduction->submit();
  }
  // done setting the forces, close boxes below

  // Get reconfiguration if present
  if ( msg->reconfig ) {
    DebugM(3,"Reconfiguring\n");
    configure(msg->newaid, msg->newgdef, msg->newgridobjid);
  }

  // send another round of data if requested

  if(msg->resendCoordinates) {
    DebugM(3,"Sending requested data right away\n");
    sendData();
  }

  groupTotalForce.resize(numGroupsRequested);
  for ( int i=0; i<numGroupsRequested; ++i ) groupTotalForce[i] = 0;

  if(setForces) {
    ResizeArrayIter<PatchElem> ap(patchList);
    Force **f = forcePtrs;
    FullAtom **t = atomPtrs;
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      CompAtom *x;
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&((*ap).r));
      f[(*ap).patchID] = 0;
      t[(*ap).patchID] = 0;
    }
  }

  delete msg;
  DebugM(3,"Done processing results\n");
}

void ComputeGlobal::doWork()
{
  DebugM(2,"doWork\n");

  ResizeArrayIter<PatchElem> ap(patchList);
  FullAtom **t = atomPtrs;

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    CompAtom *x = (*ap).positionBox->open();
    t[(*ap).patchID] = (*ap).p->getAtomList().begin();
  }

  if(!firsttime) sendData();
  else {
    if ( hasPatchZero ) {
      ComputeGlobalDataMsg *msg = new ComputeGlobalDataMsg;
      msg->lat.add(patchList[0].p->lattice);
      msg->step = -1;
      msg->count = 1;
      msg->patchcount = 0;
      comm->sendComputeGlobalData(msg);
    }
    firsttime = 0;
    comm->enableComputeGlobalResults();
  }
  DebugM(2,"done with doWork\n");
}

void ComputeGlobal::sendData()
{
  DebugM(2,"sendData\n");
  // Get positions from patches
  AtomMap *atomMap = AtomMap::Object();
  const Lattice & lattice = patchList[0].p->lattice;
  FullAtom **t = atomPtrs;

  ComputeGlobalDataMsg *msg = new  ComputeGlobalDataMsg;

  msg->step = patchList[0].p->flags.step;
  msg->count = 0;
  msg->patchcount = 0;

  AtomIDList::iterator a = aid.begin();
  AtomIDList::iterator a_e = aid.end();
  for ( ; a != a_e; ++a ) {
    LocalID localID = atomMap->localID(*a);
    if ( localID.pid == notUsed || ! t[localID.pid] ) continue;
    msg->aid.add(*a);
    msg->count++;
    FullAtom &atom = t[localID.pid][localID.index];
    Position x_orig = atom.position;
    Transform trans = atom.transform;
    msg->p.add(lattice.reverse_transform(x_orig,trans));
  }

  // calculate group centers of mass
  AtomIDList::iterator g_i, g_e;
  g_i = gdef.begin(); g_e = gdef.end();
  for ( ; g_i != g_e; ++g_i ) {
    Vector com(0,0,0);
    BigReal mass = 0.;
    for ( ; *g_i != -1; ++g_i ) {
      LocalID localID = atomMap->localID(*g_i);
      if ( localID.pid == notUsed || ! t[localID.pid] ) continue;
      msg->count++;
      FullAtom &atom = t[localID.pid][localID.index];
      Position x_orig = atom.position;
      Transform trans = atom.transform;
      com += lattice.reverse_transform(x_orig,trans) * atom.mass;
      mass += atom.mass;
    }
    DebugM(1,"Adding center of mass "<<com<<"\n");
    msg->gcom.add(com);
    msg->gmass.add(mass);
  }

  if (numActiveGridObjects > 0) {
    computeGridObjects(msg);
  }

  msg->fid.swap(fid);
  msg->tf.swap(totalForce);
  fid.resize(0);
  totalForce.resize(0);

  if ( gfcount ) msg->gtf.swap(groupTotalForce);
  msg->count += ( msg->fid.size() + gfcount );
  gfcount = 0;

  DebugM(3,"Sending data (" << msg->p.size() << " positions, "
         << msg->gcom.size() << " groups, " << msg->gridobjvalue.size()
         << " grid objects) on client\n");
  if ( hasPatchZero ) { msg->count++;  msg->lat.add(lattice); }
  if ( msg->count || msg->patchcount ) comm->sendComputeGlobalData(msg);
  else delete msg;
  comm->enableComputeGlobalResults();
}

template<class T> void ComputeGlobal::computeGridForceGrid(FullAtomList::iterator aii,
                                                           FullAtomList::iterator aei,
                                                           ForceList::iterator gfii,
                                                           Lattice const &lattice,
                                                           int gridIndex,
                                                           T *grid,
                                                           BigReal &gridObjValue)
{
  ForceList::iterator gfi = gfii;
  FullAtomList::iterator ai = aii;
  FullAtomList::iterator ae = aei;
  Molecule *mol = Node::Object()->molecule;
  for ( ; ai != ae; ai++, gfi++) {
    *gfi = Vector(0.0, 0.0, 0.0);
    if (! mol->is_atom_gridforced(ai->id, gridIndex)) {
      continue;
    }
    Real scale;
    Charge charge;
    Vector dV;
    float V;
    mol->get_gridfrc_params(scale, charge, ai->id, gridIndex);
    Position pos = grid->wrap_position(ai->position, lattice);
    DebugM(1, "id = " << ai->id << ", scale = " << scale
           << ", charge = " << charge << ", position = " << pos << "\n");
    if (grid->compute_VdV(pos, V, dV)) {
      // out-of-bounds atom
      continue;
    }
    // ignore global gfScale
    *gfi = -charge * scale * dV;
    gridObjValue += charge * scale * V;
    DebugM(1, "id = " << ai->id << ", force = " << *gfi << "\n");
  }
  DebugM(3, "gridObjValue = " << gridObjValue << "\n" << endi);
}

void ComputeGlobal::computeGridObjects(ComputeGlobalDataMsg *msg)
{
  DebugM(3,"computeGridObjects\n" << endi);
  Molecule *mol = Node::Object()->molecule;
  const Lattice &lattice = patchList[0].p->lattice;

  if (mol->numGridforceGrids < 1) {
    NAMD_bug("No grids loaded in memory but ComputeGlobal has been requested to use them.");
  }

  msg->gridobjindex.resize(numActiveGridObjects);
  msg->gridobjindex.setall(-1);
  msg->gridobjvalue.resize(numActiveGridObjects);
  msg->gridobjvalue.setall(0.0);

  size_t ig = 0, gridobjcount = 0;

  // loop over home patches
  ResizeArrayIter<PatchElem> ap(patchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {

    msg->patchcount++;

    int const numAtoms = ap->p->getNumAtoms();
    ForceList **gridForces = gridForcesPtrs[ap->p->getPatchID()];

    gridobjcount = 0;
    for (ig = 0; ig < numGridObjects; ig++) {

      DebugM(2,"Processing grid index " << ig << "\n" << endi);

      // Only process here objects requested by the GlobalMasters
      if (!gridObjActive[ig]) {
        DebugM(2,"Skipping grid index " << ig << "; it is handled by "
               "ComputeGridForce\n" << endi);
        continue;
      }

      ForceList *gridForcesGrid = gridForces[ig];
      gridForcesGrid->resize(numAtoms);

      ForceList::iterator gfi = gridForcesGrid->begin();
      FullAtomList::iterator ai = ap->p->getAtomList().begin();
      FullAtomList::iterator ae = ap->p->getAtomList().end();

      DebugM(2, "computeGridObjects(): patch = " << ap->p->getPatchID()
             << ", grid index = " << ig << "\n" << endi);
      GridforceGrid *grid = mol->get_gridfrc_grid(ig);

      msg->gridobjindex[gridobjcount] = ig;
      BigReal &gridobjvalue = msg->gridobjvalue[gridobjcount];

      if (grid->get_grid_type() == GridforceGrid::GridforceGridTypeFull) {

        GridforceFullMainGrid *g = dynamic_cast<GridforceFullMainGrid *>(grid);
        computeGridForceGrid(ai, ae, gfi, ap->p->lattice, ig, g, gridobjvalue);

      } else if (grid->get_grid_type() == GridforceGrid::GridforceGridTypeLite) {

        GridforceLiteGrid *g = dynamic_cast<GridforceLiteGrid *>(grid);
        computeGridForceGrid(ai, ae, gfi, ap->p->lattice, ig, g, gridobjvalue);
      }

      gridobjcount++;
    }
  }

  for (gridobjcount = 0; gridobjcount < numActiveGridObjects; gridobjcount++) {
    DebugM(3, "Total gridObjValue[" << msg->gridobjindex[gridobjcount]
           << "] = " << msg->gridobjvalue[gridobjcount] << "\n");
  }

  DebugM(2,"computeGridObjects done\n");
}

void ComputeGlobal::applyGridObjectForces(ComputeGlobalResultsMsg *msg,
                                          Force *extForce_in,
                                          Tensor *extVirial_in)
{
  if (msg->gridobjforce.size() == 0) return;

  if (msg->gridobjforce.size() != numActiveGridObjects) {
    NAMD_bug("ComputeGlobal received a different number of grid forces than active grids.");
  }

  Molecule *mol = Node::Object()->molecule;
  const Lattice &lattice = patchList[0].p->lattice;
  AtomMap *atomMap = AtomMap::Object();
  Force &extForce = *extForce_in;
  Tensor &extVirial = *extVirial_in;

  // map applied forces from the message
  BigRealList gridObjForces;
  gridObjForces.resize(numGridObjects);
  gridObjForces.setall(0.0);
  BigRealList::iterator gridobjforce_i = msg->gridobjforce.begin();
  BigRealList::iterator gridobjforce_e = msg->gridobjforce.end();
  int ig;
  for (ig = 0; gridobjforce_i != gridobjforce_e ;
       gridobjforce_i++, ig++) {
    if (!gridObjActive[ig]) continue;
    gridObjForces[ig] = *gridobjforce_i;
  }

  // loop over home patches
  ResizeArrayIter<PatchElem> ap(patchList);
  for (ap = ap.begin(); ap != ap.end(); ap++) {

    ForceList **gridForces = gridForcesPtrs[ap->p->getPatchID()];

    for (ig = 0; ig < numGridObjects; ig++) {

      if (!gridObjActive[ig]) continue;

      DebugM(2, "gof  = " << gridObjForces[ig] << "\n" << endi);

      ForceList *gridForcesGrid = gridForces[ig];

      FullAtomList::iterator ai = ap->p->getAtomList().begin();
      FullAtomList::iterator ae = ap->p->getAtomList().end();
      Force *f = ap->r->f[Results::normal];
      ForceList::iterator gfi = gridForcesGrid->begin();

      for ( ; ai != ae; ai++, gfi++) {
        if (! mol->is_atom_gridforced(ai->id, ig)) {
          *gfi = Vector(0.0, 0.0, 0.0);
          continue;
        }
	LocalID localID = atomMap->localID(ai->id);
        // forces were stored; flipping sign to get gradients
        Vector const gridforceatom(-1.0 * (*gfi) * gridObjForces[ig]);
        DebugM(2, "id = " << ai->id
               << ", pid = " << localID.pid
               << ", index = " << localID.index
               << ", force = " << gridforceatom << "\n" << endi);
        f[localID.index] += gridforceatom;
        extForce += gridforceatom;
        Position x_orig = ai->position;
        Transform transform = ai->transform;
        Position x_virial = lattice.reverse_transform(x_orig, transform);
        extVirial += outer(gridforceatom, x_virial);
      }
    }
  }
  // extForce and extVirial are being communicated by calling function
}

// This function is called by each HomePatch after force
// evaluation. It stores the indices and forces of the requested
// atoms here, to be sent to GlobalMasterServer during the next
// time step. The total force is the sum of three components:
// "normal", "nbond" and "slow", the latter two may be calculated
// less frequently, so their most recent values are stored in
// "f_saved" and used here. If we don't do full electrostatics,
// there's no "slow" part.
void ComputeGlobal::saveTotalForces(HomePatch *homePatch)
{
  if ( ! forceSendEnabled ) NAMD_bug("ComputeGlobal::saveTotalForces called unexpectedly");
  if ( ! forceSendActive ) return;

  if ( Node::Object()->simParameters->accelMDOn && Node::Object()->simParameters->accelMDDebugOn && Node::Object()->simParameters->accelMDdihe ) {
    int num=homePatch->numAtoms;
    FullAtomList &atoms = homePatch->atom;
    ForceList &af=homePatch->f[Results::amdf];

    for (int i=0; i<num; ++i) {
      int index = atoms[i].id;
      if (index < endRequested && isRequested[index] & 1) {
        fid.add(index);
        totalForce.add(af[i]);
      }
    }
    return;
  }

  int fixedAtomsOn = Node::Object()->simParameters->fixedAtomsOn;
  int num=homePatch->numAtoms;
  FullAtomList &atoms = homePatch->atom;
  ForceList &f1=homePatch->f[Results::normal], &f2=homePatch->f_saved[Results::nbond],
            &f3=homePatch->f_saved[Results::slow];
  Force f_sum;
  
  for (int i=0; i<num; ++i) {
    int index = atoms[i].id;
    char reqflag;
    if (index < endRequested && (reqflag = isRequested[index])) {
     f_sum = f1[i]+f2[i];
     if (dofull)
       f_sum += f3[i];
     if ( fixedAtomsOn && atoms[i].atomFixed ) f_sum = 0.;
     if ( reqflag  & 1 ) {  // individual atom
      fid.add(index);
      totalForce.add(f_sum);
     }
     if ( reqflag  & 2 ) {  // part of group
       intpair *gpend = gpair.end();
       intpair *gpi = std::lower_bound(gpair.begin(),gpend,intpair(index,0));
       if ( gpi == gpend || gpi->first != index )
         NAMD_bug("ComputeGlobal::saveTotalForces gpair corrupted.");
       do {
         ++gfcount;
         groupTotalForce[gpi->second] += f_sum;
       } while ( ++gpi != gpend && gpi->first == index );
     }
    }
  }
}
