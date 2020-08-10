/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/Sequencer.C,v $
 * $Author: jim $
 * $Date: 2016/08/26 19:40:32 $
 * $Revision: 1.1230 $
 *****************************************************************************/

// The UPPER_BOUND macro is used to eliminate all of the per atom
// computation done for the numerical integration in Sequencer::integrate()
// other than the actual force computation and atom migration.
// The idea is to "turn off" the integration for doing performance
// profiling in order to get an upper bound on the speedup available
// by moving the integration parts to the GPU.
//
// Define it in the Make.config file, i.e. CXXOPTS += -DUPPER_BOUND
// or simply uncomment the line below.
//
//#define UPPER_BOUND

//for gbis debugging; print net force on each atom
#define PRINT_FORCES 0

#include "InfoStream.h"
#include "Node.h"
#include "SimParameters.h"
#include "Sequencer.h"
#include "HomePatch.h"
#include "ReductionMgr.h"
#include "CollectionMgr.h"
#include "BroadcastObject.h"
#include "Output.h"
#include "Controller.h"
#include "Broadcasts.h"
#include "Molecule.h"
#include "NamdOneTools.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include "Random.h"
#include "PatchMap.inl"
#include "ComputeMgr.h"
#include "ComputeGlobal.h"
#include "NamdEventsProfiling.h"

#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
//
// Define NL_DEBUG below to activate D_*() macros in integrate_SOA()
// for debugging.
//
//#define NL_DEBUG
#include "Debug.h"

#if USE_HPM
#define START_HPM_STEP  1000
#define STOP_HPM_STEP   1500
#endif

#define SPECIAL_PATCH_ID  91

Sequencer::Sequencer(HomePatch *p) :
	simParams(Node::Object()->simParameters),
	patch(p),
	collection(CollectionMgr::Object()),
	ldbSteps(0)
{
    broadcast = new ControllerBroadcasts(& patch->ldObjHandle);
    reduction = ReductionMgr::Object()->willSubmit(
                  simParams->accelMDOn ? REDUCTIONS_AMD : REDUCTIONS_BASIC );
    min_reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_MINIMIZER,1);
    if (simParams->pressureProfileOn) {
      int ntypes = simParams->pressureProfileAtomTypes;
      int nslabs = simParams->pressureProfileSlabs;
      pressureProfileReduction = 
        ReductionMgr::Object()->willSubmit(
		REDUCTIONS_PPROF_INTERNAL, 3*nslabs*ntypes);
    } else {
      pressureProfileReduction = NULL;
    }
    if (simParams->multigratorOn) {
      multigratorReduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_MULTIGRATOR,MULTIGRATOR_REDUCTION_MAX_RESERVED);
    } else {
      multigratorReduction = NULL;
    }
    ldbCoordinator = (LdbCoordinator::Object());
    random = new Random(simParams->randomSeed);
    random->split(patch->getPatchID()+1,PatchMap::Object()->numPatches()+1);

    // Is soluteScaling enabled?
    if (simParams->soluteScalingOn) {
      // If so, we must "manually" perform charge scaling on startup because
      // Sequencer will not get a scripting task for initial charge scaling.
      // Subsequent rescalings will take place through a scripting task.
      rescaleSoluteCharges(simParams->soluteScalingFactorCharge);
    }

    rescaleVelocities_numTemps = 0;
    stochRescale_count = 0;
    berendsenPressure_count = 0;
//    patch->write_tip4_props();
}

Sequencer::~Sequencer(void)
{
    delete broadcast;
    delete reduction;
    delete min_reduction;
    if (pressureProfileReduction) delete pressureProfileReduction;
    delete random;
    if (multigratorReduction) delete multigratorReduction;
}

// Invoked by thread
void Sequencer::threadRun(Sequencer* arg)
{
    LdbCoordinator::Object()->startWork(arg->patch->ldObjHandle);
    arg->algorithm();
}

// Invoked by Node::run() via HomePatch::runSequencer()
void Sequencer::run(void)
{
    // create a Thread and invoke it
    DebugM(4, "::run() - this = " << this << "\n" );
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),SEQ_STK_SZ);
    CthSetStrategyDefault(thread);
    priority = PATCH_PRIORITY(patch->getPatchID());
    awaken();
}

void Sequencer::suspend(void)
{
    LdbCoordinator::Object()->pauseWork(patch->ldObjHandle);
    CthSuspend();
    LdbCoordinator::Object()->startWork(patch->ldObjHandle);
}

// Defines sequence of operations on a patch.  e.g. when
// to push out information for Compute objects to consume
// when to migrate atoms, when to add forces to velocity update.
void Sequencer::algorithm(void)
{
  int scriptTask;
  int scriptSeq = 0;
  // Blocking receive for the script barrier.
  while ( (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END ) {
    switch ( scriptTask ) {
      case SCRIPT_OUTPUT:
	submitCollections(FILE_OUTPUT);
	break;
      case SCRIPT_FORCEOUTPUT:
	submitCollections(FORCE_OUTPUT);
	break;
      case SCRIPT_MEASURE:
	submitCollections(EVAL_MEASURE);
	break;
      case SCRIPT_REINITVELS:
	reinitVelocities();
	break;
      case SCRIPT_RESCALEVELS:
	rescaleVelocitiesByFactor(simParams->scriptArg1);
	break;
      case SCRIPT_RESCALESOLUTECHARGES:
        rescaleSoluteCharges(simParams->soluteScalingFactorCharge);
        break;
      case SCRIPT_RELOADCHARGES:
	reloadCharges();
	break;
      case SCRIPT_CHECKPOINT:
        patch->checkpoint();
        checkpoint_berendsenPressure_count = berendsenPressure_count;
	break;
      case SCRIPT_REVERT:
        patch->revert();
        berendsenPressure_count = checkpoint_berendsenPressure_count;
        pairlistsAreValid = 0;
	break;
      case SCRIPT_CHECKPOINT_STORE:
      case SCRIPT_CHECKPOINT_LOAD:
      case SCRIPT_CHECKPOINT_SWAP:
      case SCRIPT_CHECKPOINT_FREE:
        patch->exchangeCheckpoint(scriptTask,berendsenPressure_count);
	break;
      case SCRIPT_ATOMSENDRECV:
      case SCRIPT_ATOMSEND:
      case SCRIPT_ATOMRECV:
        patch->exchangeAtoms(scriptTask);
        break;
      case SCRIPT_MINIMIZE:
	minimize();
	break;
      case SCRIPT_RUN:
      case SCRIPT_CONTINUE:
  //
  // DJH: Call a cleaned up version of integrate().
  //
  // We could test for simulation options and call a more basic version
  // of integrate() where we can avoid performing most tests.
  //
	integrate(scriptTask);
	break;
      default:
        NAMD_bug("Unknown task in Sequencer::algorithm");
    }
  }
  submitCollections(END_OF_RUN);
  terminate();
}


extern int eventEndOfTimeStep;

void Sequencer::integrate(int scriptTask) {
    char traceNote[24];
    char tracePrefix[20];
    sprintf(tracePrefix, "p:%d,s:",patch->patchID);
//    patch->write_tip4_props();

    //
    // DJH: Copy all data into SOA (structure of arrays)
    // from AOS (array of structures) data structure.
    //
    //patch->copy_all_to_SOA();

#ifdef TIMER_COLLECTION
    TimerSet& t = patch->timerSet;
#endif
    TIMER_INIT_WIDTH(t, KICK, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, MAXMOVE, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, DRIFT, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, PISTON, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, SUBMITHALF, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, VELBBK1, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, VELBBK2, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, RATTLE1, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, SUBMITFULL, simParams->timerBinWidth);
    TIMER_INIT_WIDTH(t, SUBMITCOLLECT, simParams->timerBinWidth);

    int &step = patch->flags.step;
    step = simParams->firstTimestep;

    // drag switches
    const Bool rotDragOn = simParams->rotDragOn;
    const Bool movDragOn = simParams->movDragOn;

    const int commOnly = simParams->commOnly;

    int &maxForceUsed = patch->flags.maxForceUsed;
    int &maxForceMerged = patch->flags.maxForceMerged;
    maxForceUsed = Results::normal;
    maxForceMerged = Results::normal;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    const BigReal timestep = simParams->dt;

    // what MTS method?
    const int staleForces = ( simParams->MTSAlgorithm == NAIVE );

    const int nonbondedFrequency = simParams->nonbondedFrequency;
    slowFreq = nonbondedFrequency;
    const BigReal nbondstep = timestep * (staleForces?1:nonbondedFrequency);
    int &doNonbonded = patch->flags.doNonbonded;
    doNonbonded = (step >= numberOfSteps) || !(step%nonbondedFrequency);
    if ( nonbondedFrequency == 1 ) maxForceMerged = Results::nbond;
    if ( doNonbonded ) maxForceUsed = Results::nbond;

    // Do we do full electrostatics?
    const int dofull = ( simParams->fullElectFrequency ? 1 : 0 );
    const int fullElectFrequency = simParams->fullElectFrequency;
    if ( dofull ) slowFreq = fullElectFrequency;
    const BigReal slowstep = timestep * (staleForces?1:fullElectFrequency);
    int &doFullElectrostatics = patch->flags.doFullElectrostatics;
    doFullElectrostatics = (dofull && ((step >= numberOfSteps) || !(step%fullElectFrequency)));
    if ( dofull && (fullElectFrequency == 1) && !(simParams->mollyOn) )
					maxForceMerged = Results::slow;
    if ( doFullElectrostatics ) maxForceUsed = Results::slow;

//#ifndef UPPER_BOUND
    const Bool accelMDOn = simParams->accelMDOn;
    const Bool accelMDdihe = simParams->accelMDdihe;
    const Bool accelMDdual = simParams->accelMDdual;
    if ( accelMDOn && (accelMDdihe || accelMDdual)) maxForceUsed = Results::amdf;

    // Is adaptive tempering on?
    const Bool adaptTempOn = simParams->adaptTempOn;
    adaptTempT = simParams->initialTemp;
    if (simParams->langevinOn)
        adaptTempT = simParams->langevinTemp;
    else if (simParams->rescaleFreq > 0)
        adaptTempT = simParams->rescaleTemp;
        

    int &doMolly = patch->flags.doMolly;
    doMolly = simParams->mollyOn && doFullElectrostatics;
    // BEGIN LA
    int &doLoweAndersen = patch->flags.doLoweAndersen;
    doLoweAndersen = simParams->loweAndersenOn && doNonbonded;
    // END LA

    int &doGBIS = patch->flags.doGBIS;
    doGBIS = simParams->GBISOn;

    int &doLCPO = patch->flags.doLCPO;
    doLCPO = simParams->LCPOOn;

    int zeroMomentum = simParams->zeroMomentum;
    
    // Do we need to return forces to TCL script or Colvar module?
    int doTcl = simParams->tclForcesOn;
	int doColvars = simParams->colvarsOn;
//#endif
    ComputeGlobal *computeGlobal = Node::Object()->computeMgr->computeGlobalObject;

    // Bother to calculate energies?
    int &doEnergy = patch->flags.doEnergy;
    int energyFrequency = simParams->outputEnergies;
#ifndef UPPER_BOUND
    const int reassignFreq = simParams->reassignFreq;
#endif

    int &doVirial = patch->flags.doVirial;
    doVirial = 1;

  if ( scriptTask == SCRIPT_RUN ) {

#ifndef UPPER_BOUND
//    printf("Doing initial rattle\n");
#ifndef UPPER_BOUND
D_MSG("rattle1()");
    TIMER_START(t, RATTLE1);
    rattle1(0.,0);  // enforce rigid bond constraints on initial positions
    TIMER_STOP(t, RATTLE1);
#endif

    if (simParams->lonepairs || simParams->singleTopology) {
      patch->atomMapper->registerIDsFullAtom(
		patch->atom.begin(),patch->atom.end());
    }

    if ( !commOnly && ( reassignFreq>0 ) && ! (step%reassignFreq) ) {
       reassignVelocities(timestep,step);
    }
#endif

    doEnergy = ! ( step % energyFrequency );
#ifndef UPPER_BOUND
    if ( accelMDOn && !accelMDdihe ) doEnergy=1;
    //Update energy every timestep for adaptive tempering
    if ( adaptTempOn ) doEnergy=1;
#endif
D_MSG("runComputeObjects()");
    runComputeObjects(1,step<numberOfSteps); // must migrate here!
#ifndef UPPER_BOUND
    rescaleaccelMD(step, doNonbonded, doFullElectrostatics); // for accelMD 
    adaptTempUpdate(step); // update adaptive tempering temperature
#endif

#ifndef UPPER_BOUND
    if ( staleForces || doTcl || doColvars ) {
      if ( doNonbonded ) saveForce(Results::nbond);
      if ( doFullElectrostatics ) saveForce(Results::slow);
    }
    if ( ! commOnly ) {
D_MSG("newtonianVelocities()");
      TIMER_START(t, KICK);
      newtonianVelocities(-0.5,timestep,nbondstep,slowstep,0,1,1);
      TIMER_STOP(t, KICK);
    }
    minimizationQuenchVelocity();
#ifndef UPPER_BOUND
D_MSG("rattle1()");
    TIMER_START(t, RATTLE1);
    rattle1(-timestep,0);
    TIMER_STOP(t, RATTLE1);
#endif
D_MSG("submitHalfstep()");
    TIMER_START(t, SUBMITHALF);
    submitHalfstep(step);
    TIMER_STOP(t, SUBMITHALF);
    if ( ! commOnly ) {
D_MSG("newtonianVelocities()");
      TIMER_START(t, KICK);
      newtonianVelocities(1.0,timestep,nbondstep,slowstep,0,1,1);
      TIMER_STOP(t, KICK);
    }
D_MSG("rattle1()");
    TIMER_START(t, RATTLE1);
    rattle1(timestep,1);
    TIMER_STOP(t, RATTLE1);
    if (doTcl || doColvars)  // include constraint forces
      computeGlobal->saveTotalForces(patch);
D_MSG("submitHalfstep()");
    TIMER_START(t, SUBMITHALF);
    submitHalfstep(step);
    TIMER_STOP(t, SUBMITHALF);
    if ( zeroMomentum && doFullElectrostatics ) submitMomentum(step);
    if ( ! commOnly ) {
D_MSG("newtonianVelocities()");
      TIMER_START(t, KICK);
      newtonianVelocities(-0.5,timestep,nbondstep,slowstep,0,1,1);
      TIMER_STOP(t, KICK);
    }
#endif
D_MSG("submitReductions()");
    TIMER_START(t, SUBMITFULL);
    submitReductions(step);
    TIMER_STOP(t, SUBMITFULL);
#ifndef UPPER_BOUND
    if(0){ // if(traceIsOn()){
        traceUserEvent(eventEndOfTimeStep);
        sprintf(traceNote, "%s%d",tracePrefix,step); 
        traceUserSuppliedNote(traceNote);
    }
#endif
    rebalanceLoad(step);

  } // scriptTask == SCRIPT_RUN

#ifndef UPPER_BOUND
  bool doMultigratorRattle = false;
#endif

  //
  // DJH: There are a lot of mod operations below and elsewhere to
  // test step number against the frequency of something happening.
  // Mod and integer division are expensive!
  // Might be better to replace with counters and test equality.
  //
#if 0
    for(int i = 0; i < NamdProfileEvent::EventsCount; i++)
	CkPrintf("-------------- [%d] %s -------------\n", i, NamdProfileEventStr[i]);
#endif

#if defined(NAMD_NVTX_ENABLED) || defined(NAMD_CMK_TRACE_ENABLED)
  int& eon = patch->flags.event_on;
  int epid = (simParams->beginEventPatchID <= patch->getPatchID()
      && patch->getPatchID() <= simParams->endEventPatchID);
  int beginStep = simParams->beginEventStep;
  int endStep = simParams->endEventStep;
  bool controlProfiling = patch->getPatchID() == 0;
#endif
 
    for ( ++step; step <= numberOfSteps; ++step )
    {
#if defined(NAMD_NVTX_ENABLED) || defined(NAMD_CMK_TRACE_ENABLED)
      eon = epid && (beginStep < step && step <= endStep);

      if (controlProfiling && step == beginStep) {
        NAMD_PROFILE_START();
      }
      if (controlProfiling && step == endStep) {
        NAMD_PROFILE_STOP();
      }
      char buf[32];
      sprintf(buf, "%s: %d", NamdProfileEventStr[NamdProfileEvent::INTEGRATE_1], patch->getPatchID());
      NAMD_EVENT_START_EX(eon, NamdProfileEvent::INTEGRATE_1, buf);
#endif
#ifndef UPPER_BOUND

      rescaleVelocities(step);
      tcoupleVelocities(timestep,step);
      stochRescaleVelocities(timestep,step);
      berendsenPressure(step);

      if ( ! commOnly ) {
        TIMER_START(t, KICK);
        newtonianVelocities(0.5,timestep,nbondstep,slowstep,staleForces,doNonbonded,doFullElectrostatics); 
        TIMER_STOP(t, KICK);
      }

      // We do RATTLE here if multigrator thermostat was applied in the previous step
      if (doMultigratorRattle) rattle1(timestep, 1);

      /* reassignment based on half-step velocities
         if ( !commOnly && ( reassignFreq>0 ) && ! (step%reassignFreq) ) {
         addVelocityToPosition(0.5*timestep);
         reassignVelocities(timestep,step);
         addVelocityToPosition(0.5*timestep);
         rattle1(0.,0);
         rattle1(-timestep,0);
         addVelocityToPosition(-1.0*timestep);
         rattle1(timestep,0);
         } */

      TIMER_START(t, MAXMOVE);
      maximumMove(timestep);
      TIMER_STOP(t, MAXMOVE);

      NAMD_EVENT_STOP(eon, NamdProfileEvent::INTEGRATE_1);  // integrate 1

      if ( simParams->langevinPistonOn || (simParams->langevinOn && simParams->langevin_useBAOAB) ) {
        if ( ! commOnly ) {
          TIMER_START(t, DRIFT);
          addVelocityToPosition(0.5*timestep);
          TIMER_STOP(t, DRIFT);
        }
        // We add an Ornstein-Uhlenbeck integration step for the case of BAOAB (Langevin)
        langevinVelocities(timestep);

        // There is a blocking receive inside of langevinPiston()
        // that might suspend the current thread of execution,
        // so split profiling around this conditional block.
        langevinPiston(step);

        if ( ! commOnly ) {
          TIMER_START(t, DRIFT);
          addVelocityToPosition(0.5*timestep);
          TIMER_STOP(t, DRIFT);
        }
      } else {
        // If Langevin is not used, take full time step directly instread of two half steps      
        if ( ! commOnly ) {
          TIMER_START(t, DRIFT);
          addVelocityToPosition(timestep); 
          TIMER_STOP(t, DRIFT);
        }
      }

      NAMD_EVENT_START(eon, NamdProfileEvent::INTEGRATE_2);

      // impose hard wall potential for Drude bond length
      hardWallDrude(timestep, 1);

      minimizationQuenchVelocity();
#endif // UPPER_BOUND

      doNonbonded = !(step%nonbondedFrequency);
      doFullElectrostatics = (dofull && !(step%fullElectFrequency));

#ifndef UPPER_BOUND 
      if ( zeroMomentum && doFullElectrostatics ) {
        // There is a blocking receive inside of correctMomentum().
        correctMomentum(step,slowstep);
      }

      // There are NO sends in submitHalfstep() just local summation 
      // into the Reduction struct.
      TIMER_START(t, SUBMITHALF);
      submitHalfstep(step);
      TIMER_STOP(t, SUBMITHALF);

      doMolly = simParams->mollyOn && doFullElectrostatics;
      // BEGIN LA
      doLoweAndersen = simParams->loweAndersenOn && doNonbonded;
      // END LA

      maxForceUsed = Results::normal;
      if ( doNonbonded ) maxForceUsed = Results::nbond;
      if ( doFullElectrostatics ) maxForceUsed = Results::slow;
      if ( accelMDOn && (accelMDdihe || accelMDdual))  maxForceUsed = Results::amdf;

      // Migrate Atoms on stepsPerCycle
      doEnergy = ! ( step % energyFrequency );
      if ( accelMDOn && !accelMDdihe ) doEnergy=1;
      if ( adaptTempOn ) doEnergy=1; 

      // Multigrator 
      if (simParams->multigratorOn) {
        doVirial = (!(step % energyFrequency) || ((simParams->outputPressure > 0) && !(step % simParams->outputPressure)) 
          || !(step % simParams->multigratorPressureFreq));
        doKineticEnergy = (!(step % energyFrequency) || !(step % simParams->multigratorTemperatureFreq));
        doMomenta = (simParams->outputMomenta > 0) && !(step % simParams->outputMomenta);
      } else {
        doVirial = 1;
        doKineticEnergy = 1;
        doMomenta = 1;
      }
#endif
      NAMD_EVENT_STOP(eon, NamdProfileEvent::INTEGRATE_2);  // integrate 2

      // The current thread of execution will suspend in runComputeObjects().
      runComputeObjects(!(step%stepsPerCycle),step<numberOfSteps);

      NAMD_EVENT_START(eon, NamdProfileEvent::INTEGRATE_3);

#ifndef UPPER_BOUND
      rescaleaccelMD(step, doNonbonded, doFullElectrostatics); // for accelMD
     
      if ( staleForces || doTcl || doColvars ) {
        if ( doNonbonded ) saveForce(Results::nbond);
        if ( doFullElectrostatics ) saveForce(Results::slow);
      }

      // reassignment based on full-step velocities
      if ( !commOnly && ( reassignFreq>0 ) && ! (step%reassignFreq) ) {
        reassignVelocities(timestep,step);
        newtonianVelocities(-0.5,timestep,nbondstep,slowstep,staleForces,doNonbonded,doFullElectrostatics);
        rattle1(-timestep,0);
      }

      if ( ! commOnly ) {
        TIMER_START(t, VELBBK1);
        langevinVelocitiesBBK1(timestep);
        TIMER_STOP(t, VELBBK1);
        TIMER_START(t, KICK);
        newtonianVelocities(1.0,timestep,nbondstep,slowstep,staleForces,doNonbonded,doFullElectrostatics);
        TIMER_STOP(t, KICK);
        TIMER_START(t, VELBBK2);
        langevinVelocitiesBBK2(timestep);
        TIMER_STOP(t, VELBBK2);
      }

      // add drag to each atom's positions
      if ( ! commOnly && movDragOn ) addMovDragToPosition(timestep);
      if ( ! commOnly && rotDragOn ) addRotDragToPosition(timestep);

      TIMER_START(t, RATTLE1);
      rattle1(timestep,1);
      TIMER_STOP(t, RATTLE1);
      if (doTcl || doColvars)  // include constraint forces
        computeGlobal->saveTotalForces(patch);

      TIMER_START(t, SUBMITHALF);
      submitHalfstep(step);
      TIMER_STOP(t, SUBMITHALF);
      if ( zeroMomentum && doFullElectrostatics ) submitMomentum(step);

      if ( ! commOnly ) {
        TIMER_START(t, KICK);
        newtonianVelocities(-0.5,timestep,nbondstep,slowstep,staleForces,doNonbonded,doFullElectrostatics);
        TIMER_STOP(t, KICK);
      }

	// rattle2(timestep,step);
#endif

        TIMER_START(t, SUBMITFULL);
	submitReductions(step);
        TIMER_STOP(t, SUBMITFULL);
        TIMER_START(t, SUBMITCOLLECT);
	submitCollections(step);
        TIMER_STOP(t, SUBMITCOLLECT);
#ifndef UPPER_BOUND
       //Update adaptive tempering temperature
        adaptTempUpdate(step);

      // Multigrator temperature and pressure steps
      multigratorTemperature(step, 1);
      multigratorPressure(step, 1);
      multigratorPressure(step, 2);
      multigratorTemperature(step, 2);
      doMultigratorRattle = (simParams->multigratorOn && !(step % simParams->multigratorTemperatureFreq));

      NAMD_EVENT_STOP(eon, NamdProfileEvent::INTEGRATE_3); // integrate 3
#endif

#if CYCLE_BARRIER
        cycleBarrier(!((step+1) % stepsPerCycle), step);
#elif PME_BARRIER
        cycleBarrier(doFullElectrostatics, step);
#elif  STEP_BARRIER
        cycleBarrier(1, step);
#endif

#ifndef UPPER_BOUND
	 if(Node::Object()->specialTracing || simParams->statsOn){
		 int bstep = simParams->traceStartStep;
		 int estep = bstep + simParams->numTraceSteps;
		 if(step == bstep || step == estep){
			 traceBarrier(step);
		 }			 
	 }

#ifdef MEASURE_NAMD_WITH_PAPI	 
	 if(simParams->papiMeasure) {
		 int bstep = simParams->papiMeasureStartStep;
		 int estep = bstep + simParams->numPapiMeasureSteps;
		 if(step == bstep || step==estep) {
			 papiMeasureBarrier(step);
		 }
	 }
#endif
	  
        if(0){ // if(traceIsOn()){
            traceUserEvent(eventEndOfTimeStep);
            sprintf(traceNote, "%s%d",tracePrefix,step); 
            traceUserSuppliedNote(traceNote);
        }
#endif // UPPER_BOUND
	rebalanceLoad(step);

#if PME_BARRIER
	// a step before PME
        cycleBarrier(dofull && !((step+1)%fullElectFrequency),step);
#endif

#if USE_HPM
        if(step == START_HPM_STEP)
          (CProxy_Node(CkpvAccess(BOCclass_group).node)).startHPM();

        if(step == STOP_HPM_STEP)
          (CProxy_Node(CkpvAccess(BOCclass_group).node)).stopHPM();
#endif

    }

  TIMER_DONE(t);
#ifdef TIMER_COLLECTION
  if (patch->patchID == SPECIAL_PATCH_ID) {
    printf("Timer collection reporting in microseconds for "
        "Patch %d\n", patch->patchID);
    TIMER_REPORT(t);
  }
#endif // TIMER_COLLECTION
    //
    // DJH: Copy updates of SOA back into AOS.
    //
    //patch->copy_updates_to_AOS();
}

// add moving drag to each atom's position
void Sequencer::addMovDragToPosition(BigReal timestep) {
  FullAtom *atom = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Molecule *molecule = Node::Object()->molecule;   // need its methods
  const BigReal movDragGlobVel = simParams->movDragGlobVel;
  const BigReal dt = timestep / TIMEFACTOR;   // MUST be as in the integrator!
  Vector movDragVel, dragIncrement;
  for ( int i = 0; i < numAtoms; ++i )
  {
    // skip if fixed atom or zero drag attribute
    if ( (simParams->fixedAtomsOn && atom[i].atomFixed) 
	 || !(molecule->is_atom_movdragged(atom[i].id)) ) continue;
    molecule->get_movdrag_params(movDragVel, atom[i].id);
    dragIncrement = movDragGlobVel * movDragVel * dt;
    atom[i].position += dragIncrement;
  }
}

// add rotating drag to each atom's position
void Sequencer::addRotDragToPosition(BigReal timestep) {
  FullAtom *atom = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Molecule *molecule = Node::Object()->molecule;   // need its methods
  const BigReal rotDragGlobVel = simParams->rotDragGlobVel;
  const BigReal dt = timestep / TIMEFACTOR;   // MUST be as in the integrator!
  BigReal rotDragVel, dAngle;
  Vector atomRadius;
  Vector rotDragAxis, rotDragPivot, dragIncrement;
  for ( int i = 0; i < numAtoms; ++i )
  {
    // skip if fixed atom or zero drag attribute
    if ( (simParams->fixedAtomsOn && atom[i].atomFixed) 
	 || !(molecule->is_atom_rotdragged(atom[i].id)) ) continue;
    molecule->get_rotdrag_params(rotDragVel, rotDragAxis, rotDragPivot, atom[i].id);
    dAngle = rotDragGlobVel * rotDragVel * dt;
    rotDragAxis /= rotDragAxis.length();
    atomRadius = atom[i].position - rotDragPivot;
    dragIncrement = cross(rotDragAxis, atomRadius) * dAngle;
    atom[i].position += dragIncrement;
  }
}

void Sequencer::minimize() {
    //
    // DJH: Copy all data into SOA (structure of arrays)
    // from AOS (array of structures) data structure.
    //
    //patch->copy_all_to_SOA();

  const int numberOfSteps = simParams->N;
  const int stepsPerCycle = simParams->stepsPerCycle;
  int &step = patch->flags.step;
  step = simParams->firstTimestep;

  int &maxForceUsed = patch->flags.maxForceUsed;
  int &maxForceMerged = patch->flags.maxForceMerged;
  maxForceUsed = Results::normal;
  maxForceMerged = Results::normal;
  int &doNonbonded = patch->flags.doNonbonded;
  doNonbonded = 1;
  maxForceUsed = Results::nbond;
  maxForceMerged = Results::nbond;
  const int dofull = ( simParams->fullElectFrequency ? 1 : 0 );
  int &doFullElectrostatics = patch->flags.doFullElectrostatics;
  doFullElectrostatics = dofull;
  if ( dofull ) {
    maxForceMerged = Results::slow;
    maxForceUsed = Results::slow;
  }
  int &doMolly = patch->flags.doMolly;
  doMolly = simParams->mollyOn && doFullElectrostatics;
  // BEGIN LA
  int &doLoweAndersen = patch->flags.doLoweAndersen;
  doLoweAndersen = 0;
  // END LA

  int &doGBIS = patch->flags.doGBIS;
  doGBIS = simParams->GBISOn;

    int &doLCPO = patch->flags.doLCPO;
    doLCPO = simParams->LCPOOn;

    int doTcl = simParams->tclForcesOn;
	int doColvars = simParams->colvarsOn;
    ComputeGlobal *computeGlobal = Node::Object()->computeMgr->computeGlobalObject;

  int &doEnergy = patch->flags.doEnergy;
  doEnergy = 1;

  patch->rattle1(0.,0,0);  // enforce rigid bond constraints on initial positions

  if (simParams->lonepairs || simParams->singleTopology) {
    patch->atomMapper->registerIDsFullAtom(
		patch->atom.begin(),patch->atom.end());
  }

  runComputeObjects(1,step<numberOfSteps); // must migrate here!

  if ( doTcl || doColvars ) {
    if ( doNonbonded ) saveForce(Results::nbond);
    if ( doFullElectrostatics ) saveForce(Results::slow);
    computeGlobal->saveTotalForces(patch);
  }

  BigReal fmax2 = TIMEFACTOR * TIMEFACTOR * TIMEFACTOR * TIMEFACTOR;

  submitMinimizeReductions(step,fmax2);
  rebalanceLoad(step);

  int downhill = 1;  // start out just fixing bad contacts
  int minSeq = 0;
  for ( ++step; step <= numberOfSteps; ++step ) {

   // Blocking receive for the minimization coefficient.
   BigReal c = broadcast->minimizeCoefficient.get(minSeq++);

   if ( downhill ) {
    if ( c ) minimizeMoveDownhill(fmax2);
    else {
      downhill = 0;
      fmax2 *= 10000.;
    }
   }
   if ( ! downhill ) {
    if ( ! c ) {  // new direction

      // Blocking receive for the minimization coefficient.
      c = broadcast->minimizeCoefficient.get(minSeq++);

      newMinimizeDirection(c);  // v = c * v + f

      // Blocking receive for the minimization coefficient.
      c = broadcast->minimizeCoefficient.get(minSeq++);

    }  // same direction
    newMinimizePosition(c);  // x = x + c * v
   }

    runComputeObjects(!(step%stepsPerCycle),step<numberOfSteps);
    if ( doTcl || doColvars ) {
      if ( doNonbonded ) saveForce(Results::nbond);
      if ( doFullElectrostatics ) saveForce(Results::slow);
      computeGlobal->saveTotalForces(patch);
    }
    submitMinimizeReductions(step,fmax2);
    submitCollections(step, 1);  // write out zeros for velocities
    rebalanceLoad(step);
  }
  quenchVelocities();  // zero out bogus velocity

    //
    // DJH: Copy updates of SOA back into AOS.
    //
    //patch->copy_updates_to_AOS();
}

// x = x + 0.1 * unit(f) for large f
void Sequencer::minimizeMoveDownhill(BigReal fmax2) {

  FullAtom *a = patch->atom.begin();
  Force *f1 = patch->f[Results::normal].begin();  // includes nbond and slow
  int numAtoms = patch->numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    if ( simParams->fixedAtomsOn && a[i].atomFixed ) continue;
    Force f = f1[i];
    if ( f.length2() > fmax2 ) {
      a[i].position += ( 0.1 * f.unit() );
      int hgs = a[i].hydrogenGroupSize;  // 0 if not parent
      for ( int j=1; j<hgs; ++j ) {
        a[++i].position += ( 0.1 * f.unit() );
      }
    }
  }

  patch->rattle1(0.,0,0);
}

// v = c * v + f
void Sequencer::newMinimizeDirection(BigReal c) {
  FullAtom *a = patch->atom.begin();
  Force *f1 = patch->f[Results::normal].begin(); // includes nbond and slow
  const bool fixedAtomsOn = simParams->fixedAtomsOn;
  const bool drudeHardWallOn = simParams->drudeHardWallOn;
  int numAtoms = patch->numAtoms;
  BigReal maxv2 = 0.;

  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].velocity *= c;
    a[i].velocity += f1[i];
    if ( drudeHardWallOn && i && (0.05 < a[i].mass) && ((a[i].mass < 1.0)) ) { // drude particle
      a[i].velocity = a[i-1].velocity;
    }
    if ( fixedAtomsOn && a[i].atomFixed ) a[i].velocity = 0;
    BigReal v2 = a[i].velocity.length2();
    if ( v2 > maxv2 ) maxv2 = v2;
  }

  { Tensor virial; patch->minimize_rattle2( 0.1 * TIMEFACTOR / sqrt(maxv2), &virial); }

  maxv2 = 0.;
  for ( int i = 0; i < numAtoms; ++i ) {
    if ( drudeHardWallOn && i && (0.05 < a[i].mass) && ((a[i].mass < 1.0)) ) { // drude particle
      a[i].velocity = a[i-1].velocity;
    }
    if ( fixedAtomsOn && a[i].atomFixed ) a[i].velocity = 0;
    BigReal v2 = a[i].velocity.length2();
    if ( v2 > maxv2 ) maxv2 = v2;
  }

  min_reduction->max(0,maxv2);
  min_reduction->submit();

  // prevent hydrogens from being left behind
  BigReal fmax2 = 0.01 * TIMEFACTOR * TIMEFACTOR * TIMEFACTOR * TIMEFACTOR;
  // int adjustCount = 0;
  int hgs;
  for ( int i = 0; i < numAtoms; i += hgs ) {
    hgs = a[i].hydrogenGroupSize;
    BigReal minChildVel = a[i].velocity.length2();
    if ( minChildVel < fmax2 ) continue;
    int adjustChildren = 1;
    for ( int j = i+1; j < (i+hgs); ++j ) {
      if ( a[j].velocity.length2() > minChildVel ) adjustChildren = 0;
    }
    if ( adjustChildren ) {
      // if ( hgs > 1 ) ++adjustCount;
      for ( int j = i+1; j < (i+hgs); ++j ) {
        if (a[i].mass < 0.01) continue;  // lone pair
        a[j].velocity = a[i].velocity;
      }
    }
  }
  // if (adjustCount) CkPrintf("Adjusting %d hydrogen groups\n", adjustCount);

}

// x = x + c * v
void Sequencer::newMinimizePosition(BigReal c) {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].position += c * a[i].velocity;
  }

  if ( simParams->drudeHardWallOn ) {
    for ( int i = 1; i < numAtoms; ++i ) {
      if ( (0.05 < a[i].mass) && ((a[i].mass < 1.0)) ) { // drude particle
        a[i].position -= a[i-1].position;
      }
    }
  }

  patch->rattle1(0.,0,0);

  if ( simParams->drudeHardWallOn ) {
    for ( int i = 1; i < numAtoms; ++i ) {
      if ( (0.05 < a[i].mass) && ((a[i].mass < 1.0)) ) { // drude particle
        a[i].position += a[i-1].position;
      }
    }
  }
}

// v = 0
void Sequencer::quenchVelocities() {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].velocity = 0;
  }
}

void Sequencer::submitMomentum(int step) {

  FullAtom *a = patch->atom.begin();
  const int numAtoms = patch->numAtoms;

  Vector momentum = 0;
  BigReal mass = 0;
if ( simParams->zeroMomentumAlt ) {
  for ( int i = 0; i < numAtoms; ++i ) {
    momentum += a[i].mass * a[i].velocity;
    mass += 1.;
  }
} else {
  for ( int i = 0; i < numAtoms; ++i ) {
    momentum += a[i].mass * a[i].velocity;
    mass += a[i].mass;
  }
}

  ADD_VECTOR_OBJECT(reduction,REDUCTION_HALFSTEP_MOMENTUM,momentum);
  reduction->item(REDUCTION_MOMENTUM_MASS) += mass;
}

void Sequencer::correctMomentum(int step, BigReal drifttime) {

  //
  // DJH: This test should be done in SimParameters.
  //
  if ( simParams->fixedAtomsOn )
    NAMD_die("Cannot zero momentum when fixed atoms are present.");

  // Blocking receive for the momentum correction vector.
  const Vector dv = broadcast->momentumCorrection.get(step);

  const Vector dx = dv * ( drifttime / TIMEFACTOR );

  FullAtom *a = patch->atom.begin();
  const int numAtoms = patch->numAtoms;

if ( simParams->zeroMomentumAlt ) {
  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].velocity += dv * a[i].recipMass;
    a[i].position += dx * a[i].recipMass;
  }
} else {
  for ( int i = 0; i < numAtoms; ++i ) {
    a[i].velocity += dv;
    a[i].position += dx;
  }
}

}

// --------- For Multigrator ---------
void Sequencer::scalePositionsVelocities(const Tensor& posScale, const Tensor& velScale) {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Position origin = patch->lattice.origin();
  if ( simParams->fixedAtomsOn ) {
    NAMD_bug("Sequencer::scalePositionsVelocities, fixed atoms not implemented");
  }
  if ( simParams->useGroupPressure ) {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      Position pos_cm(0.0, 0.0, 0.0);
      Velocity vel_cm(0.0, 0.0, 0.0);
      BigReal m_cm = 0.0;
      for (int j=0;j < hgs;++j) {
        m_cm += a[i+j].mass;
        pos_cm += a[i+j].mass*a[i+j].position;
        vel_cm += a[i+j].mass*a[i+j].velocity;
      }
      pos_cm /= m_cm;
      vel_cm /= m_cm;
      pos_cm -= origin;
      Position dpos = posScale*pos_cm;
      Velocity dvel = velScale*vel_cm;
      for (int j=0;j < hgs;++j) {
        a[i+j].position += dpos;
        a[i+j].velocity += dvel;
      }
    }
  } else {
    for ( int i = 0; i < numAtoms; i++) {
      a[i].position += posScale*(a[i].position-origin);
      a[i].velocity = velScale*a[i].velocity;
    }
  }
}

void Sequencer::multigratorPressure(int step, int callNumber) {
// Calculate new positions, momenta, and volume using positionRescaleFactor and 
// velocityRescaleTensor values returned from Controller::multigratorPressureCalcScale()
  if (simParams->multigratorOn && !(step % simParams->multigratorPressureFreq)) {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;

    // Blocking receive (get) scaling factors from Controller
    Tensor scaleTensor    = (callNumber == 1) ? broadcast->positionRescaleFactor.get(step) : broadcast->positionRescaleFactor2.get(step);
    Tensor velScaleTensor = (callNumber == 1) ? broadcast->velocityRescaleTensor.get(step) : broadcast->velocityRescaleTensor2.get(step);
    Tensor posScaleTensor = scaleTensor;
    posScaleTensor -= Tensor::identity();
    if (simParams->useGroupPressure) {
      velScaleTensor -= Tensor::identity();
    }

    // Scale volume
    patch->lattice.rescale(scaleTensor);
    // Scale positions and velocities
    scalePositionsVelocities(posScaleTensor, velScaleTensor);

    if (!patch->flags.doFullElectrostatics) NAMD_bug("Sequencer::multigratorPressure, doFullElectrostatics must be true");

    // Calculate new forces
    // NOTE: We should not need to migrate here since any migration should have happened in the
    // previous call to runComputeObjects inside the MD loop in Sequencer::integrate()
    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;
    runComputeObjects(0 /*!(step%stepsPerCycle)*/, step<numberOfSteps, 1);

    reduction->item(REDUCTION_ATOM_CHECKSUM) += numAtoms;
    reduction->item(REDUCTION_MARGIN_VIOLATIONS) += patch->marginViolations;

    // Virials etc.
    Tensor virialNormal;
    Tensor momentumSqrSum;
    BigReal kineticEnergy = 0;
    if ( simParams->pairInteractionOn ) {
      if ( simParams->pairInteractionSelf ) {
        for ( int i = 0; i < numAtoms; ++i ) {
          if ( a[i].partition != 1 ) continue;
          kineticEnergy += a[i].mass * a[i].velocity.length2();
          virialNormal.outerAdd(a[i].mass, a[i].velocity, a[i].velocity);
        }
      }
    } else {
      for ( int i = 0; i < numAtoms; ++i ) {
        if (a[i].mass < 0.01) continue;
        kineticEnergy += a[i].mass * a[i].velocity.length2();
        virialNormal.outerAdd(a[i].mass, a[i].velocity, a[i].velocity);
      }
    }
    if (!simParams->useGroupPressure) momentumSqrSum = virialNormal;
    kineticEnergy *= 0.5;
    reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY) += kineticEnergy;
    ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NORMAL, virialNormal);

    if ( simParams->fixedAtomsOn ) {
      Tensor fixVirialNormal;
      Tensor fixVirialNbond;
      Tensor fixVirialSlow;
      Vector fixForceNormal = 0;
      Vector fixForceNbond = 0;
      Vector fixForceSlow = 0;

      calcFixVirial(fixVirialNormal, fixVirialNbond, fixVirialSlow, fixForceNormal, fixForceNbond, fixForceSlow);

      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NORMAL, fixVirialNormal);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, fixVirialNbond);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, fixVirialSlow);
      ADD_VECTOR_OBJECT(reduction, REDUCTION_EXT_FORCE_NORMAL, fixForceNormal);
      ADD_VECTOR_OBJECT(reduction, REDUCTION_EXT_FORCE_NBOND, fixForceNbond);
      ADD_VECTOR_OBJECT(reduction, REDUCTION_EXT_FORCE_SLOW, fixForceSlow);
    }

    // Internal virial and group momentum
    Tensor intVirialNormal;
    Tensor intVirialNormal2;
    Tensor intVirialNbond;
    Tensor intVirialSlow;
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
        v_cm += a[j].mass * a[j].velocity;
      }
      if (simParams->useGroupPressure) momentumSqrSum.outerAdd(1.0/m_cm, v_cm, v_cm);
      x_cm /= m_cm;
      v_cm /= m_cm;
      if (simParams->fixedAtomsOn) NAMD_bug("Sequencer::multigratorPressure, simParams->fixedAtomsOn not implemented yet");
      if ( simParams->pairInteractionOn ) {
        if ( simParams->pairInteractionSelf ) {
          NAMD_bug("Sequencer::multigratorPressure, this part needs to be implemented correctly");
          for ( j = i; j < (i+hgs); ++j ) {
            if ( a[j].partition != 1 ) continue;
            BigReal mass = a[j].mass;
            Vector v = a[j].velocity;
            Vector dv = v - v_cm;
            intVirialNormal2.outerAdd (mass, v, dv);
            Vector dx = a[j].position - x_cm;
            intVirialNormal.outerAdd(1.0, patch->f[Results::normal][j], dx);
            intVirialNbond.outerAdd(1.0, patch->f[Results::nbond][j], dx);
            intVirialSlow.outerAdd(1.0, patch->f[Results::slow][j], dx);
          }
        }
      } else {
        for ( j = i; j < (i+hgs); ++j ) {
          BigReal mass = a[j].mass;
          Vector v = a[j].velocity;
          Vector dv = v - v_cm;
          intVirialNormal2.outerAdd(mass, v, dv);
          Vector dx = a[j].position - x_cm;
          intVirialNormal.outerAdd(1.0, patch->f[Results::normal][j], dx);
          intVirialNbond.outerAdd(1.0, patch->f[Results::nbond][j], dx);
          intVirialSlow.outerAdd(1.0, patch->f[Results::slow][j], dx);
        }
      }
    }

    ADD_TENSOR_OBJECT(reduction, REDUCTION_INT_VIRIAL_NORMAL, intVirialNormal);
    ADD_TENSOR_OBJECT(reduction, REDUCTION_INT_VIRIAL_NORMAL, intVirialNormal2);
    ADD_TENSOR_OBJECT(reduction, REDUCTION_INT_VIRIAL_NBOND, intVirialNbond);
    ADD_TENSOR_OBJECT(reduction, REDUCTION_INT_VIRIAL_SLOW, intVirialSlow);
    ADD_TENSOR_OBJECT(reduction, REDUCTION_MOMENTUM_SQUARED, momentumSqrSum);

    reduction->submit();
  }
}

void Sequencer::scaleVelocities(const BigReal velScale) {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  for ( int i = 0; i < numAtoms; i++) {
    a[i].velocity *= velScale;
  }
}

BigReal Sequencer::calcKineticEnergy() {
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  BigReal kineticEnergy = 0.0;
  if ( simParams->pairInteractionOn ) {
    if ( simParams->pairInteractionSelf ) {
      for (int i = 0; i < numAtoms; ++i ) {
        if ( a[i].partition != 1 ) continue;
        kineticEnergy += a[i].mass * a[i].velocity.length2();
      }
    }
  } else {
    for (int i = 0; i < numAtoms; ++i ) {
      kineticEnergy += a[i].mass * a[i].velocity.length2();
    }
  }
  kineticEnergy *= 0.5;
  return kineticEnergy;
}

void Sequencer::multigratorTemperature(int step, int callNumber) {
  if (simParams->multigratorOn && !(step % simParams->multigratorTemperatureFreq)) {
    // Blocking receive (get) velocity scaling factor.
    BigReal velScale = (callNumber == 1) ? broadcast->velocityRescaleFactor.get(step) : broadcast->velocityRescaleFactor2.get(step);
    scaleVelocities(velScale);
    // Calculate new kineticEnergy
    BigReal kineticEnergy = calcKineticEnergy();
    multigratorReduction->item(MULTIGRATOR_REDUCTION_KINETIC_ENERGY) += kineticEnergy;
    if (callNumber == 1 && !(step % simParams->multigratorPressureFreq)) {
      // If this is a pressure cycle, calculate new momentum squared sum
      FullAtom *a = patch->atom.begin();
      int numAtoms = patch->numAtoms;
      Tensor momentumSqrSum;
      if (simParams->useGroupPressure) {
        int hgs;
        for ( int i = 0; i < numAtoms; i += hgs ) {
          hgs = a[i].hydrogenGroupSize;
          int j;
          BigReal m_cm = 0;
          Position x_cm(0,0,0);
          Velocity v_cm(0,0,0);
          for ( j = i; j < (i+hgs); ++j ) {
            m_cm += a[j].mass;
            x_cm += a[j].mass * a[j].position;
            v_cm += a[j].mass * a[j].velocity;
          }
          momentumSqrSum.outerAdd(1.0/m_cm, v_cm, v_cm);
        }
      } else {
        for ( int i = 0; i < numAtoms; i++) {
          momentumSqrSum.outerAdd(a[i].mass, a[i].velocity, a[i].velocity);
        }
      }
      ADD_TENSOR_OBJECT(multigratorReduction, MULTIGRATOR_REDUCTION_MOMENTUM_SQUARED, momentumSqrSum);
    }
    // Submit reductions (kineticEnergy and, if applicable, momentumSqrSum)
    multigratorReduction->submit();

  }
}
// --------- End Multigrator ---------

//
// DJH: Calls one or more addForceToMomentum which in turn calls HomePatch
// versions. We should inline to reduce the number of function calls.
//
void Sequencer::newtonianVelocities(BigReal stepscale, const BigReal timestep, 
                                    const BigReal nbondstep, 
                                    const BigReal slowstep, 
                                    const int staleForces, 
                                    const int doNonbonded,
                                    const int doFullElectrostatics)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::NEWTONIAN_VELOCITIES);

  // Deterministic velocity update, account for multigrator
  if (staleForces || (doNonbonded && doFullElectrostatics)) {
    addForceToMomentum3(stepscale*timestep, Results::normal, 0,
                        stepscale*nbondstep, Results::nbond, staleForces,
                        stepscale*slowstep, Results::slow, staleForces);
  } else {
    addForceToMomentum(stepscale*timestep);
    if (staleForces || doNonbonded)
      addForceToMomentum(stepscale*nbondstep, Results::nbond, staleForces);
    if (staleForces || doFullElectrostatics)
      addForceToMomentum(stepscale*slowstep, Results::slow, staleForces);
  }
}

void Sequencer::langevinVelocities(BigReal dt_fs)
{
// This routine is used for the BAOAB integrator,
// Ornstein-Uhlenbeck exact solve for the O-part.
// See B. Leimkuhler and C. Matthews, AMRX (2012)
// Routine originally written by JPhillips, with fresh errors by CMatthews June2012

  if ( simParams->langevinOn && simParams->langevin_useBAOAB )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMANN*(simParams->langevinTemp);
    if (simParams->adaptTempOn && simParams->adaptTempLangevin)
    {
        kbT = BOLTZMANN*adaptTempT;
    }

    int lesReduceTemp = simParams->lesOn && simParams->lesReduceTemp;
    BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

    for ( int i = 0; i < numAtoms; ++i )
    {
      BigReal dt_gamma = dt * a[i].langevinParam;
      if ( ! dt_gamma ) continue;

      BigReal f1 = exp( -dt_gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kbT * 
                         ( a[i].partition ? tempFactor : 1.0 ) * 
                         a[i].recipMass );
      a[i].velocity *= f1;
      a[i].velocity += f2 * random->gaussian_vector();
    }
  }
}

void Sequencer::langevinVelocitiesBBK1(BigReal dt_fs)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::LANGEVIN_VELOCITIES_BBK1);
  if ( simParams->langevinOn && !simParams->langevin_useBAOAB )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    int i;

    if (simParams->drudeOn) {
      for (i = 0;  i < numAtoms;  i++) {

        if (i < numAtoms-1 &&
            a[i+1].mass < 1.0 && a[i+1].mass > 0.05) {
          //printf("*** Found Drude particle %d\n", a[i+1].id);
          // i+1 is a Drude particle with parent i

          // convert from Cartesian coordinates to (COM,bond) coordinates
          BigReal m = a[i+1].mass / (a[i].mass + a[i+1].mass);  // mass ratio
          Vector v_bnd = a[i+1].velocity - a[i].velocity;  // vel of bond
          Vector v_com = a[i].velocity + m * v_bnd;  // vel of COM
          BigReal dt_gamma;

          // use Langevin damping factor i for v_com
          dt_gamma = dt * a[i].langevinParam;
          if (dt_gamma != 0.0) {
            v_com *= ( 1. - 0.5 * dt_gamma );
          }

          // use Langevin damping factor i+1 for v_bnd
          dt_gamma = dt * a[i+1].langevinParam;
          if (dt_gamma != 0.0) {
            v_bnd *= ( 1. - 0.5 * dt_gamma );
          }

          // convert back
          a[i].velocity = v_com - m * v_bnd;
          a[i+1].velocity = v_bnd + a[i].velocity;

          i++;  // +1 from loop, we've updated both particles
        }
        else {
          BigReal dt_gamma = dt * a[i].langevinParam;
          if ( ! dt_gamma ) continue;

          a[i].velocity *= ( 1. - 0.5 * dt_gamma );
        }

      } // end for
    } // end if drudeOn
    else {

      //
      // DJH: The conditional inside loop prevents vectorization and doesn't
      // avoid much work since addition and multiplication are cheap.
      //
      for ( i = 0; i < numAtoms; ++i )
      {
        BigReal dt_gamma = dt * a[i].langevinParam;
        if ( ! dt_gamma ) continue;

        a[i].velocity *= ( 1. - 0.5 * dt_gamma );
      }

    } // end else

  } // end if langevinOn
}


void Sequencer::langevinVelocitiesBBK2(BigReal dt_fs)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::LANGEVIN_VELOCITIES_BBK2);
  if ( simParams->langevinOn && !simParams->langevin_useBAOAB ) 
  {
    //
    // DJH: This call is expensive. Avoid calling when gammas don't differ.
    // Set flag in SimParameters and make this call conditional.
    // 
    TIMER_START(patch->timerSet, RATTLE1);
    rattle1(dt_fs,1);  // conserve momentum if gammas differ
    TIMER_STOP(patch->timerSet, RATTLE1);

    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    BigReal kbT = BOLTZMANN*(simParams->langevinTemp);
    if (simParams->adaptTempOn && simParams->adaptTempLangevin)
    {
        kbT = BOLTZMANN*adaptTempT;
    }
    int lesReduceTemp = simParams->lesOn && simParams->lesReduceTemp;
    BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;
    int i;

    if (simParams->drudeOn) {
      BigReal kbT_bnd = BOLTZMANN*(simParams->drudeTemp);  // drude bond Temp

      for (i = 0;  i < numAtoms;  i++) {

        if (i < numAtoms-1 &&
            a[i+1].mass < 1.0 && a[i+1].mass > 0.05) {
          //printf("*** Found Drude particle %d\n", a[i+1].id);
          // i+1 is a Drude particle with parent i

          // convert from Cartesian coordinates to (COM,bond) coordinates
          BigReal m = a[i+1].mass / (a[i].mass + a[i+1].mass);  // mass ratio
          Vector v_bnd = a[i+1].velocity - a[i].velocity;  // vel of bond
          Vector v_com = a[i].velocity + m * v_bnd;  // vel of COM
          BigReal dt_gamma;

          // use Langevin damping factor i for v_com
          dt_gamma = dt * a[i].langevinParam;
          if (dt_gamma != 0.0) {
            BigReal mass = a[i].mass + a[i+1].mass;
            v_com += random->gaussian_vector() *
              sqrt( 2 * dt_gamma * kbT *
                  ( a[i].partition ? tempFactor : 1.0 ) / mass );
            v_com /= ( 1. + 0.5 * dt_gamma );
          }

          // use Langevin damping factor i+1 for v_bnd
          dt_gamma = dt * a[i+1].langevinParam;
          if (dt_gamma != 0.0) {
            BigReal mass = a[i+1].mass * (1. - m);
            v_bnd += random->gaussian_vector() *
              sqrt( 2 * dt_gamma * kbT_bnd *
                  ( a[i+1].partition ? tempFactor : 1.0 ) / mass );
            v_bnd /= ( 1. + 0.5 * dt_gamma );
          }

          // convert back
          a[i].velocity = v_com - m * v_bnd;
          a[i+1].velocity = v_bnd + a[i].velocity;

          i++;  // +1 from loop, we've updated both particles
        }
        else { 
          BigReal dt_gamma = dt * a[i].langevinParam;
          if ( ! dt_gamma ) continue;

          a[i].velocity += random->gaussian_vector() *
            sqrt( 2 * dt_gamma * kbT *
                ( a[i].partition ? tempFactor : 1.0 ) * a[i].recipMass );
          a[i].velocity /= ( 1. + 0.5 * dt_gamma );
        }

      } // end for
    } // end if drudeOn
    else {

      //
      // DJH: For case using same gamma (the Langevin parameter),
      // no partitions (e.g. FEP), and no adaptive tempering (adaptTempMD),
      // we can precompute constants. Then by lifting the RNG from the
      // loop (filling up an array of random numbers), we can vectorize
      // loop and simplify arithmetic to just addition and multiplication.
      //
      for ( i = 0; i < numAtoms; ++i )
      {
        BigReal dt_gamma = dt * a[i].langevinParam;
        if ( ! dt_gamma ) continue;

        a[i].velocity += random->gaussian_vector() *
          sqrt( 2 * dt_gamma * kbT *
              ( a[i].partition ? tempFactor : 1.0 ) * a[i].recipMass );
        a[i].velocity /= ( 1. + 0.5 * dt_gamma );
      }

    } // end else

  } // end if langevinOn
}


void Sequencer::berendsenPressure(int step)
{
  if ( simParams->berendsenPressureOn ) {
  berendsenPressure_count += 1;
  const int freq = simParams->berendsenPressureFreq;
  if ( ! (berendsenPressure_count % freq ) ) {
   berendsenPressure_count = 0;
   FullAtom *a = patch->atom.begin();
   int numAtoms = patch->numAtoms;
   // Blocking receive for the updated lattice scaling factor.
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   patch->lattice.rescale(factor);
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      int j;
      hgs = a[i].hydrogenGroupSize;
      if ( simParams->fixedAtomsOn && a[i].groupFixed ) {
        for ( j = i; j < (i+hgs); ++j ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
        }
        continue;
      }
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
      }
      x_cm /= m_cm;
      Position new_x_cm = x_cm;
      patch->lattice.rescale(new_x_cm,factor);
      Position delta_x_cm = new_x_cm - x_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
          continue;
        }
        a[j].position += delta_x_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < numAtoms; ++i )
    {
      if ( simParams->fixedAtomsOn && a[i].atomFixed ) {
        a[i].position = patch->lattice.apply_transform(
				a[i].fixedPosition,a[i].transform);
        continue;
      }
      patch->lattice.rescale(a[i].position,factor);
    }
   }
  }
  } else {
    berendsenPressure_count = 0;
  }
}

void Sequencer::langevinPiston(int step)
{
  if ( simParams->langevinPistonOn && ! ( (step-1-slowFreq/2) % slowFreq ) )
  {
    //
    // DJH: Loops below simplify if we lift out special cases of fixed atoms
    // and pressure excluded atoms and make them their own branch.
    //
   FullAtom *a = patch->atom.begin();
   int numAtoms = patch->numAtoms;
   // Blocking receive for the updated lattice scaling factor.
   Tensor factor = broadcast->positionRescaleFactor.get(step);
   TIMER_START(patch->timerSet, PISTON);
   // JCP FIX THIS!!!
   Vector velFactor(1/factor.xx,1/factor.yy,1/factor.zz);
   patch->lattice.rescale(factor);
   Molecule *mol = Node::Object()->molecule;
   if ( simParams->useGroupPressure )
   {
    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      int j;
      hgs = a[i].hydrogenGroupSize;
      if ( simParams->fixedAtomsOn && a[i].groupFixed ) {
        for ( j = i; j < (i+hgs); ++j ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
        }
        continue;
      }
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
        v_cm += a[j].mass * a[j].velocity;
      }
      x_cm /= m_cm;
      Position new_x_cm = x_cm;
      patch->lattice.rescale(new_x_cm,factor);
      Position delta_x_cm = new_x_cm - x_cm;
      v_cm /= m_cm;
      Velocity delta_v_cm;
      delta_v_cm.x = ( velFactor.x - 1 ) * v_cm.x;
      delta_v_cm.y = ( velFactor.y - 1 ) * v_cm.y;
      delta_v_cm.z = ( velFactor.z - 1 ) * v_cm.z;
      for ( j = i; j < (i+hgs); ++j ) {
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
          a[j].position = patch->lattice.apply_transform(
				a[j].fixedPosition,a[j].transform);
          continue;
        }
        if ( mol->is_atom_exPressure(a[j].id) ) continue;
        a[j].position += delta_x_cm;
        a[j].velocity += delta_v_cm;
      }
    }
   }
   else
   {
    for ( int i = 0; i < numAtoms; ++i )
    {
      if ( simParams->fixedAtomsOn && a[i].atomFixed ) {
        a[i].position = patch->lattice.apply_transform(
				a[i].fixedPosition,a[i].transform);
        continue;
      }
      if ( mol->is_atom_exPressure(a[i].id) ) continue;
      patch->lattice.rescale(a[i].position,factor);
      a[i].velocity.x *= velFactor.x;
      a[i].velocity.y *= velFactor.y;
      a[i].velocity.z *= velFactor.z;
    }
   }
   TIMER_STOP(patch->timerSet, PISTON);
  }
}

void Sequencer::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    ++rescaleVelocities_numTemps;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      // Blocking receive for the velcity scaling factor.
      BigReal factor = broadcast->velocityRescaleFactor.get(step);
      for ( int i = 0; i < numAtoms; ++i )
      {
        a[i].velocity *= factor;
      }
      rescaleVelocities_numTemps = 0;
    }
  }
}

void Sequencer::rescaleaccelMD (int step, int doNonbonded, int doFullElectrostatics)
{
   if (!simParams->accelMDOn) return;
   if ((step < simParams->accelMDFirstStep) || ( simParams->accelMDLastStep >0 && step > simParams->accelMDLastStep)) return;

   // Blocking receive for the Accelerated MD scaling factors.
   Vector accelMDfactor = broadcast->accelMDRescaleFactor.get(step);
   const BigReal factor_dihe = accelMDfactor[0];
   const BigReal factor_tot  = accelMDfactor[1];
   const int numAtoms = patch->numAtoms;

   if (simParams->accelMDdihe && factor_tot <1 )
       NAMD_die("accelMD broadcasting error!\n");
   if (!simParams->accelMDdihe && !simParams->accelMDdual && factor_dihe <1 )
       NAMD_die("accelMD broadcasting error!\n");

   if (simParams->accelMDdihe && factor_dihe < 1) {
        for (int i = 0; i < numAtoms; ++i)
          if (patch->f[Results::amdf][i][0] || patch->f[Results::amdf][i][1] || patch->f[Results::amdf][i][2])
              patch->f[Results::normal][i] += patch->f[Results::amdf][i]*(factor_dihe - 1);         
   }

   if ( !simParams->accelMDdihe && factor_tot < 1) {
        for (int i = 0; i < numAtoms; ++i)
            patch->f[Results::normal][i] *= factor_tot;
        if (doNonbonded) {
            for (int i = 0; i < numAtoms; ++i)
                 patch->f[Results::nbond][i] *= factor_tot;
        }
        if (doFullElectrostatics) {
            for (int i = 0; i < numAtoms; ++i)
                 patch->f[Results::slow][i] *= factor_tot;
        }
   }

   if (simParams->accelMDdual && factor_dihe < 1) {
        for (int i = 0; i < numAtoms; ++i)
           if (patch->f[Results::amdf][i][0] || patch->f[Results::amdf][i][1] || patch->f[Results::amdf][i][2])
               patch->f[Results::normal][i] += patch->f[Results::amdf][i]*(factor_dihe - factor_tot);
   }

}

void Sequencer::adaptTempUpdate(int step)
{
   //check if adaptive tempering is enabled and in the right timestep range
   if (!simParams->adaptTempOn) return;
   if ( (step < simParams->adaptTempFirstStep ) || 
     ( simParams->adaptTempLastStep > 0 && step > simParams->adaptTempLastStep )) {
        if (simParams->langevinOn) // restore langevin temperature
            adaptTempT = simParams->langevinTemp;
        return;
   }
   // Get Updated Temperature
   if ( !(step % simParams->adaptTempFreq ) && (step > simParams->firstTimestep ))
     // Blocking receive for the updated adaptive tempering temperature.
     adaptTempT = broadcast->adaptTemperature.get(step);
}

void Sequencer::reassignVelocities(BigReal timestep, int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    BigReal newTemp = simParams->reassignTemp;
    newTemp += ( step / reassignFreq ) * simParams->reassignIncr;
    if ( simParams->reassignIncr > 0.0 ) {
      if ( newTemp > simParams->reassignHold && simParams->reassignHold > 0.0 )
        newTemp = simParams->reassignHold;
    } else {
      if ( newTemp < simParams->reassignHold )
        newTemp = simParams->reassignHold;
    }
    BigReal kbT = BOLTZMANN * newTemp;

    int lesReduceTemp = simParams->lesOn && simParams->lesReduceTemp;
    BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

    for ( int i = 0; i < numAtoms; ++i )
    {
      a[i].velocity = ( ( simParams->fixedAtomsOn && 
            a[i].atomFixed && a[i].mass > 0.) ? Vector(0,0,0) :
          sqrt(kbT * (a[i].partition ? tempFactor : 1.0) * a[i].recipMass) * 
          random->gaussian_vector() );
    }
  } else {
    NAMD_bug("Sequencer::reassignVelocities called improperly!");
  }
}

void Sequencer::reinitVelocities(void)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  BigReal newTemp = simParams->initialTemp;
  BigReal kbT = BOLTZMANN * newTemp;

  int lesReduceTemp = simParams->lesOn && simParams->lesReduceTemp;
  BigReal tempFactor = lesReduceTemp ? 1.0 / simParams->lesFactor : 1.0;

  for ( int i = 0; i < numAtoms; ++i )
  {
    a[i].velocity = ( ( (simParams->fixedAtomsOn && a[i].atomFixed) || 
          a[i].mass <= 0.) ? Vector(0,0,0) :
        sqrt(kbT * (a[i].partition ? tempFactor : 1.0) * a[i].recipMass) * 
        random->gaussian_vector() );
    if ( simParams->drudeOn && i+1 < numAtoms && a[i+1].mass < 1.0 && a[i+1].mass > 0.05 ) {
      a[i+1].velocity = a[i].velocity;  // zero is good enough
      ++i;
    }
  }
}

void Sequencer::rescaleVelocitiesByFactor(BigReal factor)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  for ( int i = 0; i < numAtoms; ++i )
  {
    a[i].velocity *= factor;
  }
}

void Sequencer::reloadCharges()
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Molecule *molecule = Node::Object()->molecule;
  for ( int i = 0; i < numAtoms; ++i )
  {
    a[i].charge = molecule->atomcharge(a[i].id);
  }
}

// REST2 solute charge scaling
void Sequencer::rescaleSoluteCharges(BigReal factor)
{
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  Molecule *molecule = Node::Object()->molecule;
  BigReal sqrt_factor = sqrt(factor);
  // apply scaling to the original charge (stored in molecule)
  // of just the marked solute atoms
  for ( int i = 0; i < numAtoms; ++i ) {
    if (molecule->get_ss_type(a[i].id)) {
      a[i].charge = sqrt_factor * molecule->atomcharge(a[i].id);
    }
  }
}

void Sequencer::tcoupleVelocities(BigReal dt_fs, int step)
{
  if ( simParams->tCoupleOn )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    // Blocking receive for the temperature coupling coefficient.
    BigReal coefficient = broadcast->tcoupleCoefficient.get(step);
    Molecule *molecule = Node::Object()->molecule;
    BigReal dt = dt_fs * 0.001;  // convert to ps
    coefficient *= dt;
    for ( int i = 0; i < numAtoms; ++i )
    {
      BigReal f1 = exp( coefficient * a[i].langevinParam );
      a[i].velocity *= f1;
    }
  }
}

/** Rescale velocities with the scale factor sent from the Controller.
 *
 *  \param dt_fs The integration increment, in fs (not currently used)
 *  \param step The current timestep
 */
void Sequencer::stochRescaleVelocities(BigReal dt_fs, int step)
{
  if ( simParams->stochRescaleOn )
  {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    ++stochRescale_count;
    if ( stochRescale_count == simParams->stochRescaleFreq )
    {
      // Blocking receive for the temperature coupling coefficient.
      BigReal coefficient = broadcast->stochRescaleCoefficient.get(step);
      DebugM(4, "stochastically rescaling velocities at step " << step << " by " << coefficient << "\n");
      for ( int i = 0; i < numAtoms; ++i )
      {
        a[i].velocity *= coefficient;
      }
      stochRescale_count = 0;
    }
  }
}

void Sequencer::saveForce(const int ftag)
{
  patch->saveForce(ftag);
}

//
// DJH: Need to change division by TIMEFACTOR into multiplication by
// reciprocal of TIMEFACTOR. Done several times for each iteration of
// the integrate() loop.
//

void Sequencer::addForceToMomentum(
    BigReal timestep, const int ftag, const int useSaved
    ) {
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::ADD_FORCE_TO_MOMENTUM);
#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif
  const BigReal dt = timestep / TIMEFACTOR;
  FullAtom *atom_arr  = patch->atom.begin();
  ForceList *f_use = (useSaved ? patch->f_saved : patch->f);
  const Force *force_arr = f_use[ftag].const_begin();
  patch->addForceToMomentum(atom_arr, force_arr, dt, patch->numAtoms);
}

void Sequencer::addForceToMomentum3(
    const BigReal timestep1, const int ftag1, const int useSaved1,
    const BigReal timestep2, const int ftag2, const int useSaved2,
    const BigReal timestep3, const int ftag3, const int useSaved3
    ) {
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::ADD_FORCE_TO_MOMENTUM);
#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif
  const BigReal dt1 = timestep1 / TIMEFACTOR;
  const BigReal dt2 = timestep2 / TIMEFACTOR;
  const BigReal dt3 = timestep3 / TIMEFACTOR;
  ForceList *f_use1 = (useSaved1 ? patch->f_saved : patch->f);
  ForceList *f_use2 = (useSaved2 ? patch->f_saved : patch->f);
  ForceList *f_use3 = (useSaved3 ? patch->f_saved : patch->f);
  FullAtom *atom_arr  = patch->atom.begin();
  const Force *force_arr1 = f_use1[ftag1].const_begin();
  const Force *force_arr2 = f_use2[ftag2].const_begin();
  const Force *force_arr3 = f_use3[ftag3].const_begin();
  patch->addForceToMomentum3 (atom_arr, force_arr1, force_arr2, force_arr3,
      dt1, dt2, dt3, patch->numAtoms);
}

void Sequencer::addVelocityToPosition(BigReal timestep)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::ADD_VELOCITY_TO_POSITION);
#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif
  const BigReal dt = timestep / TIMEFACTOR;
  FullAtom *atom_arr  = patch->atom.begin();
  patch->addVelocityToPosition(atom_arr, dt, patch->numAtoms);
}

void Sequencer::hardWallDrude(BigReal dt, int pressure)
{
  if ( simParams->drudeHardWallOn ) {
    Tensor virial;
    Tensor *vp = ( pressure ? &virial : 0 );
    if ( patch->hardWallDrude(dt, vp, pressureProfileReduction) ) {
      iout << iERROR << "Constraint failure in HardWallDrude(); "
        << "simulation may become unstable.\n" << endi;
      Node::Object()->enableEarlyExit();
      terminate();
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
  }
}

void Sequencer::rattle1(BigReal dt, int pressure)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on, NamdProfileEvent::RATTLE1);
  if ( simParams->rigidBonds != RIGID_NONE ) {
    Tensor virial;
    Tensor *vp = ( pressure ? &virial : 0 );
    if ( patch->rattle1(dt, vp, pressureProfileReduction) ) {
      iout << iERROR << 
        "Constraint failure; simulation has become unstable.\n" << endi;
      Node::Object()->enableEarlyExit();
      terminate();
    }
#if 0
    printf("virial = %g %g %g  %g %g %g  %g %g %g\n",
        virial.xx, virial.xy, virial.xz,
        virial.yx, virial.yy, virial.yz,
        virial.zx, virial.zy, virial.zz);
#endif
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
#if 0
    {
      const FullAtom *a = patch->atom.const_begin();
      for (int n=0;  n < patch->numAtoms;  n++) {
        printf("pos[%d] =  %g %g %g\n", n,
            a[n].position.x, a[n].position.y, a[n].position.z);
      }
      for (int n=0;  n < patch->numAtoms;  n++) {
        printf("vel[%d] =  %g %g %g\n", n,
            a[n].velocity.x, a[n].velocity.y, a[n].velocity.z);
      }
      if (pressure) {
        for (int n=0;  n < patch->numAtoms;  n++) {
          printf("force[%d] =  %g %g %g\n", n,
              patch->f[Results::normal][n].x,
              patch->f[Results::normal][n].y,
              patch->f[Results::normal][n].z);
        }
      }
    }
#endif
  }
}

// void Sequencer::rattle2(BigReal dt, int step)
// {
//   if ( simParams->rigidBonds != RIGID_NONE ) {
//     Tensor virial;
//     patch->rattle2(dt, &virial);
//     ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
//     // we need to add to alt and int virial because not included in forces
// #ifdef ALTVIRIAL
//     ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
// #endif
//     ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,virial);
//   }
// }

void Sequencer::maximumMove(BigReal timestep)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on, NamdProfileEvent::MAXIMUM_MOVE);

  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;
  if ( simParams->maximumMove ) {
    const BigReal dt = timestep / TIMEFACTOR;
    const BigReal maxvel = simParams->maximumMove / dt;
    const BigReal maxvel2 = maxvel * maxvel;
    for ( int i=0; i<numAtoms; ++i ) {
      if ( a[i].velocity.length2() > maxvel2 ) {
	a[i].velocity *= ( maxvel / a[i].velocity.length() );
      }
    }
  } else {
    const BigReal dt = timestep / TIMEFACTOR;
    const BigReal maxvel = simParams->cutoff / dt;
    const BigReal maxvel2 = maxvel * maxvel;
    int killme = 0;
    for ( int i=0; i<numAtoms; ++i ) {
      killme = killme || ( a[i].velocity.length2() > maxvel2 );
    }
    if ( killme ) {
      killme = 0;
      for ( int i=0; i<numAtoms; ++i ) {
        if ( a[i].velocity.length2() > maxvel2 ) {
          ++killme;
          iout << iERROR << "Atom " << (a[i].id + 1) << " velocity is "
            << ( PDBVELFACTOR * a[i].velocity ) << " (limit is "
            << ( PDBVELFACTOR * maxvel ) << ", atom "
            << i << " of " << numAtoms << " on patch "
            << patch->patchID << " pe " << CkMyPe() << ")\n" << endi;
        }
      }
      iout << iERROR << 
        "Atoms moving too fast; simulation has become unstable ("
        << killme << " atoms on patch " << patch->patchID
        << " pe " << CkMyPe() << ").\n" << endi;
      Node::Object()->enableEarlyExit();
      terminate();
    }
  }
}

void Sequencer::minimizationQuenchVelocity(void)
{
  if ( simParams->minimizeOn ) {
    FullAtom *a = patch->atom.begin();
    int numAtoms = patch->numAtoms;
    for ( int i=0; i<numAtoms; ++i ) {
      a[i].velocity = 0.;
    }
  }
}

void Sequencer::submitHalfstep(int step)
{
  NAMD_EVENT_RANGE_2(patch->flags.event_on, NamdProfileEvent::SUBMIT_HALFSTEP);

  // velocity-dependent quantities *** ONLY ***
  // positions are not at half-step when called
  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

#if CMK_BLUEGENEL
  CmiNetworkProgressAfter (0);
#endif

  // For non-Multigrator doKineticEnergy = 1 always
  Tensor momentumSqrSum;
  if (doKineticEnergy || patch->flags.doVirial)
  {
    BigReal kineticEnergy = 0;
    Tensor virial;
    if ( simParams->pairInteractionOn ) {
      if ( simParams->pairInteractionSelf ) {
        for ( int i = 0; i < numAtoms; ++i ) {
          if ( a[i].partition != 1 ) continue;
          kineticEnergy += a[i].mass * a[i].velocity.length2();
          virial.outerAdd(a[i].mass, a[i].velocity, a[i].velocity);
        }
      }
    } else {
      for ( int i = 0; i < numAtoms; ++i ) {
        if (a[i].mass < 0.01) continue;
        kineticEnergy += a[i].mass * a[i].velocity.length2();
        virial.outerAdd(a[i].mass, a[i].velocity, a[i].velocity);
      }
    }
    
    if (simParams->multigratorOn && !simParams->useGroupPressure) {
      momentumSqrSum = virial;
    }
    kineticEnergy *= 0.5 * 0.5;
    reduction->item(REDUCTION_HALFSTEP_KINETIC_ENERGY) += kineticEnergy;
    virial *= 0.5;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,virial);
#ifdef ALTVIRIAL
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,virial);
#endif
  }
 
  if (pressureProfileReduction) {
    int nslabs = simParams->pressureProfileSlabs;
    const Lattice &lattice = patch->lattice;
    BigReal idz = nslabs/lattice.c().z;
    BigReal zmin = lattice.origin().z - 0.5*lattice.c().z;
    int useGroupPressure = simParams->useGroupPressure;

    // Compute kinetic energy partition, possibly subtracting off
    // internal kinetic energy if group pressure is enabled.
    // Since the regular pressure is 1/2 mvv and the internal kinetic
    // term that is subtracted off for the group pressure is
    // 1/2 mv (v-v_cm), the group pressure kinetic contribution is
    // 1/2 m * v * v_cm.  The factor of 1/2 is because submitHalfstep
    // gets called twice per timestep.
    int hgs;
    for (int i=0; i<numAtoms; i += hgs) {
      int j, ppoffset;
      hgs = a[i].hydrogenGroupSize;
      int partition = a[i].partition;

      BigReal m_cm = 0;
      Velocity v_cm(0,0,0);
      for (j=i; j< i+hgs; ++j) {
        m_cm += a[j].mass;
        v_cm += a[j].mass * a[j].velocity;
      }
      v_cm /= m_cm;
      for (j=i; j < i+hgs; ++j) {
        BigReal mass = a[j].mass;
        if (! (useGroupPressure && j != i)) {
          BigReal z = a[j].position.z;
          int slab = (int)floor((z-zmin)*idz);
          if (slab < 0) slab += nslabs;
          else if (slab >= nslabs) slab -= nslabs;
          ppoffset = 3*(slab + partition*nslabs);
        }
        BigReal wxx, wyy, wzz;
        if (useGroupPressure) {
          wxx = 0.5*mass * a[j].velocity.x * v_cm.x;
          wyy = 0.5*mass * a[j].velocity.y * v_cm.y;
          wzz = 0.5*mass * a[j].velocity.z * v_cm.z;
        } else {
          wxx = 0.5*mass * a[j].velocity.x * a[j].velocity.x;
          wyy = 0.5*mass * a[j].velocity.y * a[j].velocity.y;
          wzz = 0.5*mass * a[j].velocity.z * a[j].velocity.z;
        }
        pressureProfileReduction->item(ppoffset  ) += wxx;
        pressureProfileReduction->item(ppoffset+1) += wyy;
        pressureProfileReduction->item(ppoffset+2) += wzz;
      }
    }
  } 

  // For non-Multigrator doKineticEnergy = 1 always
  if (doKineticEnergy || patch->flags.doVirial)
  {
    BigReal intKineticEnergy = 0;
    Tensor intVirialNormal;

    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {

#if CMK_BLUEGENEL
      CmiNetworkProgress ();
#endif

      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        v_cm += a[j].mass * a[j].velocity;
      }
      if (simParams->multigratorOn && simParams->useGroupPressure) {
        momentumSqrSum.outerAdd(1.0/m_cm, v_cm, v_cm);
      }
      v_cm /= m_cm;
      if ( simParams->pairInteractionOn ) {
        if ( simParams->pairInteractionSelf ) {
          for ( j = i; j < (i+hgs); ++j ) {
            if ( a[j].partition != 1 ) continue;
            BigReal mass = a[j].mass;
            Vector v = a[j].velocity;
            Vector dv = v - v_cm;
            intKineticEnergy += mass * (v * dv);
            intVirialNormal.outerAdd (mass, v, dv);
          }
        }
      } else {
        for ( j = i; j < (i+hgs); ++j ) {
          BigReal mass = a[j].mass;
          Vector v = a[j].velocity;
          Vector dv = v - v_cm;
          intKineticEnergy += mass * (v * dv);
          intVirialNormal.outerAdd(mass, v, dv);
        }
      }
    }

    intKineticEnergy *= 0.5 * 0.5;
    reduction->item(REDUCTION_INT_HALFSTEP_KINETIC_ENERGY) += intKineticEnergy;
    intVirialNormal *= 0.5;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
    if ( simParams->multigratorOn) {
      momentumSqrSum *= 0.5;
      ADD_TENSOR_OBJECT(reduction,REDUCTION_MOMENTUM_SQUARED,momentumSqrSum);
    }
  }

}

void Sequencer::calcFixVirial(Tensor& fixVirialNormal, Tensor& fixVirialNbond, Tensor& fixVirialSlow,
  Vector& fixForceNormal, Vector& fixForceNbond, Vector& fixForceSlow) {

  FullAtom *a = patch->atom.begin();
  int numAtoms = patch->numAtoms;

  for ( int j = 0; j < numAtoms; j++ ) {
    if ( simParams->fixedAtomsOn && a[j].atomFixed ) {
      Vector dx = a[j].fixedPosition;
      // all negative because fixed atoms cancels these forces
      fixVirialNormal.outerAdd(-1.0, patch->f[Results::normal][j], dx);
      fixVirialNbond.outerAdd(-1.0, patch->f[Results::nbond][j], dx);
      fixVirialSlow.outerAdd(-1.0, patch->f[Results::slow][j], dx);
      fixForceNormal -= patch->f[Results::normal][j];
      fixForceNbond -= patch->f[Results::nbond][j];
      fixForceSlow -= patch->f[Results::slow][j];
    }
  }
}

void Sequencer::submitReductions(int step)
{
#ifndef UPPER_BOUND
  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::SUBMIT_REDUCTIONS);
  FullAtom *a = patch->atom.begin();
#endif
  int numAtoms = patch->numAtoms;

#if CMK_BLUEGENEL
  CmiNetworkProgressAfter(0);
#endif

  reduction->item(REDUCTION_ATOM_CHECKSUM) += numAtoms;
  reduction->item(REDUCTION_MARGIN_VIOLATIONS) += patch->marginViolations;

#ifndef UPPER_BOUND
  // For non-Multigrator doKineticEnergy = 1 always
  if (doKineticEnergy || doMomenta || patch->flags.doVirial)
  {
    BigReal kineticEnergy = 0;
    Vector momentum = 0;
    Vector angularMomentum = 0;
    Vector o = patch->lattice.origin();
    int i;
    if ( simParams->pairInteractionOn ) {
      if ( simParams->pairInteractionSelf ) {
        for (i = 0; i < numAtoms; ++i ) {
          if ( a[i].partition != 1 ) continue;
          kineticEnergy += a[i].mass * a[i].velocity.length2();
          momentum += a[i].mass * a[i].velocity;
          angularMomentum += cross(a[i].mass,a[i].position-o,a[i].velocity);
        }
      }
    } else {
      for (i = 0; i < numAtoms; ++i ) {
        kineticEnergy += a[i].mass * a[i].velocity.length2();
        momentum += a[i].mass * a[i].velocity;
        angularMomentum += cross(a[i].mass,a[i].position-o,a[i].velocity);
      }
      if (simParams->drudeOn) {
        BigReal drudeComKE = 0.;
        BigReal drudeBondKE = 0.;

        for (i = 0;  i < numAtoms;  i++) {
          if (i < numAtoms-1 &&
              a[i+1].mass < 1.0 && a[i+1].mass > 0.05) {
            // i+1 is a Drude particle with parent i

            // convert from Cartesian coordinates to (COM,bond) coordinates
            BigReal m_com = (a[i].mass + a[i+1].mass);  // mass of COM
            BigReal m = a[i+1].mass / m_com;  // mass ratio
            BigReal m_bond = a[i+1].mass * (1. - m);  // mass of bond
            Vector v_bond = a[i+1].velocity - a[i].velocity;  // vel of bond
            Vector v_com = a[i].velocity + m * v_bond;  // vel of COM

            drudeComKE += m_com * v_com.length2();
            drudeBondKE += m_bond * v_bond.length2();

            i++;  // +1 from loop, we've updated both particles
          }
          else {
            drudeComKE += a[i].mass * a[i].velocity.length2();
          }
        } // end for

        drudeComKE *= 0.5;
        drudeBondKE *= 0.5;
        reduction->item(REDUCTION_DRUDECOM_CENTERED_KINETIC_ENERGY)
          += drudeComKE;
        reduction->item(REDUCTION_DRUDEBOND_CENTERED_KINETIC_ENERGY)
          += drudeBondKE;
      } // end drudeOn

    } // end else

    kineticEnergy *= 0.5;
    reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY) += kineticEnergy;
    ADD_VECTOR_OBJECT(reduction,REDUCTION_MOMENTUM,momentum);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_ANGULAR_MOMENTUM,angularMomentum);  
  }

#ifdef ALTVIRIAL
  // THIS IS NOT CORRECTED FOR PAIR INTERACTIONS
  {
    Tensor altVirial;
    for ( int i = 0; i < numAtoms; ++i ) {
      altVirial.outerAdd(1.0, patch->f[Results::normal][i], a[i].position);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NORMAL,altVirial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < numAtoms; ++i ) {
      altVirial.outerAdd(1.0, patch->f[Results::nbond][i], a[i].position);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_NBOND,altVirial);
  }
  {
    Tensor altVirial;
    for ( int i = 0; i < numAtoms; ++i ) {
      altVirial.outerAdd(1.0, patch->f[Results::slow][i], a[i].position);
    }
    ADD_TENSOR_OBJECT(reduction,REDUCTION_ALT_VIRIAL_SLOW,altVirial);
  }
#endif

  // For non-Multigrator doKineticEnergy = 1 always
  if (doKineticEnergy || patch->flags.doVirial)
  {
    BigReal intKineticEnergy = 0;
    Tensor intVirialNormal;
    Tensor intVirialNbond;
    Tensor intVirialSlow;

    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
#if CMK_BLUEGENEL
      CmiNetworkProgress();
#endif
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      Velocity v_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
        v_cm += a[j].mass * a[j].velocity;
      }
      x_cm /= m_cm;
      v_cm /= m_cm;
      int fixedAtomsOn = simParams->fixedAtomsOn;
      if ( simParams->pairInteractionOn ) {
        int pairInteractionSelf = simParams->pairInteractionSelf;
        for ( j = i; j < (i+hgs); ++j ) {
          if ( a[j].partition != 1 &&
               ( pairInteractionSelf || a[j].partition != 2 ) ) continue;
          // net force treated as zero for fixed atoms
          if ( fixedAtomsOn && a[j].atomFixed ) continue;
          BigReal mass = a[j].mass;
          Vector v = a[j].velocity;
          Vector dv = v - v_cm;
          intKineticEnergy += mass * (v * dv);
          Vector dx = a[j].position - x_cm;
          intVirialNormal.outerAdd(1.0, patch->f[Results::normal][j], dx);
          intVirialNbond.outerAdd(1.0, patch->f[Results::nbond][j], dx);
          intVirialSlow.outerAdd(1.0, patch->f[Results::slow][j], dx);
        }
      } else {
        for ( j = i; j < (i+hgs); ++j ) {
          // net force treated as zero for fixed atoms
          if ( fixedAtomsOn && a[j].atomFixed ) continue;
          BigReal mass = a[j].mass;
          Vector v = a[j].velocity;
          Vector dv = v - v_cm;
          intKineticEnergy += mass * (v * dv);
          Vector dx = a[j].position - x_cm;
          intVirialNormal.outerAdd(1.0, patch->f[Results::normal][j], dx);
          intVirialNbond.outerAdd(1.0, patch->f[Results::nbond][j], dx);
          intVirialSlow.outerAdd(1.0, patch->f[Results::slow][j], dx);
        }
      }
    }

    intKineticEnergy *= 0.5;
    reduction->item(REDUCTION_INT_CENTERED_KINETIC_ENERGY) += intKineticEnergy;
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NBOND,intVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_SLOW,intVirialSlow);
  }

  if (pressureProfileReduction && simParams->useGroupPressure) {
    // subtract off internal virial term, calculated as for intVirial.
    int nslabs = simParams->pressureProfileSlabs;
    const Lattice &lattice = patch->lattice;
    BigReal idz = nslabs/lattice.c().z;
    BigReal zmin = lattice.origin().z - 0.5*lattice.c().z;
    int useGroupPressure = simParams->useGroupPressure;

    int hgs;
    for (int i=0; i<numAtoms; i += hgs) {
      int j;
      hgs = a[i].hydrogenGroupSize;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      for (j=i; j< i+hgs; ++j) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
      }
      x_cm /= m_cm;
      
      BigReal z = a[i].position.z;
      int slab = (int)floor((z-zmin)*idz);
      if (slab < 0) slab += nslabs;
      else if (slab >= nslabs) slab -= nslabs;
      int partition = a[i].partition;
      int ppoffset = 3*(slab + nslabs*partition);
      for (j=i; j < i+hgs; ++j) {
        BigReal mass = a[j].mass;
        Vector dx = a[j].position - x_cm;
        const Vector &fnormal = patch->f[Results::normal][j];
        const Vector &fnbond  = patch->f[Results::nbond][j];
        const Vector &fslow   = patch->f[Results::slow][j];
        BigReal wxx = (fnormal.x + fnbond.x + fslow.x) * dx.x;
        BigReal wyy = (fnormal.y + fnbond.y + fslow.y) * dx.y;
        BigReal wzz = (fnormal.z + fnbond.z + fslow.z) * dx.z;
        pressureProfileReduction->item(ppoffset  ) -= wxx;
        pressureProfileReduction->item(ppoffset+1) -= wyy;
        pressureProfileReduction->item(ppoffset+2) -= wzz;
      }
    }
  }

  // For non-Multigrator doVirial = 1 always
  if (patch->flags.doVirial)
  {
    if ( simParams->fixedAtomsOn ) {
      Tensor fixVirialNormal;
      Tensor fixVirialNbond;
      Tensor fixVirialSlow;
      Vector fixForceNormal = 0;
      Vector fixForceNbond = 0;
      Vector fixForceSlow = 0;

      calcFixVirial(fixVirialNormal, fixVirialNbond, fixVirialSlow, fixForceNormal, fixForceNbond, fixForceSlow);

      ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,fixVirialNormal);
      ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,fixVirialNbond);
      ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,fixVirialSlow);
      ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,fixForceNormal);
      ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NBOND,fixForceNbond);
      ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_SLOW,fixForceSlow);
    }
  }
#endif // UPPER_BOUND

  reduction->submit();
#ifndef UPPER_BOUND
  if (pressureProfileReduction) pressureProfileReduction->submit();
#endif
}

void Sequencer::submitMinimizeReductions(int step, BigReal fmax2)
{
  FullAtom *a = patch->atom.begin();
  Force *f1 = patch->f[Results::normal].begin();
  Force *f2 = patch->f[Results::nbond].begin();
  Force *f3 = patch->f[Results::slow].begin();
  const bool fixedAtomsOn = simParams->fixedAtomsOn;
  const bool drudeHardWallOn = simParams->drudeHardWallOn;
  const double drudeBondLen = simParams->drudeBondLen;
  const double drudeBondLen2 = drudeBondLen * drudeBondLen;
  const double drudeStep = 0.1/(TIMEFACTOR*TIMEFACTOR);
  const double drudeMove = 0.01;
  const double drudeStep2 = drudeStep * drudeStep;
  const double drudeMove2 = drudeMove * drudeMove;
  int numAtoms = patch->numAtoms;

  reduction->item(REDUCTION_ATOM_CHECKSUM) += numAtoms;

  for ( int i = 0; i < numAtoms; ++i ) {
    f1[i] += f2[i] + f3[i];  // add all forces
    if ( drudeHardWallOn && i && (a[i].mass > 0.05) && ((a[i].mass < 1.0)) ) { // drude particle
      if ( ! fixedAtomsOn || ! a[i].atomFixed ) {
        if ( drudeStep2 * f1[i].length2() > drudeMove2 ) {
          a[i].position += drudeMove * f1[i].unit();
        } else {
          a[i].position += drudeStep * f1[i];
        }
        if ( (a[i].position - a[i-1].position).length2() > drudeBondLen2 ) {
          a[i].position = a[i-1].position + drudeBondLen * (a[i].position - a[i-1].position).unit();
        }
      }
      Vector netf = f1[i-1] + f1[i];
      if ( fixedAtomsOn && a[i-1].atomFixed ) netf = 0;
      f1[i-1] = netf;
      f1[i] = 0.;
    }
    if ( fixedAtomsOn && a[i].atomFixed ) f1[i] = 0;
  }

  f2 = f3 = 0;  // included in f1

  BigReal maxv2 = 0.;

  for ( int i = 0; i < numAtoms; ++i ) {
    BigReal v2 = a[i].velocity.length2();
    if ( v2 > 0. ) {
      if ( v2 > maxv2 ) maxv2 = v2;
    } else {
      v2 = f1[i].length2();
      if ( v2 > maxv2 ) maxv2 = v2;
    }
  }

  if ( fmax2 > 10. * TIMEFACTOR * TIMEFACTOR * TIMEFACTOR * TIMEFACTOR )
  { Tensor virial; patch->minimize_rattle2( 0.1 * TIMEFACTOR / sqrt(maxv2), &virial, true /* forces */); }

  BigReal fdotf = 0;
  BigReal fdotv = 0;
  BigReal vdotv = 0;
  int numHuge = 0;
  for ( int i = 0; i < numAtoms; ++i ) {
    if ( simParams->fixedAtomsOn && a[i].atomFixed ) continue;
    if ( drudeHardWallOn && (a[i].mass > 0.05) && ((a[i].mass < 1.0)) ) continue; // drude particle
    Force f = f1[i];
    BigReal ff = f * f;
    if ( ff > fmax2 ) {
      if (simParams->printBadContacts) {
        CkPrintf("STEP(%i) MIN_HUGE[%i] f=%e kcal/mol/A\n",patch->flags.sequence,patch->pExt[i].id,ff);
      }
      ++numHuge;
      // pad scaling so minimizeMoveDownhill() doesn't miss them
      BigReal fmult = 1.01 * sqrt(fmax2/ff);
      f *= fmult;  ff = f * f;
      f1[i] *= fmult;
    }
    fdotf += ff;
    fdotv += f * a[i].velocity;
    vdotv += a[i].velocity * a[i].velocity;
  }

  reduction->item(REDUCTION_MIN_F_DOT_F) += fdotf;
  reduction->item(REDUCTION_MIN_F_DOT_V) += fdotv;
  reduction->item(REDUCTION_MIN_V_DOT_V) += vdotv;
  reduction->item(REDUCTION_MIN_HUGE_COUNT) += numHuge;

  {
    Tensor intVirialNormal;
    Tensor intVirialNbond;
    Tensor intVirialSlow;

    int hgs;
    for ( int i = 0; i < numAtoms; i += hgs ) {
      hgs = a[i].hydrogenGroupSize;
      int j;
      BigReal m_cm = 0;
      Position x_cm(0,0,0);
      for ( j = i; j < (i+hgs); ++j ) {
        m_cm += a[j].mass;
        x_cm += a[j].mass * a[j].position;
      }
      x_cm /= m_cm;
      for ( j = i; j < (i+hgs); ++j ) {
        BigReal mass = a[j].mass;
	// net force treated as zero for fixed atoms
        if ( simParams->fixedAtomsOn && a[j].atomFixed ) continue;
        Vector dx = a[j].position - x_cm;
        intVirialNormal.outerAdd(1.0, patch->f[Results::normal][j], dx);
        intVirialNbond.outerAdd(1.0, patch->f[Results::nbond][j], dx);
        intVirialSlow.outerAdd(1.0, patch->f[Results::slow][j], dx);
      }
    }

    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NORMAL,intVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_NBOND,intVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_INT_VIRIAL_SLOW,intVirialSlow);
  }

  if ( simParams->fixedAtomsOn ) {
    Tensor fixVirialNormal;
    Tensor fixVirialNbond;
    Tensor fixVirialSlow;
    Vector fixForceNormal = 0;
    Vector fixForceNbond = 0;
    Vector fixForceSlow = 0;

    calcFixVirial(fixVirialNormal, fixVirialNbond, fixVirialSlow, fixForceNormal, fixForceNbond, fixForceSlow);

    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NORMAL,fixVirialNormal);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_NBOND,fixVirialNbond);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,fixVirialSlow);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NORMAL,fixForceNormal);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_NBOND,fixForceNbond);
    ADD_VECTOR_OBJECT(reduction,REDUCTION_EXT_FORCE_SLOW,fixForceSlow);
  }

  reduction->submit();
}

void Sequencer::submitCollections(int step, int zeroVel)
{
  //
  // DJH: Copy updates of SOA back into AOS.
  // Do we need to update everything or is it safe to just update
  // positions and velocities separately, as needed?
  //
  //patch->copy_updates_to_AOS();

  NAMD_EVENT_RANGE_2(patch->flags.event_on,
      NamdProfileEvent::SUBMIT_COLLECTIONS);
  int prec = Output::coordinateNeeded(step);
  if ( prec ) {
    collection->submitPositions(step,patch->atom,patch->lattice,prec);
  }
  if ( Output::velocityNeeded(step) ) {
    collection->submitVelocities(step,zeroVel,patch->atom);
  }
  if ( Output::forceNeeded(step) ) {
    int maxForceUsed = patch->flags.maxForceUsed;
    if ( maxForceUsed > Results::slow ) maxForceUsed = Results::slow;
    collection->submitForces(step,patch->atom,maxForceUsed,patch->f);
  }
}

void Sequencer::runComputeObjects(int migration, int pairlists, int pressureStep)
{
  if ( migration ) pairlistsAreValid = 0;
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
  if ( pairlistsAreValid &&
       ( patch->flags.doFullElectrostatics || ! simParams->fullElectFrequency )
                         && ( pairlistsAge > (
#else
  if ( pairlistsAreValid && ( pairlistsAge > (
#endif
         (simParams->stepsPerCycle - 1) / simParams->pairlistsPerCycle ) ) ) {
    pairlistsAreValid = 0;
  }
  if ( ! simParams->usePairlists ) pairlists = 0;
  patch->flags.usePairlists = pairlists || pairlistsAreValid;
  patch->flags.savePairlists =
	pairlists && ! pairlistsAreValid;

  if ( simParams->singleTopology ) patch->reposition_all_alchpairs();
  if ( simParams->lonepairs ) patch->reposition_all_lonepairs();

  //
  // DJH: Copy updates of SOA back into AOS.
  // The positionsReady() routine starts force computation and atom migration.
  //
  // We could reduce amount of copying here by checking migration status
  // and copying velocities only when migrating. Some types of simulation
  // always require velocities, such as Lowe-Anderson.
  //
  //patch->copy_updates_to_AOS();

  patch->positionsReady(migration);  // updates flags.sequence

  int seq = patch->flags.sequence;
  int basePriority = ( (seq & 0xffff) << 15 )
                     + PATCH_PRIORITY(patch->getPatchID());
  if ( patch->flags.doGBIS && patch->flags.doNonbonded) {
    priority = basePriority + GB1_COMPUTE_HOME_PRIORITY;
    suspend(); // until all deposit boxes close
    patch->gbisComputeAfterP1();
    priority = basePriority + GB2_COMPUTE_HOME_PRIORITY;
    suspend();
    patch->gbisComputeAfterP2();
    priority = basePriority + COMPUTE_HOME_PRIORITY;
    suspend();
  } else {
    priority = basePriority + COMPUTE_HOME_PRIORITY;
    suspend(); // until all deposit boxes close
  }

  //
  // DJH: Copy all data into SOA from AOS.
  //
  // We need everything copied after atom migration.
  // When doing force computation without atom migration,
  // all data except forces will already be up-to-date in SOA
  // (except maybe for some special types of simulation).
  //
  //patch->copy_all_to_SOA();

  //
  // DJH: Copy forces to SOA.
  // Force available after suspend() has returned.
  //
  //patch->copy_forces_to_SOA();

  if ( patch->flags.savePairlists && patch->flags.doNonbonded ) {
    pairlistsAreValid = 1;
    pairlistsAge = 0;
  }
  // For multigrator, do not age pairlist during pressure step
  // NOTE: for non-multigrator pressureStep = 0 always
  if ( pairlistsAreValid && !pressureStep ) ++pairlistsAge;

  if (simParams->lonepairs) {
    {
      Tensor virial;
      patch->redistrib_lonepair_forces(Results::normal, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NORMAL, virial);
    }
    if (patch->flags.doNonbonded) {
      Tensor virial;
      patch->redistrib_lonepair_forces(Results::nbond, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, virial);
    }
    if (patch->flags.doFullElectrostatics) {
      Tensor virial;
      patch->redistrib_lonepair_forces(Results::slow, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, virial);
    }
  } else if (simParams->watmodel == WAT_TIP4) {
    {
      Tensor virial;
      patch->redistrib_tip4p_forces(Results::normal, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NORMAL, virial);
    }
    if (patch->flags.doNonbonded) {
      Tensor virial;
      patch->redistrib_tip4p_forces(Results::nbond, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, virial);
    }
    if (patch->flags.doFullElectrostatics) {
      Tensor virial;
      patch->redistrib_tip4p_forces(Results::slow, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, virial);
    }
  } else if (simParams->watmodel == WAT_SWM4) {
    {
      Tensor virial;
      patch->redistrib_swm4_forces(Results::normal, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NORMAL, virial);
    }
    if (patch->flags.doNonbonded) {
      Tensor virial;
      patch->redistrib_swm4_forces(Results::nbond, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_NBOND, virial);
    }
    if (patch->flags.doFullElectrostatics) {
      Tensor virial;
      patch->redistrib_swm4_forces(Results::slow, &virial);
      ADD_TENSOR_OBJECT(reduction, REDUCTION_VIRIAL_SLOW, virial);
    }
  }

  if (simParams->singleTopology) {
    patch->redistrib_alchpair_forces(Results::normal);
    if (patch->flags.doNonbonded) {
      patch->redistrib_alchpair_forces(Results::nbond);
    }
    if (patch->flags.doFullElectrostatics) {
      patch->redistrib_alchpair_forces(Results::slow);
    }
  }

  if ( patch->flags.doMolly ) {
    Tensor virial;
    patch->mollyMollify(&virial);
    ADD_TENSOR_OBJECT(reduction,REDUCTION_VIRIAL_SLOW,virial);
  }


  // BEGIN LA
  if (patch->flags.doLoweAndersen) {
      patch->loweAndersenFinish();
  }
  // END LA
#ifdef NAMD_CUDA_XXX
  int numAtoms = patch->numAtoms;
  FullAtom *a = patch->atom.begin();
  for ( int i=0; i<numAtoms; ++i ) {
    CkPrintf("%d %g %g %g\n", a[i].id,
        patch->f[Results::normal][i].x +
        patch->f[Results::nbond][i].x +
        patch->f[Results::slow][i].x,
        patch->f[Results::normal][i].y + 
        patch->f[Results::nbond][i].y +
        patch->f[Results::slow][i].y,
        patch->f[Results::normal][i].z +
        patch->f[Results::nbond][i].z +
        patch->f[Results::slow][i].z);
    CkPrintf("%d %g %g %g\n", a[i].id,
        patch->f[Results::normal][i].x,
        patch->f[Results::nbond][i].x,
        patch->f[Results::slow][i].x);
    CkPrintf("%d %g %g %g\n", a[i].id,
        patch->f[Results::normal][i].y,
        patch->f[Results::nbond][i].y,
        patch->f[Results::slow][i].y);
    CkPrintf("%d %g %g %g\n", a[i].id,
        patch->f[Results::normal][i].z,
        patch->f[Results::nbond][i].z,
        patch->f[Results::slow][i].z);
  }
#endif

#if PRINT_FORCES
  int numAtoms = patch->numAtoms;
  FullAtom *a = patch->atom.begin();
  for ( int i=0; i<numAtoms; ++i ) {
    float fxNo = patch->f[Results::normal][i].x;
    float fxNb = patch->f[Results::nbond][i].x;
    float fxSl = patch->f[Results::slow][i].x;
    float fyNo = patch->f[Results::normal][i].y;
    float fyNb = patch->f[Results::nbond][i].y;
    float fySl = patch->f[Results::slow][i].y;
    float fzNo = patch->f[Results::normal][i].z;
    float fzNb = patch->f[Results::nbond][i].z;
    float fzSl = patch->f[Results::slow][i].z;
    float fx = fxNo+fxNb+fxSl;
    float fy = fyNo+fyNb+fySl;
    float fz = fzNo+fzNb+fzSl;

		float f = sqrt(fx*fx+fy*fy+fz*fz);
    int id = patch->pExt[i].id;
    int seq = patch->flags.sequence;
    float x = patch->p[i].position.x;
    float y = patch->p[i].position.y;
    float z = patch->p[i].position.z;
    //CkPrintf("FORCE(%04i)[%04i] = <% .4e, % .4e, % .4e> <% .4e, % .4e, % .4e> <% .4e, % .4e, % .4e> <<% .4e, % .4e, % .4e>>\n", seq,id,
    CkPrintf("FORCE(%04i)[%04i] = % .9e % .9e % .9e\n", seq,id,
    //CkPrintf("FORCE(%04i)[%04i] = <% .4e, % .4e, % .4e> <% .4e, % .4e, % .4e> <% .4e, % .4e, % .4e>\n", seq,id,
//fxNo,fyNo,fzNo,
//fxNb,fyNb,fzNb,
//fxSl,fySl,fzSl,
fx,fy,fz
);
	}
#endif
}

void Sequencer::rebalanceLoad(int timestep) {
  if ( ! ldbSteps ) {
    ldbSteps = LdbCoordinator::Object()->getNumStepsToRun();
  }
  if ( ! --ldbSteps ) {
    patch->submitLoadStats(timestep);
    ldbCoordinator->rebalance(this,patch->getPatchID());
    pairlistsAreValid = 0;
  }
}

void Sequencer::cycleBarrier(int doBarrier, int step) {
#if USE_BARRIER
	if (doBarrier)
          // Blocking receive for the cycle barrier.
	  broadcast->cycleBarrier.get(step);
#endif
}

void Sequencer::traceBarrier(int step){
        // Blocking receive for the trace barrier.
	broadcast->traceBarrier.get(step);
}

#ifdef MEASURE_NAMD_WITH_PAPI
void Sequencer::papiMeasureBarrier(int step){
        // Blocking receive for the PAPI measure barrier.
	broadcast->papiMeasureBarrier.get(step);
}
#endif

void Sequencer::terminate() {
  LdbCoordinator::Object()->pauseWork(patch->ldObjHandle);
  CthFree(thread);
  CthSuspend();
}

