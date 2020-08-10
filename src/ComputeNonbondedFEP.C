/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Common operations for ComputeNonbonded classes
*/

#include "ComputeNonbondedInl.h"
#include "ComputeNonbondedAlch.h"

/* ******************************************** */
/* vdW energy, force and lambda2 energy for FEP */
/* ******************************************** */
inline void fep_vdw_forceandenergies (BigReal A, BigReal B, BigReal r2, 
    BigReal myVdwShift, BigReal myVdwShift2, BigReal switchdist2, 
    BigReal cutoff2, BigReal switchfactor, Bool vdwForceSwitching,
    BigReal myVdwLambda, BigReal myVdwLambda2, Bool alchWCAOn,
    BigReal myRepLambda, BigReal myRepLambda2, BigReal* alch_vdw_energy,
    BigReal* alch_vdw_force, BigReal* alch_vdw_energy_2) {
  BigReal U, F, U_2, dU, dU_2, switchmul, switchmul2;
  /*
   * The variable naming here is unavoidably awful. A number immediately after
   * a variable implies a distance to that power. The suffix "_2" indicates an
   * evaluation at alchLambda2. The prefix "inv" means the inverse (previously
   * this also used underscores - it was awful).
   */
  if (alchWCAOn) {
    // WCA-on, auxilliary, lambda-dependent cutoff based on Rmin
    //
    // Avoid divide by zero - corectly zeroes interaction below.
    const BigReal Rmin2 = (B <= 0.0 ? 0.0 : powf(2.0*A/B, 1.f/3));
    if (myRepLambda < 1.0) {
      // modified repulsive soft-core
      // All of the scaling is baked into the shift and cutoff values. There
      // are no linear coupling terms out front.
      const BigReal WCAshift = Rmin2*(1 - myRepLambda)*(1 - myRepLambda);
      if (r2 <= Rmin2 - WCAshift) {
        const BigReal epsilon = B*B/(4.0*A);
        vdw_forceandenergy(A, B, r2 + WCAshift, &U, &F);
        *alch_vdw_energy = U + epsilon;
        *alch_vdw_force = F;
      } else {
        *alch_vdw_energy = 0.0;
        *alch_vdw_force = 0.0;
      }
    } else { // myRepLambda == 1.0
      if (vdwForceSwitching) {
        // force and potential switching
        if (r2 <= Rmin2) {
          // normal LJ, potential w/two shifts
          const BigReal epsilon = B*B/(4.0*A);
          vdw_fswitch_shift(A, B, switchdist2, cutoff2, &dU);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = U + (1 - myVdwLambda)*epsilon + myVdwLambda*dU;
          *alch_vdw_force = F;
        } else if (r2 <= switchdist2) {
          // normal LJ potential w/shift
          vdw_fswitch_shift(A, B, switchdist2, cutoff2, &dU);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = myVdwLambda*(U + dU);
          *alch_vdw_force = myVdwLambda*F;
        } else { // r2 > switchdist
          // normal LJ potential with linear coupling   
          vdw_fswitch_forceandenergy(A, B, r2, switchdist2, cutoff2, &U, &F);
          *alch_vdw_energy = myVdwLambda*U;
          *alch_vdw_force = myVdwLambda*F;
        }
      } else {
        // potential switching - also correct for no switching
        if (r2 <= Rmin2) {
          // normal LJ potential w/shift
          const BigReal epsilon = B*B/(4.0*A);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = U + (1 - myVdwLambda)*epsilon;
          *alch_vdw_force = F; 
        } else { // r2 > Rmin2
          // normal LJ potential with linear coupling
          vdw_switch(r2, switchdist2, cutoff2, switchfactor, &switchmul, \
              &switchmul2);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = myVdwLambda*switchmul*U;
          *alch_vdw_force = myVdwLambda*(switchmul*F + switchmul2*U);
        }
      }
    }
    if (myRepLambda2 < 1.0) {
      // Same as above, but energy only. This is done separately bc the
      // cutoffs are lambda dependent.
      const BigReal WCAshift_2 = Rmin2*(1 - myRepLambda2)*(1 - myRepLambda2);
      if (r2 <= Rmin2 - WCAshift_2) {
        const BigReal epsilon = B*B/(4.0*A);
        vdw_energy(A, B, r2 + WCAshift_2, &U_2);
        *alch_vdw_energy_2 = U_2 + epsilon;
      } else {
        *alch_vdw_energy_2 = 0.0;
      }
    } else { // myRepLambda2 == 1.0
      if (vdwForceSwitching) {
        // force and potential switching
        if (r2 <= Rmin2) {
          // normal LJ, potential w/two shifts
          const BigReal epsilon = B*B/(4.0*A);
          vdw_fswitch_shift(A, B, switchdist2, cutoff2, &dU_2);
          vdw_energy(A, B, r2, &U_2);
          *alch_vdw_energy_2 = \
              U_2 + (1 - myVdwLambda2)*epsilon + myVdwLambda2*dU_2;
        } else if (r2 <= switchdist2) {
          // normal LJ potential w/shift
          vdw_fswitch_shift(A, B, switchdist2, cutoff2, &dU_2);
          vdw_energy(A, B, r2, &U_2);
          *alch_vdw_energy_2 = myVdwLambda2*(U_2 + dU_2);
        } else { // r2 > switchdist
          vdw_fswitch_energy(A, B, r2, switchdist2, cutoff2, &U_2);
          *alch_vdw_energy_2 = myVdwLambda2*U_2;
        }
      } else {
        // potential switching - also correct for no switching
        if (r2 <= Rmin2) {
          // normal LJ potential w/shift
          const BigReal epsilon = B*B/(4.0*A);
          vdw_energy(A, B, r2, &U_2);
          *alch_vdw_energy_2 = U_2 + (1 - myVdwLambda2)*epsilon;
        } else { // r2 > Rmin2
          // normal LJ potential with linear coupling
          vdw_switch(r2, switchdist2, cutoff2, switchfactor, &switchmul, \
              &switchmul2);
          vdw_energy(A, B, r2, &U_2);
          *alch_vdw_energy_2 = myVdwLambda2*switchmul*U_2;
        }
      }
    }
  } else { // WCA-off
    if (vdwForceSwitching) {
      // force and potential switching
      if (r2 <= switchdist2) {
        // normal LJ potential w/shift
        vdw_forceandenergy(A, B, r2 + myVdwShift, &U, &F);
        vdw_fswitch_shift(A, B, switchdist2 + myVdwShift, cutoff2, &dU);
        vdw_energy(A, B, r2 + myVdwShift2, &U_2);
        vdw_fswitch_shift(A, B, switchdist2 + myVdwShift2, cutoff2, &dU_2);
        *alch_vdw_energy = myVdwLambda*(U + dU);
        *alch_vdw_energy_2 = myVdwLambda2*(U_2 + dU_2);
        *alch_vdw_force = myVdwLambda*F;
      } else { // r2 > switchdist2
        vdw_fswitch_forceandenergy(A, B, r2 + myVdwShift, \
            switchdist2 + myVdwShift, cutoff2, &U, &F);
        vdw_fswitch_energy(A, B, r2 + myVdwShift2, switchdist2 + myVdwShift2, \
            cutoff2, &U_2);
        *alch_vdw_energy = myVdwLambda*U;
        *alch_vdw_energy_2 = myVdwLambda2*U_2;
        *alch_vdw_force = myVdwLambda*F;
      }
    } else {
      // potential switching - also correct for no switching
      vdw_switch(r2, switchdist2, cutoff2, switchfactor, &switchmul, \
          &switchmul2);
      vdw_forceandenergy(A, B, r2 + myVdwShift, &U, &F);
      vdw_energy(A, B, r2 + myVdwShift2, &U_2); 
      *alch_vdw_energy = myVdwLambda*switchmul*U;
      *alch_vdw_energy_2 = myVdwLambda2*switchmul*U_2;
      *alch_vdw_force = myVdwLambda*(switchmul*F + switchmul2*U);
    }
  }
}

#define FEPFLAG
#define CALCENERGY

#define NBTYPE NBPAIR
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#define NBTYPE NBSELF
#include "ComputeNonbondedBase.h"
#define FULLELECT
#include "ComputeNonbondedBase.h"
#define MERGEELECT
#include "ComputeNonbondedBase.h"
#undef MERGEELECT
#define SLOWONLY
#include "ComputeNonbondedBase.h"
#undef SLOWONLY
#undef FULLELECT
#undef  NBTYPE

#undef CALCENERGY
#undef FEPFLAG


