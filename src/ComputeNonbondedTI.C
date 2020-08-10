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

/* ********************************** */
/* vdW energy, force and dU/dl for TI */
/* ********************************** */
inline void ti_vdw_force_energy_dUdl (BigReal A, BigReal B, BigReal r2, 
    BigReal myVdwShift, BigReal switchdist2, BigReal cutoff2,
    BigReal switchfactor, Bool vdwForceSwitching, BigReal myVdwLambda,
    BigReal alchVdwShiftCoeff, Bool alchWCAOn, BigReal myRepLambda,
    BigReal* alch_vdw_energy, BigReal* alch_vdw_force,
    BigReal* alch_vdw_dUdl) {
  BigReal U, F, dU, switchmul, switchmul2;
  /*
   * The variable naming here is unavoidably awful. A number immediately after
   * a variable implies a distance to that power. The suffix "_2" indicates an
   * evaluation at alchLambda2. The prefix "inv" means the inverse (previously
   * this also used underscores - it was awful).
   */
  if (alchWCAOn) {
    /* THIS BRANCH IS PARTIALLY INCOMPLETE AND BLOCKED INSIDE SimParameters
     *
     * The problem is that WCA creates TWO energy terms and thus two different
     * TI derivatives, but the reduction code is currently only set up for
     * one derivative. There is no such problem with FEP because the energy
     * is not decomposed. In principle, the contribution to the free energy
     * from each derivative is zero while the other is changing (e.g. while
     * repulsion is changing dispersion is not), but the derivatives are still
     * non-zero and it is NAMD convention to report them even as such.
     * However, even if we were willing to break this convention, both
     * derivatives still compute to the free energy at the boundary where
     * repulsion is completely on and dispersion is exactly zero - two
     * derivatives are thus required for a complete working implementation.
     *
     */

    // WCA-on, auxilliary, lambda-dependent cutoff based on Rmin
    //
    // Avoid divide by zero - correctly zeroes interaction below.
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
        *alch_vdw_dUdl = Rmin2*(1 - myRepLambda)*F;
      } else {
        *alch_vdw_energy = 0.0;
        *alch_vdw_force = 0.0;
        *alch_vdw_dUdl = 0.0;
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
          *alch_vdw_dUdl = dU - epsilon;
        } else if (r2 <= switchdist2) {
          // normal LJ potential w/shift
          vdw_fswitch_shift(A, B, switchdist2, cutoff2, &dU);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = myVdwLambda*(U + dU);
          *alch_vdw_force = myVdwLambda*F;
          *alch_vdw_dUdl = U + dU;
        } else { // r2 > switchdist
          // normal LJ potential with linear coupling
          vdw_fswitch_forceandenergy(A, B, r2, switchdist2, cutoff2, &U, &F);
          *alch_vdw_energy = myVdwLambda*U;
          *alch_vdw_force = myVdwLambda*F;
          *alch_vdw_dUdl = U;
        }
      } else {
        // potential switching - also correct for no switching
        if (r2 <= Rmin2) {
          // normal LJ potential w/shift
          const BigReal epsilon = B*B/(4.0*A);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = U + (1 - myVdwLambda)*epsilon;
          *alch_vdw_force = F;
          *alch_vdw_dUdl = -epsilon;
        } else { // r2 > Rmin2
          // normal LJ potential with linear coupling
          vdw_switch(r2, switchdist2, cutoff2, switchfactor, &switchmul, \
              &switchmul2);
          vdw_forceandenergy(A, B, r2, &U, &F);
          *alch_vdw_energy = myVdwLambda*switchmul*U;
          *alch_vdw_force = myVdwLambda*(switchmul*F + switchmul2*U);
          *alch_vdw_dUdl = switchmul*U;
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
        *alch_vdw_energy = myVdwLambda*(U + dU);
        *alch_vdw_force = myVdwLambda*F;
        *alch_vdw_dUdl = U + 0.5*myVdwLambda*alchVdwShiftCoeff*F + dU;
      } else { // r2 > switchdist2
        vdw_fswitch_forceandenergy(A, B, r2 + myVdwShift, \
            switchdist2 + myVdwShift, cutoff2, &U, &F);
        *alch_vdw_energy = myVdwLambda*U;
        *alch_vdw_force = myVdwLambda*F;
        *alch_vdw_dUdl = U + 0.5*myVdwLambda*alchVdwShiftCoeff*F;
      }
    } else {
      // potential switching - also correct for no switching
      vdw_switch(r2, switchdist2, cutoff2, switchfactor, &switchmul, \
          &switchmul2);
      vdw_forceandenergy(A, B, r2 + myVdwShift, &U, &F);
      *alch_vdw_energy = myVdwLambda*switchmul*U;
      *alch_vdw_force = myVdwLambda*(switchmul*F + switchmul2*U);
      *alch_vdw_dUdl = switchmul*(U + 0.5*myVdwLambda*alchVdwShiftCoeff*F);
    }
  }
}


/*************THERMODYNAMIC INTEGRATION*************/
#define TIFLAG
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

#undef TIFLAG

