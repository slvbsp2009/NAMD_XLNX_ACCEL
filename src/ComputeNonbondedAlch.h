/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 *  Inline Lennard-Jones energy/force/utility functions for alchemy          *
 *****************************************************************************/

/* Lennard-Jones force and energy.
 *
 * If using conventional soft-core, r2 should be passed as r2 + myVdwShift and
 * the final result should be scaled by myVdwLambda.
 */
inline void vdw_forceandenergy(const BigReal A, const BigReal B,
    const BigReal r2, BigReal *U, BigReal *F) {
  const BigReal inv_r2 = 1. / r2;
  const BigReal inv_r6 = inv_r2*inv_r2*inv_r2;
  *U = inv_r6*(A*inv_r6 - B);
  *F = 6*inv_r2*(2*(*U) + B*inv_r6);
}

// Same as above, energy only.
inline void vdw_energy(const BigReal A, const BigReal B,
    const BigReal r2, BigReal *U) {
  const BigReal inv_r2 = 1. / r2;
  const BigReal inv_r6 = inv_r2*inv_r2*inv_r2;
  *U = inv_r6*(A*inv_r6 - B);
}

/* Multipliers for the vdW potential switching function and its derivative.
 * This is correct even if switching is not active.
 */
inline void vdw_switch(const BigReal r2, const BigReal switchdist2,
    const BigReal cutoff2, const BigReal switchfactor, BigReal* switchmul,
    BigReal* switchmul2) {
  if (r2 <= switchdist2) {
    // "switching off" will always fall through here.
    *switchmul = 1.0;
    *switchmul2 = 0.0;
  } else {
    const BigReal cdiff2 = (cutoff2 - r2);
    const BigReal sdiff2 = (r2 - switchdist2);
    *switchmul = switchfactor*cdiff2*cdiff2*(cdiff2 + 3*sdiff2);
    *switchmul2 = 12*switchfactor*cdiff2*sdiff2;
  }
}

/* Force shifted Lennard-Jones force and energy.
 *
 * If using conventional soft-core, r2 and switchdist2 should be passed as
 * r2 + myVdwShift and switchdist2 + myVdwShift and the final result should be
 * scaled by myVdwLambda
 */
inline void vdw_fswitch_forceandenergy(const BigReal A, const BigReal B,
    const BigReal r2, const BigReal switchdist2, const BigReal cutoff2,
    BigReal* U, BigReal* F) {
  // TODO: Some rearrangement could be done to use rsqrt() instead?
  const BigReal inv_r2 = 1. / r2;
  const BigReal inv_r6 = inv_r2*inv_r2*inv_r2;
  const BigReal inv_r3 = sqrt(inv_r6);
  const BigReal inv_cutoff6 = 1. / (cutoff2*cutoff2*cutoff2);
  const BigReal inv_cutoff3 = sqrt(inv_cutoff6);
  const BigReal switchdist6 = switchdist2*switchdist2*switchdist2;
  const BigReal k_vdwa = A / (1 - switchdist6*inv_cutoff6);
  const BigReal k_vdwb = B / (1 - sqrt(switchdist6*inv_cutoff6));
  const BigReal tmpa = inv_r6 - inv_cutoff6;
  const BigReal tmpb = inv_r3 - inv_cutoff3;
  *U = k_vdwa*tmpa*tmpa - k_vdwb*tmpb*tmpb;
  *F = 6*inv_r2*(2*k_vdwa*tmpa*inv_r6 - k_vdwb*tmpb*inv_r3);
}

// Same as above, energy only.
inline void vdw_fswitch_energy(const BigReal A, const BigReal B,
    const BigReal r2, const BigReal switchdist2, const BigReal cutoff2,
    BigReal* U) {  
  // TODO: Some rearrangement could be done to use rsqrt() instead?
  const BigReal inv_r2 = 1. / r2;
  const BigReal inv_r6 = inv_r2*inv_r2*inv_r2;
  const BigReal inv_r3 = sqrt(inv_r6);
  const BigReal inv_cutoff6 = 1. / (cutoff2*cutoff2*cutoff2);
  const BigReal inv_cutoff3 = sqrt(inv_cutoff6);
  const BigReal switchdist6 = switchdist2*switchdist2*switchdist2;
  const BigReal k_vdwa = A / (1 - switchdist6*inv_cutoff6);
  const BigReal k_vdwb = B / (1 - sqrt(switchdist6*inv_cutoff6));
  const BigReal tmpa = inv_r6 - inv_cutoff6;
  const BigReal tmpb = inv_r3 - inv_cutoff3;
  *U = k_vdwa*tmpa*tmpa - k_vdwb*tmpb*tmpb;
}

/* Energy shifts for the vdW force/potential switching functions.
 *
 * If using conventional soft-core, switchdist2 should be passed as
 * switchdist2 + myVdwShift.
 */
inline void vdw_fswitch_shift(const BigReal A, const BigReal B,
    const BigReal switchdist2, const BigReal cutoff2, BigReal* dU) {
  // TODO: Some rearrangement could be done to use rsqrt() instead.
  const BigReal cutoff6 = cutoff2*cutoff2*cutoff2;
  const BigReal switchdist6 = switchdist2*switchdist2*switchdist2;
  const BigReal v_vdwa = -A / (cutoff6*switchdist6);
  const BigReal v_vdwb = -B / sqrt(cutoff6*switchdist6);
  *dU = v_vdwa - v_vdwb; //deltaV2 from Steinbach & Brooks
}

