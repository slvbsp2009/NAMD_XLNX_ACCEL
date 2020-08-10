/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: topo_mol.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.15 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef TOPO_MOL_H
#define TOPO_MOL_H

#include "topo_defs.h"

struct topo_mol;
typedef struct topo_mol topo_mol;

topo_mol * topo_mol_create(topo_defs *defs);
void topo_mol_destroy(topo_mol *mol);

void topo_mol_error_handler(topo_mol *mol, void*, void *, void (*print_msg)(void *, void *,const char *));

int topo_mol_segment(topo_mol *mol, const char *segid);

int topo_mol_segment_first(topo_mol *mol, const char *rname);
int topo_mol_segment_last(topo_mol *mol, const char *rname);

int topo_mol_segment_auto_angles(topo_mol *mol, int autogen);
int topo_mol_segment_auto_dihedrals(topo_mol *mol, int autogen);

int topo_mol_residue(topo_mol *mol, const char *resid, const char *rname,
						const char *chain);
int topo_mol_mutate(topo_mol *mol, const char *resid, const char *rname);

int topo_mol_end(topo_mol *mol);

typedef struct topo_mol_ident_t {
  const char *segid;
  const char *resid;
  const char *aname;
} topo_mol_ident_t;

int topo_mol_patch(topo_mol *mol, const topo_mol_ident_t *targets,
			int ntargets, const char *rname, int prepend,
			int warn_angles, int warn_dihedrals, int deflt);

int topo_mol_regenerate_angles(topo_mol *mol);
int topo_mol_regenerate_dihedrals(topo_mol *mol);
int topo_mol_regenerate_resids(topo_mol *mol);

void topo_mol_delete_atom(topo_mol *mol, const topo_mol_ident_t *target);

int topo_mol_set_name(topo_mol *mol, const topo_mol_ident_t *target,
                           const char *name);

int topo_mol_set_resname(topo_mol *mol, const topo_mol_ident_t *target,
                              const char *rname);

int topo_mol_set_segid(topo_mol *mol, const topo_mol_ident_t *target,
                              const char *segid);

int topo_mol_multiply_atoms(topo_mol *mol, const topo_mol_ident_t *targets,
					int ntargets, int ncopies);

int topo_mol_set_element(topo_mol *mol, const topo_mol_ident_t *target,
					const char *element, int replace);

int topo_mol_set_chain(topo_mol *mol, const topo_mol_ident_t *target,
					const char *chain, int replace);

int topo_mol_set_xyz(topo_mol *mol, const topo_mol_ident_t *target,
					double x, double y, double z);

int topo_mol_set_vel(topo_mol *mol, const topo_mol_ident_t *target,
                                        double vx, double vy, double vz);

int topo_mol_set_mass(topo_mol *mol, const topo_mol_ident_t *target,
                      double mass);

int topo_mol_set_charge(topo_mol *mol, const topo_mol_ident_t *target,
                        double charge);

int topo_mol_set_bfactor(topo_mol *mol, const topo_mol_ident_t *target, 
                         double bfactor);

int topo_mol_guess_xyz(topo_mol *mol);

int topo_mol_add_patch(topo_mol *mol, const char *pname, int deflt);

int topo_mol_add_patchres(topo_mol *mol, const topo_mol_ident_t *targets);

int topo_mol_validate_patchres(topo_mol *mol, const char *pname, const char *segid, const char *resid);

#if defined(NEWPSFGEN)
struct topo_mol_atom_t;
int is_hydrogen(struct topo_mol_atom_t *atom);
int is_oxygen(struct topo_mol_atom_t *atom);

/* prototype of the function to assign coordinates of the drude particle
 * to the host atom
 */
int topo_mol_set_drude_xyz(topo_mol *mol, const topo_mol_ident_t *target,
                                        double x, double y, double z);
#endif

#endif

