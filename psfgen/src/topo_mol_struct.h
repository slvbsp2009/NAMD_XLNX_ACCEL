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
 *      $RCSfile: topo_mol_struct.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.26 $      $Date: 2020/03/10 04:54:55 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef TOPO_DEFS_MOL_H
#define TOPO_DEFS_MOL_H

#include "hasharray.h"
#include "memarena.h"
#include "topo_defs_struct.h"
#include "topo_mol.h"


#define NAMEMAXLEN 10 /* Max string length used for atom names and types */
#define NAMETOOLONG(X) ( strlen(X) >= NAMEMAXLEN )

#if defined(NEWPSFGEN)

/* define macro to identify islonepair atom's field to be used especially,
 * in the read-in of psf files with lonepairs and drude particles
 */
#define ISLONEPAIR 1
#define ISDRUDE 2

#endif

struct topo_mol_atom_t;

/* topo_mol_bond_t, topo_mol_angle_t, topo_mol_dihedral_t, topo_mol_improper_t
 * topo_mol_cmap_t, topo_mol_exclusion_t are fields of the topo_mol_atom_t,
 * meaning that every atom will have, for instance, two atoms in the 
 * topo_mol_bond_t. This seems to be a redundancy and should be worked out
 * carefully.
*/ 

/** Molecule bond data structure */
typedef struct topo_mol_bond_t {
  struct topo_mol_bond_t *next[2]; /* pointer to the next bonds in the linked list */
  struct topo_mol_atom_t *atom[2]; /* atoms belonging to the bond.  */
  int del;
} topo_mol_bond_t;

/** Molecule angle data structure */
typedef struct topo_mol_angle_t {
  struct topo_mol_angle_t *next[3]; /* pointer to the next angles in the linked list */
  struct topo_mol_atom_t *atom[3]; /* atoms belonging to the angles.  */
  int del;
} topo_mol_angle_t;

#if !defined(NEWPSFGEN)

/** Molecule dihedral data structure */
typedef struct topo_mol_dihedral_t {
  struct topo_mol_dihedral_t *next[4]; /* pointer to the next dihedral angles in the linked list */
  struct topo_mol_atom_t *atom[4]; /* atoms belonging to the dihedral angles */
  int del;
} topo_mol_dihedral_t;

#else

/** Molecule dihedral data structure */
typedef struct topo_mol_dihedral_t {
  struct topo_mol_dihedral_t *next; /* pointer to the next dihedral angles in the linked list */
  struct topo_mol_atom_t *atom[3]; /* atoms belonging to the dihedral angles */
  int del;
} topo_mol_dihedral_t;

#endif

#if !defined(NEWPSFGEN)
/** Molecule improper data structure */
typedef struct topo_mol_improper_t {
  struct topo_mol_improper_t *next[4]; /* pointer to the next improper angles in the linked list */
  struct topo_mol_atom_t *atom[4]; /* atoms belonging to the improper angles */
  int del;
} topo_mol_improper_t;
#else
typedef struct topo_mol_improper_t {
  struct topo_mol_improper_t *next; /* pointer to the next improper angles in the linked list */
  struct topo_mol_atom_t *atom[3]; /* atoms belonging to the improper angles */
  int del;
} topo_mol_improper_t;
#endif

/** Molecule CMAP data structure */
typedef struct topo_mol_cmap_t {
  struct topo_mol_cmap_t *next[8]; /* pointer to the next cmap in the linked list */
  struct topo_mol_atom_t *atom[8]; /* atoms belonging to the cmap */
  int del;
} topo_mol_cmap_t;

/** Molecule exclusion data structure */
typedef struct topo_mol_exclusion_t {
  struct topo_mol_exclusion_t *next[2]; /* pointer to the next exclusion in the linked list */
  struct topo_mol_atom_t *atom[2]; /* atoms belonging to the exclusion list */
  int del;
} topo_mol_exclusion_t;

/** Molecule Conformation data structure */
typedef struct topo_mol_conformation_t {
  struct topo_mol_conformation_t *next[4]; /* pointer to the next conformation in the linked list */
  struct topo_mol_atom_t *atom[4]; /* atoms belonging to the conformation list */
  int del;                         /* was the conformation deleted */
  int improper; /* flag defining an improper. Distinguished by the presence of a '*' just before
the name of the third atom. */
  double dist12;        /* distance between atoms 1 and 2 */
  double angle123;      /* angle composed by the atoms 1, 2 and 3 */
  double dihedral;      /* dihedral angle composed by the atoms all atoms in the IC card */
  double angle234;      /* angle composed by the atoms 2, 3 and 4 */
  double dist34;        /* distance between atoms 3 and 4 */
} topo_mol_conformation_t;


#define TOPO_MOL_XYZ_VOID 0 /* initializing value to assign to atoms->xyz_state */
#define TOPO_MOL_XYZ_SET 1 /* coordinates set from pdb or by command*/
#define TOPO_MOL_XYZ_GUESS 2 /* coordinates guessed */
#define TOPO_MOL_XYZ_BADGUESS 3 /* coordinates badly guessed and used to set occupancy 0 */


/** Molecule Atom data structure */
typedef struct topo_mol_atom_t {

#if !defined(NEWPSFGEN)
  struct topo_mol_atom_t *next; /* pointer to the next atom in the linked list */
#endif
  
  struct topo_mol_atom_t *copy;
  topo_mol_bond_t *bonds; /* all possible bonds associated with the atom, including itself */ 
  topo_mol_angle_t *angles; /* all possible angles associated with the atom, including itself */
  topo_mol_dihedral_t *dihedrals; /* all possible dihedrals angles associated with the atom, including itself */
  topo_mol_improper_t *impropers; /* all possible impropers angles associated with the atom, including itself */
  topo_mol_cmap_t *cmaps; /* the cmpa correction where the atom is defined, including itself */
  topo_mol_exclusion_t *exclusions; /* the exclusions where the atom is defined, including itself */
  topo_mol_conformation_t *conformations; /* the conformations where the atom is defined from the IC cards */
  char name[NAMEMAXLEN]; /* name of the atom */
  char type[NAMEMAXLEN]; /* type of the atom */
  char element[NAMEMAXLEN]; /* element of the atom */
  double mass;  /* mass of the atom */
  double charge; /* charge of the atom */
  double x,y,z; /* x, y and z coordinates of the atom */
  double vx,vy,vz; /* x, y and z componentes of the velocity vector of the atom */
  int xyz_state; /* how the coordinates were obtained? not set (TOPO_MOL_XYZ_VOID),
   set (TOPO_MOL_XYZ_SET) , guessed (TOPO_MOL_XYZ_GUESS) or badly guess (TOPO_MOL_XYZ_BADGUESS) */
  int partition; /* bfactor of the atom also used to flag the ones with coordinates guessed */
  int atomid; /* atomindex to be used to write in the files. Only at the writing time 
  are the ids attributed to the atoms */
  
#if defined(NEWPSFGEN)

  int del;/* signal that the atom was deleted */
  /* The lonepairs are identified as atoms with a non-null
   * lonepair pointer. topo_defs_lonepair completes the info of 
   * topo_defs_lonepair_t
   *
   * currently the maximum of 2 lonepairs per atom is implemented and is in 
   * drude forcefield
  */
  int isdrudlonepair; /* int flag to identify as lonepair or drude particle
                       * - for faster check.
                       */
  struct topo_mol_lonepair_t *lonepair; /* lone pairs 1 attached to the atom*/
  
  /* drude particles*/
  char dname[NAMEMAXLEN];/* drude particle name*/
  double alpha; /* particle polarizability */
  double thole; /* scale factor default = 1.3*/
  double dcharge; /* charge of the drude 
                   * KDRUDE is the force constant (in kcal/mol/Angst**2) for the 
                   * bond between. Default 500 kcal/mol/Angst**2
                   * CCELEC =332.071600 (Coulomb's constant CHARMM const value);
                   * 
                   * q = sqrt( 2*KDRUDE * alpha / CCELEC ) * sign(alpha)
                   */
  double *dxyz; /* pointer to the coordinates of the drude to be used when
                 * reading coordinates from a file already containing lone pairs
                 */
#endif
  
} topo_mol_atom_t;

#if defined(NEWPSFGEN)
/** Lonepair structure  - please see topo_defs_struct.h for more info*/
typedef struct topo_mol_lonepair_t {
  topo_mol_atom_t **atoms; /* atoms[0] == host (atom directly attached to the 
                            * lonepair)
                            * atoms[1] == atom defined after the host in the 
                            *    LONEPAIR entry
                            * atoms[2] == second atom defined in the LONEPAIR
                            * atoms[n+1] == NULL
                            */
  float distance;
  float angle; /* scale in case of the colinear*/
  float dihedral; /* ignored in case of the colinear*/
  int lptype; /* Colinear, Relatice, Bisector or Center as defined in the topo_defs.h
               * #define COLINEARLP 1 colinear lonepair
               * #define RELATIVELP 2 relative lonepair
               * #define BISECTORLP 3 bisector lonepair
               * #define CENTERLP 4 center lonepair
               */
} topo_mol_lonepair_t;

/** Anisotropy structure */
typedef struct topo_mol_anisotropy_t {
  struct topo_mol_atom_t **atoms; /* atoms involved in the definition of the anisotropy
                                   */
  struct topo_mol_anisotropy_t *next;
  float k11; /* definition of the anisotropic tensor 11 22 defined in the RTF
               * and the 33 is deducted from the previous 2
              */
  float k22;
  float k33;
  int del;
  int patch;
} topo_mol_anisotropy_t;
#endif

typedef struct topo_mol_residue_t {
  char resid[NAMEMAXLEN]; /* ID of the residue */
  char name[NAMEMAXLEN];  /* name of the residue */
  char chain[NAMEMAXLEN]; /* chain where the residue is located */

#if !defined(NEWPSFGEN)
  topo_mol_atom_t *atoms; /* linked list containing the atoms */
#else  
  topo_mol_atom_t **atomArray; /* array of atoms of the residue */
  int atomSize; /* Number of atoms in the residue. Different from the size of 
  array as the array size is  atomSize + 1(null pointer)*/
  int reordered; /* Was the original atoms' order changed by an insertion at the beginning 
  or in the middle? Patches tend to change the original order of the residues */
  int lonepairs; /* does the residue has lonepairs definined*/
  int numaniso; /* number of anisotropy definitions in the residue*/
  topo_mol_anisotropy_t *aniso; /* linked list with the anisotropy for drude ff*/
  char pres[NAMEMAXLEN];/* flag to store the patch name being applied to the residue*/
#endif

} topo_mol_residue_t;


typedef struct topo_mol_segment_t {
  char segid[NAMEMAXLEN]; /* segment ID */
  topo_mol_residue_t *residue_array; /* Residue array. Array index from hash table keys */
  hasharray *residue_hash; /* Residue hash table */

  int auto_angles; /* Flag for automatic determination of the angles */
  int auto_dihedrals; /* Flag for automatic determination of the dihedral angles */
  char pfirst[NAMEMAXLEN]; /* Patch to be applied to the first residue (Default NTER in proteins) */
  char plast[NAMEMAXLEN]; /* Patch to be applied to the last residue (Default CTER in proteins) */

} topo_mol_segment_t;

/** Patch structure - still need revision */
typedef struct topo_mol_patchres_t { 
  struct topo_mol_patchres_t *next; 
  char segid[NAMEMAXLEN];
  char resid[NAMEMAXLEN];
} topo_mol_patchres_t;

/** Patch structure  - still need revision*/
typedef struct topo_mol_patch_t {
  struct topo_mol_patch_t *next;
  char pname[NAMEMAXLEN];
  int npres;
  int deflt;
  topo_mol_patchres_t *patchresids;
} topo_mol_patch_t;

/** Main molecule structure */
struct topo_mol {
  void *newerror_handler_inter; /* tcl interpreter to out put message data */
  void *newerror_handler_vdata; /* Clientdata to out put message data (text) */
  void (*newerror_handler)(void *, void *, const char *); /* function to output messages. 
                                                   * This function is the
                                                   * newhandle_msg function defined in tcl_main.c 
                                                   * the first void is the Clientdata pointer,
                                                   * the second one is the tcl interpreter pointer
                                                  */
  
  topo_defs *defs; /* topology data */

  int npatch;
  topo_mol_patch_t *patches;
  topo_mol_patch_t *curpatch;

  topo_mol_segment_t **segment_array;  /* array of segments which indexes are 
                                        * defined from the hash table 
                                       */
  hasharray *segment_hash; /* hash table of the segment in the molecule */
  topo_mol_segment_t *buildseg; /* Current segment being built */

  memarena *arena;
  memarena *angle_arena;
  memarena *dihedral_arena;

#if defined(NEWPSFGEN)
  int drude; /* flag if the molecule contains drude particles */
  int lonepairs; /* flag the molecule contains lonepairs particles  */
#endif
};

topo_mol_bond_t * topo_mol_bond_next(
                topo_mol_bond_t *tuple, topo_mol_atom_t *atom);

topo_mol_angle_t * topo_mol_angle_next(
                topo_mol_angle_t *tuple, topo_mol_atom_t *atom);

topo_mol_dihedral_t * topo_mol_dihedral_next(
                topo_mol_dihedral_t *tuple, topo_mol_atom_t *atom);

topo_mol_improper_t * topo_mol_improper_next(
                topo_mol_improper_t *tuple, topo_mol_atom_t *atom);

topo_mol_cmap_t * topo_mol_cmap_next(
                topo_mol_cmap_t *tuple, topo_mol_atom_t *atom);

topo_mol_exclusion_t * topo_mol_exclusion_next(
                topo_mol_exclusion_t *tuple, topo_mol_atom_t *atom);

#endif

