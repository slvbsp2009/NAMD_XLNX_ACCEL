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
 *      $RCSfile: topo_defs_struct.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.13 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef TOPO_DEFS_STRUCT_H
#define TOPO_DEFS_STRUCT_H

#include "memarena.h"
#include "hasharray.h"
#include "topo_defs.h"

#define NAMEMAXLEN 10
#define NAMETOOLONG(X) ( strlen(X) >= NAMEMAXLEN )

typedef struct topo_defs_type_t {
  char name[NAMEMAXLEN];
  char element[NAMEMAXLEN];
  int id;
  double mass;
} topo_defs_type_t;

#if defined(NEWPSFGEN)  
/** LonePairs data structure */
/* Bare in mind that although different types of lone pairs are supported
 * (colinear, relative and center), they share the same basic information:
 * Distance from the host, the atoms associated with the lone pair
 * 2 in the case of colinear and 3 in the case of relative bisector and center
 * but not implemented yet - please update this upon implementation), 
 * angle (which in the colinear is the scale),
 * and dihedral (ignored in the colinear lone pair)
*/

typedef struct topo_defs_lonepair_t {
  // char name[NAMEMAXLEN]; /* Name of the lonepair. Duplicated info from atomdef,
  //                         * necessary to compare the names as the lone pairs
  //                         * have 2 definition sections : 
  //                         * 1 - in the listing of the atoms 
  //                         * 2 - in the definition of the lonepair type
  //                        */
  struct topo_defs_atom_t **atoms; /* atoms involved in the definition of the lone pair
                            * In this case I am using an array instead of defining
                            * the atoms individualy, like bonds, angles and dihedrals
                            * to protect against future implementations like Lonepair 
                            * Center that place a lone pair at the weighted center
                            * of atom selection without fixed limit of number of atoms.
                            * the application already in mentioned by Alex Mackarell is 
                            * to place a Lonepair in the center of rings, like beneze 
                           */
  /*atomsName and relAtoms serve to store information of lonepairs defined 
   * in patches belonging to different residues (rel != 0). The use of pointers
   * makes it faster to find the atoms, but this is not possible when the atoms
   * boling to different residues. In the drude forcefield the former case is 
   * the most common.
  */
  // struct topo_defs_lonepair_t *nextlonepair;
  double distance;
  double angle; /* scale in case of the colinear*/
  double dihedral; /* ignored in case of the colinear*/
  int lptype; /* type of lonepair, important to use the proper formula to place
             * the particle, values= COLINERLP, RELATIVELP, BISECTORLP or CENTERLP 
             * macros defined in topo_defs.h
            */
  int numatoms; /* number of atoms defining the lonepair*/
} topo_defs_lonepair_t;

typedef struct topo_defs_anisotropy_t {
  struct topo_defs_atom_t **atoms; /* atoms involved in the definition of the anisotropy
                                   */
  char **atomsname; /* atoms involved in the definition of the anisotropy
                               */
  int *res;
  int *rel;
  
  struct topo_defs_anisotropy_t *next;
  float k11;/* definition of the anisotropic tensor 11 22 defined in the RTF
               * and the 33 is dedeucted from the previous 2.
               * evaluate the constants K that go into the psf so they are 
               * evaluated once in the topology of the residue, and not
               * everytime a residue is defined in the molecule
              */
  float k22;
  float k33;
  int patch;
  int del;
} topo_defs_anisotropy_t;

#endif

typedef struct topo_defs_atom_t {
  struct topo_defs_atom_t *next;
  char name[NAMEMAXLEN];
  char type[NAMEMAXLEN];
  double charge;
  int res, rel; /* res - in case of being a patch applied to two different residues
                 * the declaration of the atom is precedeed by the relative residue
                 * the the atom belongs to. Example:
                 * PRES DISU
                 * ATOM 1CB  CT2   
                 * ATOM 1SG  SM
                 * GROUP                     
                 * ATOM 2SG  SM     
                 * ATOM 2CB  CT2
                 * rel - relative position of the atom: 
                 * '-'' -> previous residue
                 * '+' -> next residue
                 * '#' -> 2 residues after the current (source unkown)
                 */
  int del; /* Flag to sign the DELETE of the pacthes */
                  
#if defined(NEWPSFGEN)
int atomIndex; /* store the position of the atom in the topology and make the position
                * of the correspondent atom in the residue (in the molecule) to correspond
                * this index
                */
  int patch; /* int flag to identify atom as part of a patch as it is important
              * to get the right index
             */
  int islonepair; /* int flag to identify as lonepair - for faster check if the 
                   * particle is a lonepair.
                  */
  topo_defs_lonepair_t *lonepair; /* pointer to the lonepair data */
  
  /* drude particles*/
  char dname[NAMEMAXLEN];/* drude particle name*/
  double alpha; /* particle polarizability */
  double thole; /* scale factor default = 1.3*/
#endif
                  
} topo_defs_atom_t;

typedef struct topo_defs_bond_t {
  struct topo_defs_bond_t *next;

#if !defined(NEWPSFGEN)

  char atom1[NAMEMAXLEN];
  char atom2[NAMEMAXLEN];

#else
  
  char atomstr1[NAMEMAXLEN]; /* Store both atom name and the pointer to the 
                              * atom in cases it is not possible to find the
                              * atom, for example, in the declaration of 
                              * patches where sometimes the atoms are belong
                              * to different residues
                              */
  char atomstr2[NAMEMAXLEN];
  topo_defs_atom_t *atom1;
  topo_defs_atom_t *atom2;

#endif

  int res1, rel1;
  int res2, rel2;
  int del;
} topo_defs_bond_t;

typedef struct topo_defs_angle_t {
  struct topo_defs_angle_t *next;
  
#if !defined(NEWPSFGEN)

  char atom1[NAMEMAXLEN];
  char atom2[NAMEMAXLEN];
  char atom3[NAMEMAXLEN];

#else

  char atomstr1[NAMEMAXLEN];
  char atomstr2[NAMEMAXLEN];
  char atomstr3[NAMEMAXLEN];
  topo_defs_atom_t *atom1;
  topo_defs_atom_t *atom2;
  topo_defs_atom_t *atom3;

#endif

  int res1, rel1;
  int res2, rel2;
  int res3, rel3;
  int del;
} topo_defs_angle_t;

typedef struct topo_defs_dihedral_t {
  struct topo_defs_dihedral_t *next;
  
#if !defined(NEWPSFGEN)

  char atom1[NAMEMAXLEN];
  char atom2[NAMEMAXLEN];
  char atom3[NAMEMAXLEN];
  char atom4[NAMEMAXLEN];
  
#else

  char atomstr1[NAMEMAXLEN];
  char atomstr2[NAMEMAXLEN];
  char atomstr3[NAMEMAXLEN];
  char atomstr4[NAMEMAXLEN];
  topo_defs_atom_t *atom1;
  topo_defs_atom_t *atom2;
  topo_defs_atom_t *atom3;
  topo_defs_atom_t *atom4;

#endif
  
  int res1, rel1;
  int res2, rel2;
  int res3, rel3;
  int res4, rel4;
  int del;
} topo_defs_dihedral_t;

typedef struct topo_defs_improper_t {
  struct topo_defs_improper_t *next;
  
#if !defined(NEWPSFGEN)

  char atom1[NAMEMAXLEN];
  char atom2[NAMEMAXLEN];
  char atom3[NAMEMAXLEN];
  char atom4[NAMEMAXLEN];
  
#else

  char atomstr1[NAMEMAXLEN];
  char atomstr2[NAMEMAXLEN];
  char atomstr3[NAMEMAXLEN];
  char atomstr4[NAMEMAXLEN];
  topo_defs_atom_t *atom1;
  topo_defs_atom_t *atom2;
  topo_defs_atom_t *atom3;
  topo_defs_atom_t *atom4;

#endif
  
  int res1, rel1;
  int res2, rel2;
  int res3, rel3;
  int res4, rel4;
  int del;
} topo_defs_improper_t;

typedef struct topo_defs_cmap_t {
  struct topo_defs_cmap_t *next;
  char atoml[8][NAMEMAXLEN];
  int resl[8], rell[8];
  int del;
} topo_defs_cmap_t;

typedef struct topo_defs_exclusion_t {
  struct topo_defs_exclusion_t *next;

#if !defined(NEWPSFGEN)

  char atom1[NAMEMAXLEN];
  char atom2[NAMEMAXLEN];

#else

  char atomstr1[NAMEMAXLEN];
  char atomstr2[NAMEMAXLEN];
  topo_defs_atom_t *atom1;
  topo_defs_atom_t *atom2;

#endif

  int res1, rel1;
  int res2, rel2;
  int del;
} topo_defs_exclusion_t;

typedef struct topo_defs_conformation_t {
  struct topo_defs_conformation_t *next;

#if !defined(NEWPSFGEN)

  char atom1[NAMEMAXLEN];
  char atom2[NAMEMAXLEN];
  char atom3[NAMEMAXLEN];
  char atom4[NAMEMAXLEN];
  
#else

  char atomstr1[NAMEMAXLEN];
  char atomstr2[NAMEMAXLEN];
  char atomstr3[NAMEMAXLEN];
  char atomstr4[NAMEMAXLEN];
  topo_defs_atom_t *atom1;
  topo_defs_atom_t *atom2;
  topo_defs_atom_t *atom3;
  topo_defs_atom_t *atom4;

#endif

  int res1, rel1;
  int res2, rel2;
  int res3, rel3;
  int res4, rel4;
  int del;
  int improper;
  double dist12, angle123, dihedral, angle234, dist34;
} topo_defs_conformation_t;

typedef struct topo_defs_residue_t {
  char name[NAMEMAXLEN];
  int patch;
  topo_defs_atom_t *atoms;
  topo_defs_bond_t *bonds;
  topo_defs_angle_t *angles;
  topo_defs_dihedral_t *dihedrals;
  topo_defs_improper_t *impropers;
  topo_defs_cmap_t *cmaps;
  topo_defs_exclusion_t *exclusions;
  topo_defs_conformation_t *conformations;
  char pfirst[NAMEMAXLEN];
  char plast[NAMEMAXLEN];

#if defined(NEWPSFGEN)  
  int atomNum; //total atom number defined in the residue
  int lonepairs; /* int flag indicating if the residue contains lonepairs 
                 * used to signaling the need to update the "bonding" 
                 * information of the lonepair after it is definined in the molecule
                 * since the atoms are defined in the reverse order of the topology,
                 * the lp are defined before its hosts
                */
  topo_defs_anisotropy_t *aniso; /* anisotropy declarations in case of drude ff
                                  * Linked list instead of an array as the number is 
                                  * easy to predict
                                */
#endif
} topo_defs_residue_t;

typedef struct topo_defs_topofile_t {
/*   struct topo_defs_topofile_t *next; */
  char filename[256];
} topo_defs_topofile_t;

struct topo_defs {
  void *newerror_handler_inter; /* tcl interpreter to out put message data */
  void *newerror_handler_vdata; /* Clientdata to out put message data (text) */
  void (*newerror_handler)(void *, void *, const char *); /* function to output messages. 
                                                   * This function is the
                                                   * newhandle_msg function defined in tcl_main.c 
                                                   * the first void is the Clientdata pointer,
                                                   * the second one is the tcl interpreter pointer
                                                  */
  int auto_angles;
  int auto_dihedrals;
  int cmaps_present;
  char pfirst[NAMEMAXLEN];
  char plast[NAMEMAXLEN];

  topo_defs_topofile_t *topo_array;
  hasharray *topo_hash;

  topo_defs_type_t *type_array;
  hasharray *type_hash;

  topo_defs_residue_t *residue_array;
  hasharray *residue_hash;
  topo_defs_residue_t *buildres;
  int buildres_no_errors;
  memarena *arena;
};

#endif

