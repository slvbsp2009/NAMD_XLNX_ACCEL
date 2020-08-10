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
 *      $RCSfile: topo_defs.c,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.21 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "topo_defs_struct.h"

#if defined(_MSC_VER)
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif

#if defined(NEWPSFGEN)
#include "psfgen.h"
static topo_defs_atom_t * topo_defs_find_atom_by_name(topo_defs *defs, 
                            const char *aname, const int res, const int rel) {
  topo_defs_atom_t *atom = 0;
  
  // return 0 if the atom refers to other residue (rel != 0)
  if (rel)
    return 0;  
  if ( ! defs->buildres || ! defs->buildres->atoms )
    return 0;
  for ( atom = defs->buildres->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(aname,atom->name) && res == atom->res) break;
  }
  return atom;
}

#endif

topo_defs * topo_defs_create(void) {
  topo_defs *defs;
  if ( (defs = (topo_defs*) malloc(sizeof(topo_defs))) ) {
    defs->newerror_handler_inter = 0;
    defs->newerror_handler_vdata = 0;
    defs->newerror_handler = 0;
    defs->auto_angles = 0;
    defs->auto_dihedrals = 0;
    defs->cmaps_present = 0;
    strcpy(defs->pfirst,"");
    strcpy(defs->plast,"");
    defs->buildres = 0;
    defs->buildres_no_errors = 0;
    defs->topo_hash = hasharray_create(
	(void**) &(defs->topo_array), sizeof(topo_defs_topofile_t));
    defs->type_hash = hasharray_create(
	(void**) &(defs->type_array), sizeof(topo_defs_type_t));
    defs->residue_hash = hasharray_create(
	(void**) &(defs->residue_array), sizeof(topo_defs_residue_t));
    defs->arena = memarena_create();
    if ( ! defs->type_hash || ! defs->residue_hash ||
	! defs->arena || ! defs->topo_hash ||
	topo_defs_residue(defs,"NONE",1) ||
	topo_defs_residue(defs,"None",1) ||
	topo_defs_residue(defs,"none",1) ) {
      topo_defs_destroy(defs);
      return 0;
    }
    topo_defs_end(defs);
  }
  return defs;
}

void topo_defs_destroy(topo_defs *defs) {
  int i,n;
  struct topo_defs_atom_t *a, *a2;
  struct topo_defs_bond_t *b, *b2;
  struct topo_defs_angle_t *an, *an2;
  struct topo_defs_dihedral_t *di, *di2;
  struct topo_defs_improper_t *im, *im2;
  struct topo_defs_cmap_t *cm, *cm2;
  struct topo_defs_exclusion_t *ex, *ex2;
  struct topo_defs_conformation_t *c, *c2;
  
  if ( ! defs ) return;
  hasharray_destroy(defs->topo_hash);
  hasharray_destroy(defs->type_hash);
  n = hasharray_count(defs->residue_hash);
  for ( i=0; i<n; ++i ) {
    a = defs->residue_array[i].atoms;
    while ( a ) {
      a2 = a->next;
      free((void*)a);
      a = a2;
    }
    b = defs->residue_array[i].bonds;
    while ( b ) {
      b2 = b->next;
      free((void*)b);
      b = b2;
    }
    an = defs->residue_array[i].angles;
    while ( an ) {
      an2 = an->next;
      free((void*)an);
      an = an2;
    }
    di = defs->residue_array[i].dihedrals;
    while ( di ) {
      di2 = di->next;
      free((void*)di);
      di = di2;
    }
    im = defs->residue_array[i].impropers;
    while ( im ) {
      im2 = im->next;
      free((void*)im);
      im = im2;
    }
    cm = defs->residue_array[i].cmaps;
    while ( cm ) {
      cm2 = cm->next;
      free((void*)cm);
      cm = cm2;
    }
    ex = defs->residue_array[i].exclusions;
    while ( ex ) {
      ex2 = ex->next;
      free((void*)ex);
      ex = ex2;
    }
    c = defs->residue_array[i].conformations;
    while ( c ) {
      c2 = c->next;
      free((void*)c);
      c = c2;
    }
  }
  hasharray_destroy(defs->residue_hash);
  memarena_destroy(defs->arena);
  free((void*)defs);
}

void topo_defs_error_handler(topo_defs *defs, void *vdata, void *v, 
          void (*print_msg)(void *, void *, const char *)) {
  if ( defs ) {
    defs->newerror_handler = print_msg;
    defs->newerror_handler_inter = v;
    defs->newerror_handler_vdata = vdata;
  }
}

/* internal method */
void topo_defs_log_error(topo_defs *defs, const char *msg) {
  if (defs && msg && defs->newerror_handler) {
    defs->newerror_handler(defs->newerror_handler_vdata, defs->newerror_handler_inter, msg);
  }
}

void topo_defs_auto_angles(topo_defs *defs, int autogen) {
  if ( defs ) defs->auto_angles = ! ! autogen;
}

void topo_defs_auto_dihedrals(topo_defs *defs, int autogen) {
  if ( defs ) defs->auto_dihedrals = ! ! autogen;
}

int topo_defs_type(topo_defs *defs, const char *atype, const char *element, double mass, int id) {
  int i;
  topo_defs_type_t *newitem;
  char errmsg[64 + NAMEMAXLEN];
  if ( ! defs ) return -1;
  if ( NAMETOOLONG(atype) ) return -2;
  if ( NAMETOOLONG(element) ) return -3;
  if ( ( i = hasharray_index(defs->type_hash,atype) ) != HASHARRAY_FAIL ) {
    sprintf(errmsg,"duplicate type key %s",atype);
    topo_defs_log_error(defs,errmsg);
    newitem = &defs->type_array[i];
  } else {
    i = hasharray_insert(defs->type_hash,atype);
    if ( i == HASHARRAY_FAIL ) return -4;
    newitem = &defs->type_array[i];
    strcpy(newitem->name,atype);
  }
  newitem->id = id;
  strcpy(newitem->element,element);
  newitem->mass = mass;
  return 0;
}

int topo_defs_residue(topo_defs *defs, const char *rname, int patch) {
  int i;
  topo_defs_residue_t *newitem;
  char errmsg[64 + NAMEMAXLEN];
  if ( ! defs ) return -1;
  defs->buildres = 0;
  defs->buildres_no_errors = 0;
  if ( NAMETOOLONG(rname) ) return -2;
  if ( ( i = hasharray_index(defs->residue_hash,rname) ) != HASHARRAY_FAIL ) {
    char *oldname = defs->residue_array[i].name;
    if ( strcmp(rname,oldname) ) {
      sprintf(errmsg,"replacing residue alias %s for %s with new residue %s",rname,oldname,rname);
      topo_defs_log_error(defs,errmsg);
      hasharray_delete(defs->residue_hash,rname);
    } else {
      sprintf(errmsg,"duplicate residue key %s will be ignored",rname);
      topo_defs_log_error(defs,errmsg);
      /* newitem = &defs->residue_array[i]; */
      defs->buildres_no_errors = 1;
      return 0;
    }
  }
  i = hasharray_insert(defs->residue_hash,rname);
  if ( i == HASHARRAY_FAIL ) return -4;
  newitem = &defs->residue_array[i];
  strcpy(newitem->name,rname);
  newitem->patch = patch;
  newitem->atoms = 0;
  newitem->bonds = 0;
  newitem->angles = 0;
  newitem->dihedrals = 0;
  newitem->impropers = 0;
  newitem->cmaps = 0;
  newitem->exclusions = 0;
  newitem->conformations = 0;
  
#if defined(NEWPSFGEN)

  newitem->atomNum = 0;
  newitem->lonepairs = 0;
  newitem->aniso = 0;
  
#endif

  strcpy(newitem->pfirst,defs->pfirst);
  strcpy(newitem->plast,defs->plast);
  defs->buildres = newitem;
  return 0;
}

int topo_defs_end(topo_defs *defs) {
  if ( ! defs ) return -1;
  defs->buildres = 0;
  defs->buildres_no_errors = 0;
  return 0;
}

/** Define atoms from the topology. This function is also called when defining
  * atoms in the pactch definition.
  */
#if !defined(NEWPSFGEN)

int topo_defs_atom(topo_defs *defs, const char *rname, int del,
	const char *aname, int ares, int arel,
	const char *atype, double charge) {
    
#else

int topo_defs_atom(topo_defs *defs, const char *rname, int del,
	const char *aname, int ares, int arel,
	const char *atype, double charge, const char *dname,
  double alpha, double thole) {
    
  char tmpdname[NAMEMAXLEN];
#endif

  topo_defs_atom_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for atom");
    return -1;
  }
  if ( NAMETOOLONG(aname) ) return -2;
  if ( NAMETOOLONG(atype) ) return -3;
  if ( ares && ! defs->buildres->patch ) return -4;
  if ( arel && ! defs->buildres->patch ) return -4;
  if ( del && ! defs->buildres->patch ) return -5;
  newitem = (topo_defs_atom_t*) malloc(sizeof(topo_defs_atom_t));
  if ( ! newitem )  return -6;
  newitem->res = ares;
  newitem->rel = arel;
  newitem->del = del;
  newitem->charge = charge;
  strcpy(newitem->name,aname);
  strcpy(newitem->type,atype);


#if defined(NEWPSFGEN)
  
  /* Flag the atom as patch */ 
  if (defs->buildres->patch) {
    newitem->patch = 1;
  } else {
    newitem->patch = 0;
  }
  newitem->atomIndex = defs->buildres->atomNum;
  
  /* only increment the atom index if the this does not refers
   * to the declaration of the delete command in the patch definition
   */
  if (!del) {
    ++defs->buildres->atomNum;
  }
  
  /* defining lone pair */
  newitem->lonepair = 0;
  if (! strncasecmp("LP", atype, 2) ) {

    newitem->islonepair = 1;
  } else {

    newitem->islonepair = 0;
  }
  /* alpha identifies the atom as a drude particle host */
  if (alpha) {
    if (!dname) {
      strcpy(tmpdname,"D");
      strcat(tmpdname,aname);
      strcpy(newitem->dname,tmpdname);
    } else {
      strcpy(newitem->dname,dname);
    }
    newitem->alpha = alpha;
    newitem->thole = thole;
  } else {
    newitem->dname[0] = '\0';
    newitem->alpha = 0;
    newitem->thole = 0;
  }
#endif

newitem->next = defs->buildres->atoms;
defs->buildres->atoms = newitem;

  return 0;
}

#if defined(NEWPSFGEN)

/** Update the information of lone pairs */
/* The lonepairs are identified as atoms with a non-null
 * lonepair pointer. topo_defs_lonepair complets the info of 
 * topo_defs_lonepair_t
*/

int topo_defs_lonepair(topo_defs *defs, char **atomsname, int numatoms, 
  double distance, double angle, double dihedral, int lptype, int *res) {
  
  topo_defs_atom_t *atomlppointer = 0; /* pointer to the atom defined as lonepair*/
  topo_defs_lonepair_t *lonepair =0;  /* pointer to the lonepair data of the atom*/
  char errmsg[128];
  int i;
  
  /* Cannot use topo_defs_find_atom_by_name for patches as the 
   * command returns 0 for atoms from different residues (rel != 0).
   * But in patches like DISU, we need to search using res to find 
   * the lonepair correspondent to the correct definition.
  */
  atomlppointer = topo_defs_find_atom_by_name(defs,atomsname[0],res[0],0);
  if ( !atomlppointer) {
    sprintf(errmsg, "could not find / %s", atomsname[0]);
    topo_defs_log_error(defs,errmsg);
    return -1;
  }
  
  atomlppointer->lonepair = (topo_defs_lonepair_t*) malloc(sizeof(topo_defs_lonepair_t));
  if ( ! atomlppointer->lonepair )  return -2;
  
  lonepair = atomlppointer->lonepair;
  
  lonepair->atoms = (topo_defs_atom_t **)malloc((numatoms)*sizeof(topo_defs_atom_t*));
  for (i = 0; i < numatoms; i++) {
    if ( NAMETOOLONG(atomsname[i]) ) return -3;

    lonepair->atoms[i] = topo_defs_find_atom_by_name(defs,atomsname[i],res[i],0);
    if ( !lonepair->atoms[i] ) {
      sprintf(errmsg, "atom %s not found in lonepair %s", atomsname[i], 
      atomsname[0]);
      topo_defs_log_error(defs,errmsg);
      return -3;
    }
  }
  lonepair->distance = distance;
  /* In the colinear lp, angle == scale*/
  if (lptype != COLINEARLP) {
    angle *= (M_PI/180.0);
  }
  lonepair->angle = angle;
  lonepair->dihedral = dihedral * (M_PI/180.0);
  lonepair->lptype = lptype;
  lonepair->numatoms = numatoms;
  
  defs->buildres->lonepairs = 1;
  return 0;
}

int topo_defs_anisotropy(topo_defs *defs, char **atomsname, int *res, int *rel,
  double a11, double a22, int del) {
  
  topo_defs_anisotropy_t *newitem = NULL;

  int i, diff;
  
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for anisotropy");
    return -1;
  }
  for (i = 0; i < 4; i++) {
     if ( NAMETOOLONG(atomsname[i]) ) return -2;
     if ( res[i] && ! defs->buildres->patch ) return -4;
  }
  
  if ( del && ! defs->buildres->patch ) return -5;
  
  /* Check if this anisotropy was defined already*/
  
  for ( newitem = defs->buildres->aniso; newitem; newitem=newitem->next ) {
    char errmsg[128];
    diff = 0;
    for (i = 0; i < 4; i++) {
      if ( strcmp(newitem->atomsname[i],atomsname[i]) ) {
        diff = 1; 
        break;
      }
      if ( newitem->res[i] != res[i] ) {
        diff = 1; 
        break;
      }
      if ( newitem->rel[i] != rel[i] ) {
        diff = 1; 
        break;
      }
    }
    if (diff) continue; 
    if ( newitem->del != del ) continue;
    
    sprintf(errmsg, "duplicate anisotropy %s %s %s %s in residue %s", atomsname[0], 
          atomsname[1], atomsname[2], atomsname[3], defs->buildres->name);
    topo_defs_log_error(defs,errmsg);
    return -6;
  }
  
  newitem = (topo_defs_anisotropy_t*)malloc(sizeof(topo_defs_anisotropy_t));
  if ( ! newitem )  return -7;
  /* Define anisotropy */
  newitem->atoms = (topo_defs_atom_t**)malloc(4*sizeof(topo_defs_atom_t*));
  newitem->atomsname = (char**)malloc(4*sizeof(char*));
  newitem->res = (int*)malloc(4*sizeof(int));
  newitem->rel = (int*)malloc(4*sizeof(int));
  for (i = 0; i < 4; i++) {    
    
    newitem->atoms[i] = topo_defs_find_atom_by_name(defs,atomsname[i], res[i],rel[i]);
    newitem->atomsname[i] = (char*)malloc(sizeof(atomsname[i]));
    strcpy(newitem->atomsname[i], atomsname[i]);
    newitem->res[i] = res[i];
    newitem->rel[i] = rel[i];

  }
  newitem->del = 0;
  /* From the CHARMM code (genpsf.src)
   * translate:
   *     A11ANIS(IPT) == a11
   *    A22ANIS(IPT) == a22
   * A33 = THREE-A11ANIS(IPT)-A22ANIS(IPT)
   * K11(NANISO) = ONE/A11ANIS(IPT)
   * K22(NANISO) = ONE/A22ANIS(IPT)
   * K33(NANISO) = ONE/A33
   * K11(IANISO) = KDRUDE*K11(IANISO)
   * K22(IANISO) = KDRUDE*K22(IANISO)
   * K33(IANISO) = KDRUDE*K33(IANISO)
   * K33(IANISO) = K33(IANISO) - KDRUDE
   * K11(IANISO) = K11(IANISO) - KDRUDE - K33(IANISO)
   * K22(IANISO) = K22(IANISO) - KDRUDE - K33(IANISO)
   * KDRUDE is the force constant (in kcal/mol/Angst**2) for the bond
   * between. Default 500 kcal/mol/Angst**2
  */
  
  newitem->k33 = (K_DRUDE/(3 - a11 - a22));
  newitem->k33 -= K_DRUDE;
  newitem->k11 = (K_DRUDE/a11) - K_DRUDE - newitem->k33;
  newitem->k22 = (K_DRUDE/a22) - K_DRUDE - newitem->k33;
  newitem->patch = defs->buildres->patch;
  
  newitem->next = defs->buildres->aniso;
  defs->buildres->aniso = newitem;
  return 0;
}
#endif
int topo_defs_bond(topo_defs *defs, const char *rname, int del,
	const char *a1name, int a1res, int a1rel,
	const char *a2name, int a2res, int a2rel) {
  topo_defs_bond_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for bond");
    return -1;
  }
  if ( NAMETOOLONG(a1name) ) return -2;
  if ( NAMETOOLONG(a2name) ) return -3;
  if ( del && ! defs->buildres->patch ) return -4;
  if ( ( a1res || a2res ) && ! defs->buildres->patch ) return -4;
  for ( newitem=defs->buildres->bonds; newitem; newitem=newitem->next ) {
    char errmsg[128];
    if ( newitem->res1 != a1res ) continue;
    if ( newitem->rel1 != a1rel ) continue;
    if ( newitem->res2 != a2res ) continue;
    if ( newitem->rel2 != a2rel ) continue;
    if ( newitem->del != del ) continue;

#if !defined(NEWPSFGEN)  
    if ( strcmp(newitem->atom1,a1name) ) continue;
    if ( strcmp(newitem->atom2,a2name) ) continue;
#else
    if ( strcmp(newitem->atomstr1,a1name) ) continue;
    if ( strcmp(newitem->atomstr2,a2name) ) continue;
#endif
    
    sprintf(errmsg, "duplicate bond %s %s in residue %s", a1name, a2name, defs->buildres->name);
    topo_defs_log_error(defs,errmsg);
    return -6; /* XXX What is -6? replace with symbolic constant */
  }

  for ( newitem=defs->buildres->bonds; newitem; newitem=newitem->next ) {
    char errmsg[128];
    if ( newitem->res1 != a2res ) continue;
    if ( newitem->rel1 != a2rel ) continue;
    if ( newitem->res2 != a1res ) continue;
    if ( newitem->rel2 != a1rel ) continue;
    if ( newitem->del != del ) continue;
    
#if !defined(NEWPSFGEN)  
    if ( strcmp(newitem->atom1,a2name) ) continue;
    if ( strcmp(newitem->atom2,a1name) ) continue;
#else
    if ( strcmp(newitem->atomstr1,a2name) ) continue;
    if ( strcmp(newitem->atomstr2,a1name) ) continue;
#endif
    
    sprintf(errmsg, "duplicate bond %s %s in residue %s", a1name, a2name, defs->buildres->name);
    topo_defs_log_error(defs,errmsg);
    return -7;
  }

  if ( ( a1res == a2res ) && ( a1rel == a2rel ) && ! strcmp(a1name,a2name) ) {
    char errmsg[128];
    sprintf(errmsg, "self bond %s %s in residue %s", a1name, a2name, defs->buildres->name);
    topo_defs_log_error(defs,errmsg);
    return -8;
  }

  newitem = (topo_defs_bond_t*) malloc(sizeof(topo_defs_bond_t));
  if ( ! newitem )  return -5;
  newitem->res1 = a1res;
  newitem->rel1 = a1rel;
  newitem->res2 = a2res;
  newitem->rel2 = a2rel;
  newitem->del = del;

#if !defined(NEWPSFGEN)  
  strcpy(newitem->atom1,a1name);
  strcpy(newitem->atom2,a2name);
#else
  newitem->atom1 = topo_defs_find_atom_by_name(defs, a1name, a1res, a1rel);
  newitem->atom2 = topo_defs_find_atom_by_name(defs, a2name, a2res, a2rel);
  strcpy(newitem->atomstr1, a1name);
  strcpy(newitem->atomstr2, a2name);
#endif

  newitem->next = defs->buildres->bonds;
  defs->buildres->bonds = newitem;
  return 0;
}

int topo_defs_angle(topo_defs *defs, const char *rname, int del,
	const char *a1name, int a1res, int a1rel,
	const char *a2name, int a2res, int a2rel,
	const char *a3name, int a3res, int a3rel) {
  topo_defs_angle_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for angle");
    return -1;
  }
  if ( NAMETOOLONG(a1name) ) return -2;
  if ( NAMETOOLONG(a2name) ) return -3;
  if ( NAMETOOLONG(a3name) ) return -4;
  if ( del && ! defs->buildres->patch ) return -5;
  if ( ( a1res || a2res || a3res ) && ! defs->buildres->patch ) return -6;
  newitem = (topo_defs_angle_t*) malloc(sizeof(topo_defs_angle_t));
  if ( ! newitem )  return -7;
  newitem->res1 = a1res;
  newitem->rel1 = a1rel;
  newitem->res2 = a2res;
  newitem->rel2 = a2rel;
  newitem->res3 = a3res;
  newitem->rel3 = a3rel;
  newitem->del = del;

#if !defined(NEWPSFGEN)  

  strcpy(newitem->atom1,a1name);
  strcpy(newitem->atom2,a2name);
  strcpy(newitem->atom3,a3name);

#else

  newitem->atom1 = topo_defs_find_atom_by_name(defs, a1name, a1res, a1rel);
  newitem->atom2 = topo_defs_find_atom_by_name(defs, a2name, a2res, a2rel);
  newitem->atom3 = topo_defs_find_atom_by_name(defs, a3name, a3res, a3rel);
  strcpy(newitem->atomstr1,a1name);
  strcpy(newitem->atomstr2,a2name);
  strcpy(newitem->atomstr3,a3name);

#endif

  newitem->next = defs->buildres->angles;
  defs->buildres->angles = newitem;
  return 0;
}

int topo_defs_dihedral(topo_defs *defs, const char *rname, int del,
	const char *a1name, int a1res, int a1rel,
	const char *a2name, int a2res, int a2rel,
	const char *a3name, int a3res, int a3rel,
	const char *a4name, int a4res, int a4rel) {
  topo_defs_dihedral_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for dihedral");
    return -1;
  }
  if ( NAMETOOLONG(a1name) ) return -2;
  if ( NAMETOOLONG(a2name) ) return -3;
  if ( NAMETOOLONG(a3name) ) return -4;
  if ( NAMETOOLONG(a4name) ) return -5;
  if ( del && ! defs->buildres->patch ) return -6;
  if ( ( a1res || a2res || a3res || a4res ) &&
			! defs->buildres->patch ) return -7;
  newitem = (topo_defs_dihedral_t*) malloc(sizeof(topo_defs_dihedral_t));
  if ( ! newitem )  return -8;
  newitem->res1 = a1res;
  newitem->rel1 = a1rel;
  newitem->res2 = a2res;
  newitem->rel2 = a2rel;
  newitem->res3 = a3res;
  newitem->rel3 = a3rel;
  newitem->res4 = a4res;
  newitem->rel4 = a4rel;
  newitem->del = del;
  
#if !defined(NEWPSFGEN)  

  strcpy(newitem->atom1,a1name);
  strcpy(newitem->atom2,a2name);
  strcpy(newitem->atom3,a3name);
  strcpy(newitem->atom4,a4name);

#else

  newitem->atom1 = topo_defs_find_atom_by_name(defs, a1name, a1res, a1rel);
  newitem->atom2 = topo_defs_find_atom_by_name(defs, a2name, a2res, a2rel);
  newitem->atom3 = topo_defs_find_atom_by_name(defs, a3name, a3res, a3rel);
  newitem->atom4 = topo_defs_find_atom_by_name(defs, a4name, a4res, a4rel);
  strcpy(newitem->atomstr1,a1name);
  strcpy(newitem->atomstr2,a2name);
  strcpy(newitem->atomstr3,a3name);
  strcpy(newitem->atomstr4,a4name);

#endif

  newitem->next = defs->buildres->dihedrals;
  defs->buildres->dihedrals = newitem;
  return 0;
}

int topo_defs_improper(topo_defs *defs, const char *rname, int del,
	const char *a1name, int a1res, int a1rel,
	const char *a2name, int a2res, int a2rel,
	const char *a3name, int a3res, int a3rel,
	const char *a4name, int a4res, int a4rel) {
  topo_defs_improper_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for improper");
    return -1;
  }
  if ( NAMETOOLONG(a1name) ) return -2;
  if ( NAMETOOLONG(a2name) ) return -3;
  if ( NAMETOOLONG(a3name) ) return -4;
  if ( NAMETOOLONG(a4name) ) return -5;
  if ( del && ! defs->buildres->patch ) return -6;
  if ( ( a1res || a2res || a3res || a4res ) &&
			! defs->buildres->patch ) return -7;
  newitem = (topo_defs_improper_t*) malloc(sizeof(topo_defs_improper_t));
  if ( ! newitem )  return -8;
  newitem->res1 = a1res;
  newitem->rel1 = a1rel;
  newitem->res2 = a2res;
  newitem->rel2 = a2rel;
  newitem->res3 = a3res;
  newitem->rel3 = a3rel;
  newitem->res4 = a4res;
  newitem->rel4 = a4rel;
  newitem->del = del;

#if !defined(NEWPSFGEN)  

  strcpy(newitem->atom1,a1name);
  strcpy(newitem->atom2,a2name);
  strcpy(newitem->atom3,a3name);
  strcpy(newitem->atom4,a4name);

#else

  newitem->atom1 = topo_defs_find_atom_by_name(defs, a1name, a1res, a1rel);
  newitem->atom2 = topo_defs_find_atom_by_name(defs, a2name, a2res, a2rel);
  newitem->atom3 = topo_defs_find_atom_by_name(defs, a3name, a3res, a3rel);
  newitem->atom4 = topo_defs_find_atom_by_name(defs, a4name, a4res, a4rel);
  strcpy(newitem->atomstr1,a1name);
  strcpy(newitem->atomstr2,a2name);
  strcpy(newitem->atomstr3,a3name);
  strcpy(newitem->atomstr4,a4name);

#endif
  
  newitem->next = defs->buildres->impropers;
  defs->buildres->impropers = newitem;
  return 0;
}

int topo_defs_cmap(topo_defs *defs, const char *rname, int del,
	const char* const anamel[8], const int aresl[8], const int arell[8]) {
  int i;
  topo_defs_cmap_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for cmap");
    return -1;
  }
  for ( i=0; i<8; ++i ) {
    if ( NAMETOOLONG(anamel[i]) ) return -2-i;
  }
  if ( del && ! defs->buildres->patch ) return -10;
  if ( ( aresl[0] || aresl[1] || aresl[2] || aresl[3] ||
         aresl[4] || aresl[5] || aresl[6] || aresl[7] ) &&
			! defs->buildres->patch ) return -11;
  newitem = (topo_defs_cmap_t*) malloc(sizeof(topo_defs_cmap_t));
  if ( ! newitem )  return -12;
  for ( i=0; i<8; ++i ) {
    newitem->resl[i] = aresl[i];
    newitem->rell[i] = arell[i];
    strcpy(newitem->atoml[i],anamel[i]);
  }
  newitem->del = del;
  newitem->next = defs->buildres->cmaps;
  defs->buildres->cmaps = newitem;
  if ( ! defs->cmaps_present ) {
    topo_defs_log_error(defs,"cross-term entries present in topology definitions");
  }
  defs->cmaps_present = 1;
  return 0;
}

int topo_defs_exclusion(topo_defs *defs, const char *rname, int del,
	const char *a1name, int a1res, int a1rel,
	const char *a2name, int a2res, int a2rel) {
  topo_defs_exclusion_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for explicit exclusion");
    return -1;
  }
  if ( NAMETOOLONG(a1name) ) return -2;
  if ( NAMETOOLONG(a2name) ) return -3;
  if ( del && ! defs->buildres->patch ) return -4;
  if ( ( a1res || a2res ) && ! defs->buildres->patch ) return -4;
  newitem = (topo_defs_exclusion_t*) malloc(sizeof(topo_defs_exclusion_t));
  if ( ! newitem )  return -5;
  newitem->res1 = a1res;
  newitem->rel1 = a1rel;
  newitem->res2 = a2res;
  newitem->rel2 = a2rel;
  newitem->del = del;
  
#if !defined(NEWPSFGEN)  

  strcpy(newitem->atom1,a1name);
  strcpy(newitem->atom2,a2name);

#else

  newitem->atom1 = topo_defs_find_atom_by_name(defs, a1name, a1res, a1rel);
  newitem->atom2 = topo_defs_find_atom_by_name(defs, a2name, a2res, a2rel);
  strcpy(newitem->atomstr1,a1name);
  strcpy(newitem->atomstr2,a2name);

#endif

  newitem->next = defs->buildres->exclusions;
  defs->buildres->exclusions = newitem;
  return 0;
}

int topo_defs_conformation(topo_defs *defs, const char *rname, int del,
	const char *a1name, int a1res, int a1rel,
	const char *a2name, int a2res, int a2rel,
	const char *a3name, int a3res, int a3rel,
	const char *a4name, int a4res, int a4rel,
	double dist12, double angle123, double dihedral, int improper,
	double angle234, double dist34) {
  topo_defs_conformation_t *newitem;
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for conformation");
    return -1;
  }
  if ( NAMETOOLONG(a1name) ) return -2;
  if ( NAMETOOLONG(a2name) ) return -3;
  if ( NAMETOOLONG(a3name) ) return -4;
  if ( NAMETOOLONG(a4name) ) return -5;
  if ( del && ! defs->buildres->patch ) return -6;
  if ( ( a1res || a2res || a3res || a4res ) &&
			! defs->buildres->patch ) return -7;
  newitem = (topo_defs_conformation_t*)malloc(sizeof(topo_defs_conformation_t));
  if ( ! newitem )  return -8;
  newitem->res1 = a1res;
  newitem->rel1 = a1rel;
  newitem->res2 = a2res;
  newitem->rel2 = a2rel;
  newitem->res3 = a3res;
  newitem->rel3 = a3rel;
  newitem->res4 = a4res;
  newitem->rel4 = a4rel;
  newitem->del = del;
  newitem->improper = improper;
  newitem->dist12 = dist12;
  newitem->angle123 = angle123;
  newitem->dihedral = dihedral;
  newitem->angle234 = angle234;
  newitem->dist34 = dist34;
  
#if !defined(NEWPSFGEN)  

  strcpy(newitem->atom1,a1name);
  strcpy(newitem->atom2,a2name);
  strcpy(newitem->atom3,a3name);
  strcpy(newitem->atom4,a4name);

#else

  newitem->atom1 = topo_defs_find_atom_by_name(defs, a1name, a1res, a1rel);
  newitem->atom2 = topo_defs_find_atom_by_name(defs, a2name, a2res, a2rel);
  newitem->atom3 = topo_defs_find_atom_by_name(defs, a3name, a3res, a3rel);
  newitem->atom4 = topo_defs_find_atom_by_name(defs, a4name, a4res, a4rel);
  strcpy(newitem->atomstr1,a1name);
  strcpy(newitem->atomstr2,a2name);
  strcpy(newitem->atomstr3,a3name);
  strcpy(newitem->atomstr4,a4name);

#endif

  newitem->next = defs->buildres->conformations;
  defs->buildres->conformations = newitem;
  return 0;
}


int topo_defs_default_patching_first(topo_defs *defs, const char *pname) {
  if ( ! defs ) return -1;
  if ( NAMETOOLONG(pname) ) return -2;
  strcpy(defs->pfirst,pname);
  return 0;
}

int topo_defs_default_patching_last(topo_defs *defs, const char *pname) {
  if ( ! defs ) return -1;
  if ( NAMETOOLONG(pname) ) return -2;
  strcpy(defs->plast,pname);
  return 0;
}

int topo_defs_patching_first(topo_defs *defs, const char *rname,
	const char *pname) {
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for patching");
    return -1;
  }
  if ( NAMETOOLONG(pname) ) return -2;
  strcpy(defs->buildres->pfirst,pname);
  return 0;
}


int topo_defs_patching_last(topo_defs *defs, const char *rname,
	const char *pname) {
  if ( ! defs ) return -1;
  if ( ! defs->buildres ) {
    if ( defs->buildres_no_errors ) return 0;
    topo_defs_log_error(defs,"no residue in progress for patching");
    return -1;
  }
  if ( NAMETOOLONG(pname) ) return -2;
  strcpy(defs->buildres->plast,pname);
  return 0;
}

/* int topo_defs_add_topofile(topo_defs *defs, const char *filename) { */
/*   topo_defs_topofile_t **topofiles; */
/*   topo_defs_topofile_t *topofiletmp; */
/*   if ( ! defs ) return -1; */
/*   if ( strlen(filename)>=256 ) return -2; */
/*   topofiles = &(defs->topofiles); */
/*   topofiletmp = 0; */
/*   topofiletmp = memarena_alloc(defs->arena,sizeof(topo_defs_topofile_t)); */
/*   if ( ! topofiletmp ) return -3; */
  
/*   strcpy(topofiletmp->filename,filename); */

/*   printf("add_topo %i %s;\n", defs->ntopo, topofiletmp->filename);  */
/*   defs->ntopo++; */
/*   topofiletmp->next = *topofiles; */
/*   *topofiles = topofiletmp; */
/*   return 0; */
/* } */

int topo_defs_add_topofile(topo_defs *defs, const char *filename) {
/*   topo_defs_topofile_t **topofiles; */
/*   topo_defs_topofile_t *topofiletmp; */
/*   if ( ! defs ) return -1; */
/*   if ( strlen(filename)>=256 ) return -2; */
/*   topofiles = &(defs->topofiles); */
/*   topofiletmp = 0; */
/*   topofiletmp = memarena_alloc(defs->arena,sizeof(topo_defs_topofile_t)); */
/*   if ( ! topofiletmp ) return -3; */
  
/*   strcpy(topofiletmp->filename,filename); */

/*   printf("add_topo %i %s;\n", defs->ntopo, topofiletmp->filename);  */
/*   defs->ntopo++; */
/*   topofiletmp->next = *topofiles; */
/*   *topofiles = topofiletmp; */
/*   return 0; */

  int i;
  topo_defs_topofile_t *newitem;
  char errmsg[64 + 256];
  if ( ! defs ) return -1;
  if ( strlen(filename)>=256 ) return -2;
  if ( ( i = hasharray_index(defs->topo_hash,filename) ) != HASHARRAY_FAIL ) {
    sprintf(errmsg,"duplicate topology file %s",filename);
    topo_defs_log_error(defs,errmsg);
    newitem = &defs->topo_array[i];
  } else {
    i = hasharray_insert(defs->topo_hash,filename);
    if ( i == HASHARRAY_FAIL ) return -4;
    newitem = &defs->topo_array[i];
    strcpy(newitem->filename,filename);
  }
  return 0;
}

