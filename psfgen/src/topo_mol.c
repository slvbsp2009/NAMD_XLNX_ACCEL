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
 *      $RCSfile: topo_mol.c,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.69 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> /* for isspace() */
#include "topo_defs_struct.h"
#include "topo_mol_struct.h"

#if defined(NEWPSFGEN)
#include "psfgen.h"
#endif

#if defined(_MSC_VER)
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif

#if defined(NEWPSFGEN)

static topo_mol_atom_t * topo_mol_find_atom_by_name(topo_mol_residue_t *res, 
  const char *aname);

static int topo_mol_find_atomindex(topo_mol_residue_t *res, 
    const char *aname);
    
#endif    
topo_mol * topo_mol_create(topo_defs *defs) {
  topo_mol *mol;
  if ( ! defs ) return 0;
  if ( (mol = (topo_mol*) malloc(sizeof(topo_mol))) ) {
    mol->newerror_handler_vdata = 0;
    mol->newerror_handler_inter = 0;
    mol->newerror_handler = 0;
    mol->defs = defs;
    mol->npatch = 0;
    mol->patches = 0;
    mol->curpatch = 0;
    mol->segment_hash = hasharray_create(
	(void**) &(mol->segment_array), sizeof(topo_mol_segment_t*));
    mol->buildseg = 0;
    mol->arena = memarena_create();
    mol->angle_arena = memarena_create();
    mol->dihedral_arena = memarena_create();
    if ( ! mol->segment_hash || ! mol->arena ) {
      topo_mol_destroy(mol);
      return 0;
    }
#if defined(NEWPSFGEN)
    /* initialize the drude and lonepairs flags to 0 */
    mol->drude = 0;
    mol->lonepairs = 0;
#endif
  }
  return mol;
}

void topo_mol_destroy(topo_mol *mol) {
  int i,n;

#if defined(NEWPSFGEN)
  int j,nres;
  topo_mol_residue_t *res;
#endif
  topo_mol_segment_t *s;
  
  if ( ! mol ) return;

  n = hasharray_count(mol->segment_hash);
  for ( i=0; i<n; ++i ) {
    s = mol->segment_array[i];
    if ( ! s ) continue;
#if defined(NEWPSFGEN)
    nres = hasharray_count(s->residue_hash);
    for (j = 0; j < nres; ++j) {
      res = &(s->residue_array[j]);
      free((void*)res->atomArray);
    }
#endif
    hasharray_destroy(s->residue_hash);
  }
  hasharray_destroy(mol->segment_hash);
  memarena_destroy(mol->arena);
  memarena_destroy(mol->angle_arena);
  memarena_destroy(mol->dihedral_arena);
  free((void*)mol);
}

void topo_mol_error_handler(topo_mol *mol, void *vdata, void *v, void (*print_msg)(void *, void *,const char *)) {
  if ( mol ) {
    mol->newerror_handler = print_msg;
    mol->newerror_handler_inter = v;
    mol->newerror_handler_vdata = vdata;
  }
}

/* internal method */
static void topo_mol_log_error(topo_mol *mol, const char *msg) {
  if (mol && msg && mol->newerror_handler)
    mol->newerror_handler(mol->newerror_handler_vdata, mol->newerror_handler_inter, msg);
}

static topo_mol_segment_t * topo_mol_get_seg(topo_mol *mol,
			const topo_mol_ident_t *target) {
  int iseg;
  char errmsg[64 + 3*NAMEMAXLEN];

  if ( ! mol ) return 0;
  iseg = hasharray_index(mol->segment_hash,target->segid);
  if ( iseg == HASHARRAY_FAIL ) {
    sprintf(errmsg,"no segment %s",target->segid);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }
  return mol->segment_array[iseg];
}

static topo_mol_residue_t * topo_mol_get_res(topo_mol *mol,
			const topo_mol_ident_t *target, int irel) {
  int nres, ires;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  char errmsg[64 + 3*NAMEMAXLEN];
  seg = topo_mol_get_seg(mol,target);
  if ( ! seg ) return 0;
  nres = hasharray_count(seg->residue_hash);
  ires = hasharray_index(seg->residue_hash,target->resid);
  if ( ires == HASHARRAY_FAIL ) {
    sprintf(errmsg,"no residue %s of segment %s",
					target->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }
  if ( (ires+irel) < 0 || (ires+irel) >= nres ) {
    res = seg->residue_array + ires;
    if ( irel < 0 )
      sprintf(errmsg,"no residue %d before %s:%s of segment %s",
		-1*irel,res->name,res->resid,target->segid);
    if ( irel > 0 )
      sprintf(errmsg,"no residue %d past %s:%s of segment %s",
		irel,res->name,res->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }

  return (seg->residue_array + ires + irel);
}


static topo_mol_atom_t * topo_mol_get_atom(topo_mol *mol,
			const topo_mol_ident_t *target, int irel) {
#if defined(NEWPSFGEN)
  int i=0;
#endif
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom = 0;
  char errmsg[64 + 3*NAMEMAXLEN];
  res = topo_mol_get_res(mol,target,irel);
  if ( ! res ) return 0;
  
#if !defined(NEWPSFGEN)
  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }
#else
  for ( i = 0; i < res->atomSize; i++ ) {
    if ( ! strcmp(target->aname,res->atomArray[i]->name) ) {
      atom = res->atomArray[i];
      break;
    }
  }
#endif
  
  if ( ! atom ) {
    sprintf(errmsg,"no atom %s in residue %s:%s of segment %s",
		target->aname,res->name,res->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
  }
  return atom;
}


static topo_mol_atom_t *topo_mol_get_atom_from_res(
    const topo_mol_residue_t *res, const char *aname) {
    topo_mol_atom_t *atom = 0;
    
#if !defined(NEWPSFGEN)

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(aname,atom->name) ) break;
  }
  
#else

  int i;
  for ( i = 0; i < res->atomSize; i++ ) {
    atom = res->atomArray[i];
    if ( ! strcmp(aname,atom->name) ) return atom;
  }  
    
#endif

  return atom;
}

int topo_mol_segment(topo_mol *mol, const char *segid) {
  int i;
  topo_mol_segment_t *newitem;
  char errmsg[32 + NAMEMAXLEN];
  if ( ! mol ) return -1;
  mol->buildseg = 0;
  if ( NAMETOOLONG(segid) ) return -2;
  if ( ( i = hasharray_index(mol->segment_hash,segid) ) != HASHARRAY_FAIL ) {
    sprintf(errmsg,"duplicate segment key %s",segid);
    topo_mol_log_error(mol,errmsg);
    return -3;
  } else {
    i = hasharray_insert(mol->segment_hash,segid);
    if ( i == HASHARRAY_FAIL ) return -4;
    newitem = mol->segment_array[i] = (topo_mol_segment_t*)
		memarena_alloc(mol->arena,sizeof(topo_mol_segment_t));
    if ( ! newitem ) return -5;
  }
  strcpy(newitem->segid,segid);
  newitem->residue_hash = hasharray_create(
	(void**) &(newitem->residue_array), sizeof(topo_mol_residue_t));
  strcpy(newitem->pfirst,"");
  strcpy(newitem->plast,"");
  newitem->auto_angles = mol->defs->auto_angles;
  newitem->auto_dihedrals = mol->defs->auto_dihedrals;
  mol->buildseg = newitem;
  return 0;
}

int topo_mol_segment_first(topo_mol *mol, const char *rname) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for first patch");
    return -1;
  }
  if ( NAMETOOLONG(rname) ) return -2;
  strcpy(mol->buildseg->pfirst,rname);
  return 0;
}

int topo_mol_segment_last(topo_mol *mol, const char *rname) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for last patch");
    return -1;
  }
  if ( NAMETOOLONG(rname) ) return -2;
  strcpy(mol->buildseg->plast,rname);
  return 0;
}

int topo_mol_segment_auto_angles(topo_mol *mol, int autogen) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for auto angles");
    return -1;
  }
  mol->buildseg->auto_angles = autogen;
  return 0;
}

int topo_mol_segment_auto_dihedrals(topo_mol *mol, int autogen) {
  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for auto dihedrals");
    return -1;
  }
  mol->buildseg->auto_dihedrals = autogen;
  return 0;
}

int topo_mol_residue(topo_mol *mol, const char *resid, const char *rname,
						const char *chain) {
  int i;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *newitem;
  char errmsg[32 + NAMEMAXLEN];

  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for residue");
    return -1;
  }
  seg = mol->buildseg;
  if ( NAMETOOLONG(resid) ) return -2;
  if ( NAMETOOLONG(rname) ) return -3;
  if ( hasharray_index(seg->residue_hash,resid) != HASHARRAY_FAIL ) {
    sprintf(errmsg,"duplicate residue key %s",resid);
    topo_mol_log_error(mol,errmsg);
    return -3;
  }

  if ( hasharray_index(mol->defs->residue_hash,rname) == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown residue type %s",rname);
    topo_mol_log_error(mol,errmsg);
  }

  i = hasharray_insert(seg->residue_hash,resid);
  if ( i == HASHARRAY_FAIL ) return -4;
  newitem = &(seg->residue_array[i]);
  strcpy(newitem->resid,resid);
  strcpy(newitem->name,rname);
  strcpy(newitem->chain,chain);
  
#if !defined(NEWPSFGEN)
  newitem->atoms = 0;
#endif

  return 0;
}

int topo_mol_mutate(topo_mol *mol, const char *resid, const char *rname) {
  int ires;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  char errmsg[32 + 3*NAMEMAXLEN];

  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for mutate");
    return -1;
  }
  seg = mol->buildseg;

  if ( NAMETOOLONG(resid) ) return -2;
  if ( NAMETOOLONG(rname) ) return -3;
  ires = hasharray_index(seg->residue_hash,resid);
  if ( ires == HASHARRAY_FAIL ) {
    sprintf(errmsg,"residue %s does not exist",resid);
    topo_mol_log_error(mol,errmsg);
    return -1;
  }
  res = seg->residue_array + ires;
  sprintf(errmsg,"mutating residue %s from %s to %s",resid,res->name,rname);
  topo_mol_log_error(mol,errmsg);

  if ( hasharray_index(mol->defs->residue_hash,rname) == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown residue type %s",rname);
    topo_mol_log_error(mol,errmsg);
  }

  strcpy(res->name,rname);

  return 0;
}

#if !defined(NEWPSFGEN) 

static topo_mol_atom_t * topo_mol_unlink_atom(
		topo_mol_atom_t **atoms, const char *aname) {
  topo_mol_atom_t **atom;
  topo_mol_atom_t *oldatom;
  if ( ! atoms ) return 0;
  for ( atom = atoms ; *atom; atom = &((*atom)->next) ) {
    if ( ! strcmp(aname,(*atom)->name) ) break;
  }
  oldatom = *atom;
  if ( *atom ) *atom = ((*atom)->next);
  return oldatom;
}

#else

static topo_mol_atom_t * topo_mol_unlink_atom(
		topo_mol_residue_t *res, const char *aname) {

  topo_mol_atom_t *atom = topo_mol_find_atom_by_name(res, aname);
  if (atom) {
    /*
     * Delete the atom by point the pointer to the next element,
     * and do the same until a null pointer is reached
     * pointer1 = Atom1 pointer2=Atom2 pointer3=Atom3 pointer4 = Null
     * to delete atom2:
     * pointer1 = Atom1 pointer2=Atom3 pointer3=Null
    */
    topo_mol_atom_t **atomaux=NULL;
    topo_mol_atom_t **atoprev=NULL;
    int found = 0;
    for (atomaux = res->atomArray; *atomaux; atomaux++) {
      if (*atomaux == atom) {
        found = 1;
        atoprev = atomaux;
        ++atoprev;
        --res->atomSize;
      }
      if (found) {
       *atomaux = *atoprev;
       if (!*atoprev) {
         break;
       } else {
         /*
          * if the deleted atom is not at the end of the array, the order
          * of the array changed
         */
         res->reordered = 1;
       }
       ++atoprev;
      }
    }
    atom->del = 1;
  } 

  return atom;

}

#endif

#if !defined(NEWPSFGEN)

static topo_mol_atom_t * topo_mol_find_atom(topo_mol_atom_t **newatoms,
		topo_mol_atom_t *oldatoms, const char *aname) {
  topo_mol_atom_t *atom, **newatom;
  if ( ! oldatoms ) return 0;
  for ( atom = oldatoms; atom; atom = atom->next ) {
    if ( ! strcmp(aname,atom->name) ) break;
  }
  if ( atom && *newatoms != oldatoms ) {
    for ( newatom = newatoms; *newatom != oldatoms; newatom = &(*newatom)->next );
    *newatom = atom->next;
    atom->next = *newatoms;
    *newatoms = oldatoms;
  }
  return atom;
}

#else

static topo_mol_atom_t * topo_mol_find_atom_by_name(topo_mol_residue_t *res, 
                                                    const char *aname) {
  int i=0;
  topo_mol_atom_t *atom = 0;
  if ( ! res->atomArray ) return 0;
  for ( i = 0; i < res->atomSize; i++ ) {
    if (!res->atomArray[i]) continue;
    if ( ! strcmp(aname, res->atomArray[i]->name) ) {
      atom = res->atomArray[i];
      break;
    }
  }

  return atom;
}


/* Same as topo_mol_find_atom_by_name but returning the 
 * the atom index. The reason for two different procedures
 * is to avoid additional operations to retreave the atom
 * from the array position when the purpose is to get the atom
 */
static int topo_mol_find_atomindex(topo_mol_residue_t *res, 
                                   const char *aname) {
  int i=0;
  int atomindex = -1;
  if ( ! res->atomArray ) return 0;
  for ( i = 0; i < res->atomSize; i++ ) {
    if (!res->atomArray[i]) continue;
    if ( ! strcmp(aname, res->atomArray[i]->name) ) {
      atomindex = i;
      break;
    }
  }
  return atomindex;
}
#endif

#if !defined(NEWPSFGEN) 
static int topo_mol_add_atom(topo_mol *mol, topo_mol_atom_t **atoms,
                             topo_mol_atom_t *oldatoms, 
                             topo_defs_atom_t *atomdef) {

#else

static int topo_mol_update_lonepair(topo_mol *mol, topo_defs_atom_t *atomdef, 
                                    topo_mol_residue_t *res, 
                                    const topo_mol_ident_t *target , int patch);

// Compare two integers to be used in qsort
int cmpfunc (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}


static int topo_mol_add_atom(topo_mol *mol, topo_defs_residue_t *defres, 
                             topo_defs_atom_t *atomdef, 
                             topo_mol_residue_t *res) {
  topo_defs_residue_t *prevres;
  topo_mol_lonepair_t *lonepair;
  int signalpha;
  int atomindex;
#endif

  int idef;
  topo_mol_atom_t *atomtmp;
  topo_defs_type_t *atype;
  char errmsg[128];
  
#if !defined(NEWPSFGEN) 
  if ( ! mol || ! atoms ) return -1;
#else
  if ( ! mol ) return -1;
#endif

  idef = hasharray_index(mol->defs->type_hash,atomdef->type);
  if ( idef == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown atom type %s",atomdef->type);
    topo_mol_log_error(mol,errmsg);
    return -3;
  }
  atomtmp = 0;
  
#if !defined(NEWPSFGEN)
  if ( oldatoms ) atomtmp = topo_mol_find_atom(atoms,oldatoms,atomdef->name);
#else 
  atomindex = atomdef->atomIndex;
  /* if defining a patch */
  if ( defres ) {
    
    int newlimit = 0;
    int i, j;

    /* Check if the record for the patch is empty or if the name of the
     * current patch is different from the one in the record. This ensures
     * that the atomArray is resized only once per patch. 
     */
    if (res->pres[0] == '\0' || strcmp(res->pres,defres->name)) {  
      int newatoms = 0;
      int newatompos = 0; // flag to identify new atoms
      int *newpos = NULL;
      int indaux;
      int added, atmind;
      topo_defs_atom_t *defatom;
      int numdefatom = 0; /**number of atoms defined in the patch */
      int prevatomnum = 0;
      
      int residef = hasharray_index(mol->defs->residue_hash,res->name);
      if ( residef == HASHARRAY_FAIL ) {
        sprintf(errmsg,".something went really wrong. Unknown residue %s",res->name);
        topo_mol_log_error(mol,errmsg);
        return -3;
      }
      prevres = &(mol->defs->residue_array[residef]);
      prevatomnum = prevres->atomNum; 
      /* evaluate how many patch atoms are already defined in the residue.
       * when the atom is redefined in the patch, no new atom is added, the 
       * current definition is updated on the atom atomtmp
       */
    
      /* store the position of the atoms */
      newpos = (int*) malloc(defres->atomNum * sizeof(int));
      for (i = 0; i < defres->atomNum; i++) {
        newpos[i] = -2;
      }
      
      // get the exact number of new atoms
      for (defatom = defres->atoms; defatom ; defatom = defatom->next) {
        
        /* Skip the declaration of invalid atoms like the test_hack.tcl.
         * For some reason the old version of psfgen was accepting atoms declared
         * as FOO and without atom type in the definition of Patches. I would 
         * recommend just bail out with a message for all the atoms without type.
         */
        if (topo_mol_find_atomindex(res, defatom->name) == -1 && \
          defatom->type[0] == '\0') continue;

        // Skip the atom deletion declaration
        if (!defatom->del ) {
          indaux = topo_mol_find_atomindex(res, defatom->name);
          if ( indaux > -1) {
            if (newatompos) {
              
              for (j = 0; j < newatompos; j++) {
                /* the position of the new atom will be the 
                 * indaux (index of the existing atom) + 1 (actual position of the new atom) + 
                 * the position in the buffer of consecutive new atoms "newatoms"
                */
                newpos[newatoms] = indaux;
                ++newatoms;
              }
              newatompos = 0;
            }
            --numdefatom;
          } else if (atomdef->res == defatom->res) {
            ++newatompos;
          }
          ++numdefatom;
        } else {
          --prevatomnum;
        }
      }
      /* only continue if new atoms were defined... */
      if ( newatoms > 0 || newatompos > 0) {
        atomindex += prevatomnum;
        /* if the patch has no common atom with the patch, the reference
         * position to place the new atoms ins not found. In this case, place 
         * at the beginning of the residue
         */
        if (newpos[0] == -2) {
          for (j = 0; j < numdefatom; j++) {
            newpos[j] = -1;
          }
          newatoms = numdefatom;
        }
        if (newatoms > 0) {
          // the +1 is to account the null pointer at the end of the array
          newlimit = res->atomSize + newatoms +1;
          res->atomArray = (topo_mol_atom_t **) 
                realloc(res->atomArray, (newlimit)*sizeof(topo_mol_atom_t*));
          
          if (!res->atomArray) {
            sprintf(errmsg,"Failed to reallocate residue %s array",
              res->name);
            topo_mol_log_error(mol,errmsg);
            return -3;
          }
          for  (i = res->atomSize; i < newlimit; i++) {
            res->atomArray[i] = 0;
          }
          
          if (newatoms > 0) {
            topo_mol_atom_t ** tmp = NULL;

            qsort(newpos, newatoms, sizeof(int), cmpfunc);

            /* tmp array that will contain the end result of the insertions of
             * the empty spaces for the new atoms and the atoms not modified
             */
            tmp = (topo_mol_atom_t **) calloc(1, newlimit * sizeof(topo_mol_atom_t*));
            added = 0; /* store how many atoms were already added, so we can
                        * account with the difference in the original position
                        * and the position in the new residue
                        */
            atmind = 0; /* indexes of the unchanged atom array */          
            for (i = 0; i < res->atomSize + newatoms; i++) {
              /* Open space in the array to add the new atoms. This way multiple insertions
               * are allowed in an arbitrary order, which might be the case of patches 
               * defined by the users.
               * the condition added < newatoms ensures that the added never
               * goes beyond the upper limit of the newpos array
               */
              if (added < newatoms && newpos[added] == -1 ) {
                tmp[i] = NULL;
                added++;
                continue;
              }
              memcpy(tmp + i, res->atomArray + atmind, sizeof(topo_mol_atom_t*));
              if (added < newatoms && atmind == newpos[added]) {
                if (newpos[added] > -1) {
                  while (newpos[added] == atmind) {
                    i++;
                    tmp[i] = NULL;
                    added++; 
                  }
                } 
              }
              ++atmind;
            }
            memcpy(res->atomArray, tmp, newlimit * sizeof(topo_mol_atom_t*));
            free(tmp);
          }
          res->atomSize = newlimit -1;
          res->reordered = 1;
        }
      } else {
        newlimit = res->atomSize;
      }

      free((void*)newpos);
    } else {
      newlimit = res->atomSize;
    }
    // -2 to skip the null pointer and because the indexes go until newlimit -1 
    for  (i = newlimit -2; i > 0; i--) {
      if (!res->atomArray[i]) {
        break;
      }
    }
    atomindex = i;
    atomtmp = topo_mol_find_atom_by_name(res, atomdef->name);
    if (atomtmp && atomtmp->isdrudlonepair == ISLONEPAIR) {
      lonepair = atomtmp->lonepair;
    }
  }  
  
#endif

  if ( ! atomtmp ) {
    atomtmp = memarena_alloc(mol->arena,sizeof(topo_mol_atom_t));
    if ( ! atomtmp ) return -2;
    strcpy(atomtmp->name,atomdef->name);
    atomtmp->bonds = 0;
    atomtmp->angles = 0;
    atomtmp->dihedrals = 0;
    atomtmp->impropers = 0;
    atomtmp->cmaps = 0;
    atomtmp->exclusions = 0;
    atomtmp->conformations = 0;
    atomtmp->x = 0;
    atomtmp->y = 0;
    atomtmp->z = 0;
    atomtmp->vx = 0;
    atomtmp->vy = 0;
    atomtmp->vz = 0;
    atomtmp->xyz_state = TOPO_MOL_XYZ_VOID;
    atomtmp->partition = 0;
    atomtmp->atomid = 0;

#if !defined(NEWPSFGEN)
    atomtmp->next = *atoms;
    *atoms = atomtmp;
#else
    atomtmp->del = 0;
    res->atomArray[atomindex] = atomtmp;
    
    /* only increment the variable when applying patchres
     * as the variable is used to compute array sizes when realloc for the patch
     */
    if (!defres) {
      res->atomSize++;
    }
#endif

  } 

  atomtmp->copy = 0;
  atomtmp->charge = atomdef->charge;
  strcpy(atomtmp->type,atomdef->type);
  atype = &(mol->defs->type_array[idef]);
  strcpy(atomtmp->element,atype->element);
  atomtmp->mass = atype->mass;
  
#if defined(NEWPSFGEN)
  /* Now deal with the lonepair info*/
  if (atomdef->islonepair == ISLONEPAIR) {
    lonepair = memarena_alloc(mol->arena,sizeof(topo_mol_lonepair_t));
    lonepair->atoms = 0;
    if (atomdef->lonepair) {
      lonepair->distance = atomdef->lonepair->distance;
      lonepair->angle = atomdef->lonepair->angle;
      lonepair->dihedral = atomdef->lonepair->dihedral;
      lonepair->lptype = atomdef->lonepair->lptype;
      atomtmp->lonepair = lonepair;
      atomtmp->isdrudlonepair = ISLONEPAIR; 
    }
    /* flag the molecule as containing lonepairs */
    mol->lonepairs = 1;
  } else {
    atomtmp->lonepair = 0;
    atomtmp->isdrudlonepair = 0;
  }

  /* Update the drude info */
  strcpy(atomtmp->dname,atomdef->dname);
  atomtmp->alpha = atomdef->alpha;
  atomtmp->thole = atomdef->thole;
  atomtmp->dxyz = 0; /* coordinates of the lone pair */
  if (atomtmp->alpha) {
    /* Drude particle charge
    * KDRUDE is the force constant (in kcal/mol/Angst**2) for the bond
    * between. Default 500 kcal/mol/Angst**2
    * CCELEC =332.071600 (Coulomb's constant CHARMM const value);
    * 
    * q = sqrt( 2*KDRUDE * alpha / CCELEC ) * sign(alpha)
    */
    signalpha = (atomtmp->alpha < 0) ? -1 : (atomtmp->alpha > 0);
    atomtmp->dcharge = sqrt((2*K_DRUDE * atomtmp->alpha * signalpha) / 332.0636) * signalpha;
    
    /* subtract the drude charge and mass from the host charge and mass. These
     *  updates can be made upfront rather than during the write out process.
    */
    atomtmp->charge -= atomtmp->dcharge;
    atomtmp->mass -= 0.4;

    /* flag the molecule as containing drude particles */
    mol->drude = 1;
  } else {
    atomtmp->dcharge = 0;
  }
#endif

  return 0;
}


#if defined(NEWPSFGEN) 
static int topo_mol_update_lonepair(topo_mol *mol, topo_defs_atom_t *atomdef, 
                                    topo_mol_residue_t *res, 
                                    const topo_mol_ident_t *target, int patch) {
  topo_mol_ident_t tmpt;
  topo_mol_atom_t *tmpatom = NULL;
  topo_mol_atom_t *atomlp = NULL;
  int i, idef;
  char errmsg[128];
  
  if ( ! mol ) return -1;
  
  idef = hasharray_index(mol->defs->type_hash,atomdef->type);
  if ( idef == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown atom type %s",atomdef->type);
    topo_mol_log_error(mol,errmsg);
    return -3;
  }

  /* if residue was not reordered or if this is not a patch, the 
   * atom index of the topology correspondes to atom index in the molecule
  */
  if (!patch && !res->reordered) {
    atomlp = res->atomArray[atomdef->atomIndex];
  } else {
    if (!target) {
      /* target == 0 means that this process is done after 
       * expanding the atomArray during the topo_mol_add_atom
       * all the current lone pairs are part of the same residue
      */
      atomlp = topo_mol_find_atom_by_name(res,atomdef->name);
    } else {
      tmpt = target[atomdef->res];
      tmpt.aname = atomdef->name;
      atomlp = topo_mol_get_atom(mol,&tmpt,atomdef->rel);
    }
    
  }
  if ( !atomlp ) return -1;
  
  if (!atomlp->lonepair->atoms) {
    // atomlp->lonepair->atoms = 
    // (topo_mol_atom_t **)malloc((atomdef->lonepair->numatoms)*sizeof(topo_mol_atom_t*));
    atomlp->lonepair->atoms = 
    (topo_mol_atom_t **)memarena_alloc(mol->arena,(atomdef->lonepair->numatoms)*sizeof(topo_mol_atom_t*));
    if (!atomlp->lonepair->atoms) {
      return -1;
    }
  }
  
  for (i = 0; i < atomdef->lonepair->numatoms; i++) {
    if (!patch && !res->reordered) {
      /* In the mol, I am storing the lonepair as atom[0] so I can just have a list of of lone pairs
       * and being able to recosntruct the lonepair psf section. This is useful as the 
       * the additive forcefiel has only a few lonepairs compared with the rest.
       * The drude forcefiel, there are many, so an independet list is not necessary
      */
      atomlp->lonepair->atoms[i] = res->atomArray[atomdef->lonepair->atoms[i]->atomIndex];
      
    } else {
      
      if (!target) {
        tmpatom = topo_mol_find_atom_by_name(res,atomdef->lonepair->atoms[i]->name);
      } else {
        tmpt = target[atomdef->lonepair->atoms[i]->res];
        tmpt.aname = atomdef->lonepair->atoms[i]->name;
        tmpatom = topo_mol_get_atom(mol,&tmpt,atomdef->lonepair->atoms[i]->rel);
      }
      if ( !tmpatom ) return -1;
      
      atomlp->lonepair->atoms[i] = tmpatom;
    }
  }

  return 0;  
}



static int topo_mol_add_anisotropy(topo_mol *mol,
                      topo_defs_anisotropy_t *aniso, topo_mol_residue_t *res, 
                      const topo_mol_ident_t *target) {
                      topo_mol_ident_t tmpt;

  topo_mol_anisotropy_t *newitem;
  topo_mol_atom_t *atom;
  int i;
   
  if ( ! mol ) return -1;
  /* if residue was not reordered or if this is not a patch, the 
  * atom index of the topology correspondes to atom index in the molecule
  */
  newitem = memarena_alloc(mol->arena,sizeof(topo_mol_anisotropy_t));
  // newitem->atoms = (topo_mol_atom_t **)malloc(4*sizeof(topo_mol_atom_t*));
  newitem->atoms = memarena_alloc(mol->arena, 4*sizeof(topo_mol_atom_t*));  
  if (!newitem->atoms) {
    return -1;
  }
  for (i = 0; i < 4; i++) {
     atom = 0;
    if (!aniso->patch && !res->reordered && aniso->atoms[i]) {
      atom = res->atomArray[aniso->atoms[i]->atomIndex];
    } else {

      tmpt = target[aniso->res[i]];
      tmpt.aname = aniso->atomsname[i];
      atom = topo_mol_get_atom(mol,&tmpt,aniso->rel[i]);
    
    }
    if (!atom) {
      return -2;
    }
    newitem->atoms[i] = atom;
  }

  newitem->k11 = aniso->k11;
  newitem->k22 = aniso->k22;
  newitem->k33 = aniso->k33;
  newitem->del = 0;
  newitem->next = res->aniso;
  res->aniso = newitem;

  ++res->numaniso ;
  return 0;
}


#endif

topo_mol_bond_t * topo_mol_bond_next(
		topo_mol_bond_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  return 0;
}

topo_mol_angle_t * topo_mol_angle_next(
		topo_mol_angle_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  return 0;
}

#if !defined(NEWPSFGEN)
topo_mol_dihedral_t * topo_mol_dihedral_next(
		topo_mol_dihedral_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];

  return 0;
}
#endif

#if !defined(NEWPSFGEN)
topo_mol_improper_t * topo_mol_improper_next(
		topo_mol_improper_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];
  return 0;
}
#endif
topo_mol_cmap_t * topo_mol_cmap_next(
		topo_mol_cmap_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];
  if ( tuple->atom[4] == atom ) return tuple->next[4];
  if ( tuple->atom[5] == atom ) return tuple->next[5];
  if ( tuple->atom[6] == atom ) return tuple->next[6];
  if ( tuple->atom[7] == atom ) return tuple->next[7];
  return 0;
}

topo_mol_exclusion_t * topo_mol_exclusion_next(
		topo_mol_exclusion_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  return 0;
}


static topo_mol_conformation_t * topo_mol_conformation_next(
		topo_mol_conformation_t *tuple, topo_mol_atom_t *atom) {
  if ( tuple->atom[0] == atom ) return tuple->next[0];
  if ( tuple->atom[1] == atom ) return tuple->next[1];
  if ( tuple->atom[2] == atom ) return tuple->next[2];
  if ( tuple->atom[3] == atom ) return tuple->next[3];
  return 0;
}

static void topo_mol_destroy_atom(topo_mol_atom_t *atom) {
  topo_mol_bond_t *bondtmp;
  topo_mol_angle_t *angletmp;
  topo_mol_dihedral_t *dihetmp;
  topo_mol_improper_t *imprtmp;
  topo_mol_cmap_t *cmaptmp;
  topo_mol_exclusion_t *excltmp;
  topo_mol_conformation_t *conftmp;
  if ( ! atom ) return;

  for ( bondtmp = atom->bonds; bondtmp;
		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
    bondtmp->del = 1;
  }

  for ( angletmp = atom->angles; angletmp;
		angletmp = topo_mol_angle_next(angletmp,atom) ) {
    angletmp->del = 1;
  }
  
#if !defined(NEWPSFGEN)
  for ( dihetmp = atom->dihedrals; dihetmp;
		dihetmp = topo_mol_dihedral_next(dihetmp,atom) ) {
    dihetmp->del = 1;
  }
#else  
  for ( dihetmp = atom->dihedrals; dihetmp;
        dihetmp = dihetmp->next ) {
    dihetmp->del = 1;
  }
#endif  

#if !defined(NEWPSFGEN)
  for ( imprtmp = atom->impropers; imprtmp;
		imprtmp = topo_mol_improper_next(imprtmp,atom) ) {
    imprtmp->del = 1;
  }
#else
  for ( imprtmp = atom->impropers; imprtmp;
        imprtmp = imprtmp->next ) {
    imprtmp->del = 1;
  }
#endif  

  for ( cmaptmp = atom->cmaps; cmaptmp;
		cmaptmp = topo_mol_cmap_next(cmaptmp,atom) ) {
    cmaptmp->del = 1;
  }
  for ( excltmp = atom->exclusions; excltmp;
		excltmp = topo_mol_exclusion_next(excltmp,atom) ) {
    excltmp->del = 1;
  }
  for ( conftmp = atom->conformations; conftmp;
		conftmp = topo_mol_conformation_next(conftmp,atom) ) {
    conftmp->del = 1;
  }
  
}

#if !defined(NEWPSFGEN)

static void topo_mol_del_atom(topo_mol_residue_t *res, const char *aname) {
  if ( ! res ) return;
  topo_mol_destroy_atom(topo_mol_unlink_atom(&(res->atoms),aname));
}

#else

static void topo_mol_del_aniso(topo_mol_residue_t *res, const char *aname);


static void topo_mol_del_atom(topo_mol_residue_t *res, const char *aname) {
  if ( ! res ) return;

  topo_mol_destroy_atom(topo_mol_unlink_atom(res,aname));
  /* Delete the anisotropy if exists
   * this step is done here are we need the res reference to get the aniso data
   * Although in most cases there is no explicit delete anisotropy command, 
   * some patches delete atoms that define the anisotropy. In this cases,
   * we need to check of the deleted atom is part of an anisotropy definition
   * the only definition of "DELETE ANISOTROPY" that I found was in the 
   * toppar_drude_nucleic_acid_2017c.str, PRES DEOX, but I believe that this line
   * is just a left over from some testing. The anisotropy is redifined a few lines
   * below exactly in the same way
  */
  if (res->numaniso) {
    topo_mol_del_aniso(res,aname);
  }
}
#endif


#if defined(NEWPSFGEN)
/** Delete the any anisotropy definition in the res *res containing the atom *atomdef */
static void topo_mol_del_aniso(topo_mol_residue_t *res, const char *aname) {
  
  topo_mol_anisotropy_t *aniso;
  int i;
  for (aniso = res->aniso; aniso; aniso = aniso->next) {
    
    for (i = 0; i < 4  && !aniso->del ; i++) {
      if (!strcmp(aniso->atoms[i]->name,aname)) {
        aniso->del = 1;
        --res->numaniso;
        break;
      }
    }
  }
}

#endif  

/*
 * The add_xxx_to_residues routines exist because topo_mol_end can do
 * more intelligent error checking than what's done in the add_xxx
 * routines.  The add_xxx routines are called by topo_mol_patch, which
 * has to be more general (and more paranoid) about its input.  Returning
 * nonzero from add_xxx_to_residues is always a serious error.
 */

#if !defined(NEWPSFGEN)

static int add_bond_to_residues(topo_mol *mol, 
    const topo_mol_residue_t *res1, const char *aname1,
    const topo_mol_residue_t *res2, const char *aname2) {
      
#else 

static int add_bond_to_residues(topo_mol *mol,
    const topo_mol_ident_t *targets, 
    const topo_mol_residue_t *res1,
    const topo_mol_residue_t *res2, topo_defs_bond_t *def) {
      
#endif
 
  topo_mol_bond_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0;

#if !defined(NEWPSFGEN)

  a1 = topo_mol_get_atom_from_res(res1, aname1);
  a2 = topo_mol_get_atom_from_res(res2, aname2);
  if (!a1 || !a2) return -1;
#else

  if (res1->reordered || !def->atom1 || def->atom1->patch) {
    a1 = topo_mol_get_atom_from_res(res1, def->atomstr1);
  } else {
    a1 = res1->atomArray[def->atom1->atomIndex];
  }
  
  if (res2->reordered || !def->atom2 || def->atom2->patch) {
    a2 = topo_mol_get_atom_from_res(res2, def->atomstr2);
  } else {
    a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if (!a1 || !a2) return -1;
  if (a1->isdrudlonepair || a2->isdrudlonepair) return 0;
#endif
  
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->bonds;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->bonds;
  tuple->atom[1] = a2;
  tuple->del = 0;
  a1->bonds = tuple;
  a2->bonds = tuple;
  
  
  return 0;
}

static int topo_mol_add_bond(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_bond_t *def) {
  topo_mol_bond_t *tuple;
  topo_mol_atom_t *a1 = 0, *a2 = 0;
  topo_mol_ident_t t1, t2;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL;
#endif

  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;

#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return -3;
  }
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return -5;
  }
  if (a1->isdrudlonepair|| a2->isdrudlonepair) return 0;
#endif
  
  
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->bonds;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->bonds;
  tuple->atom[1] = a2;
  tuple->del = 0;
  a1->bonds = tuple;
  a2->bonds = tuple;

  return 0;

}

static void topo_mol_del_bond(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_bond_t *def) {
  topo_mol_bond_t *tuple;
  topo_mol_atom_t *a1 = 0, *a2 = 0;
  topo_mol_ident_t t1, t2;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL;
#endif

  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  
#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( ! a1 ) return;
  }

  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];

  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return;
  }
#endif
  
  for ( tuple = a1->bonds; tuple;
    tuple = topo_mol_bond_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2 ) tuple->del = 1;
    if ( tuple->atom[0] == a2 && tuple->atom[1] == a1 ) tuple->del = 1;
  }
  
}


static int topo_mol_add_angle(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_angle_t *def) {
  topo_mol_angle_t *tuple;
  topo_mol_atom_t *a1 = 0, *a2 = 0, *a3 = 0;
  topo_mol_ident_t t1, t2, t3;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL;
#endif

  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  
#else 

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return -3;
  }

  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];

  if (def->atom2) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return -5;
  }
  
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];

  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return -7;
  }
  if (a1->isdrudlonepair|| a2->isdrudlonepair || a3->isdrudlonepair ) return 0;
#endif

  tuple = memarena_alloc(mol->angle_arena,sizeof(topo_mol_angle_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->angles;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->angles;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->angles;
  tuple->atom[2] = a3;
  tuple->del = 0;
  a1->angles = tuple;
  a2->angles = tuple;
  a3->angles = tuple;
  return 0;
}

static void topo_mol_del_angle(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_angle_t *def) {
  topo_mol_angle_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0;
  topo_mol_ident_t t1, t2, t3;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL;
#endif

  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  
#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return;
  }
  
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];

  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( ! a2 ) return;
  }

  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];

  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return;
  }
  
#endif

  for ( tuple = a1->angles; tuple;
		tuple = topo_mol_angle_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 ) tuple->del = 1;
    if ( tuple->atom[0] == a3 && tuple->atom[1] == a2
	&& tuple->atom[2] == a1 ) tuple->del = 1;
  }
}


static int topo_mol_add_dihedral(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_dihedral_t *def) {
  topo_mol_dihedral_t *tuple;
  topo_mol_atom_t *a1 = 0, *a2 = 0, *a3 =0, *a4 =0;
  topo_mol_ident_t t1, t2, t3, t4;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL, *res4=NULL;
#endif
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  if ( ! a4 ) return -9;
  tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->dihedrals;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->dihedrals;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->dihedrals;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->dihedrals;
  tuple->atom[3] = a4;
  tuple->del = 0;
  a1->dihedrals = tuple;
  a2->dihedrals = tuple;
  a3->dihedrals = tuple;
  a4->dihedrals = tuple;

#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return -3;
  }
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  
  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return -5;
  }
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  
  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return -7;
  }
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  
  if (def->atom4 && !def->atom4->patch) {
    res4 = topo_mol_get_res(mol,&t4,def->rel4);
    if (!res4->reordered) a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if ( !a4 ) {
    t4.aname = def->atomstr4;
    a4 = topo_mol_get_atom(mol,&t4,def->rel4);
    if ( !a4 ) return -9;
  }
  
  if (a1->isdrudlonepair || a2->isdrudlonepair || 
      a3->isdrudlonepair || a4->isdrudlonepair) return 0;
  tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
  if ( ! tuple ) return -10;
  tuple->next = a1->dihedrals;
  tuple->atom[0] = a2;
  tuple->atom[1] = a3;
  tuple->atom[2] = a4;
  tuple->del = 0;
  a1->dihedrals = tuple;
#endif 
  
  
  return 0;
}

static void topo_mol_del_dihedral(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_dihedral_t *def) {
  topo_mol_dihedral_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  topo_mol_ident_t t1, t2, t3, t4;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL, *res4=NULL;
#endif
  
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  
  for ( tuple = a1->dihedrals; tuple;
		tuple = topo_mol_dihedral_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 && tuple->atom[3] == a4 ) tuple->del = 1;
    if ( tuple->atom[0] == a4 && tuple->atom[1] == a3
	&& tuple->atom[2] == a2 && tuple->atom[3] == a1 ) tuple->del = 1;
  }
  
#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return;
  }

  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  
  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( ! a2 ) return;
  }

  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];

  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return;
  }

  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];

  if (def->atom4 && !def->atom4->patch) {
    res4 = topo_mol_get_res(mol,&t4,def->rel4);
    if (!res4->reordered) a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if ( !a4) {
    t4.aname = def->atomstr4;
    a4 = topo_mol_get_atom(mol,&t4,def->rel4);
    if ( !a4 ) return;
  }
  
  for ( tuple = a1->dihedrals; tuple; tuple = a1->dihedrals->next ) {
    /* 
     * XXX this causes memory corruption because it is not 0-based indexing.
     * Fixed ~ Joao
     */
    if (tuple->atom[0] == a2 && tuple->atom[1] == a3 && tuple->atom[2] == a4) 
      tuple->del = 1;
  }
#endif

  
}

#if !defined(NEWPSFGEN)

static int add_improper_to_residues(topo_mol *mol, 
    const topo_mol_residue_t *res1, const char *aname1,
    const topo_mol_residue_t *res2, const char *aname2,
    const topo_mol_residue_t *res3, const char *aname3,
    const topo_mol_residue_t *res4, const char *aname4) {
      
#else

static int add_improper_to_residues(topo_mol *mol, 
    const topo_mol_ident_t *targets,
    const topo_mol_residue_t *res1, 
    const topo_mol_residue_t *res2,
    const topo_mol_residue_t *res3,
    const topo_mol_residue_t *res4, topo_defs_improper_t *def) {
      
#endif

  topo_mol_improper_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  
#if !defined(NEWPSFGEN)

  a1 = topo_mol_get_atom_from_res(res1, aname1);
  a2 = topo_mol_get_atom_from_res(res2, aname2);
  a3 = topo_mol_get_atom_from_res(res3, aname3);
  a4 = topo_mol_get_atom_from_res(res4, aname4);
  if (!a1 || !a2 || !a3 || !a4) return -1;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->impropers;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->impropers;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->impropers;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->impropers;
  tuple->atom[3] = a4;
  tuple->del = 0;
  a1->impropers = tuple;
  a2->impropers = tuple;
  a3->impropers = tuple;
  a4->impropers = tuple;
  
#else

  if (res1->reordered || !def->atom1 || def->atom1->patch) {
    a1 = topo_mol_get_atom_from_res(res1, def->atomstr1);
  } else {
    a1 = res1->atomArray[def->atom1->atomIndex];
  }
  
  if (res2->reordered || !def->atom2 || def->atom2->patch) {
    a2 = topo_mol_get_atom_from_res(res2, def->atomstr2);
  } else {
    a2 = res2->atomArray[def->atom2->atomIndex];
  }
  
  if (res3->reordered || !def->atom3 || def->atom3->patch) {
    a3 = topo_mol_get_atom_from_res(res3, def->atomstr3);
  } else {
    a3 = res3->atomArray[def->atom3->atomIndex];
  }
  
  if (res4->reordered || !def->atom4 || def->atom4->patch) {
    a4 = topo_mol_get_atom_from_res(res4, def->atomstr4);
  } else {
    a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if (!a1 || !a2 || !a3 || !a4) return -1;
  if (a1->isdrudlonepair || a2->isdrudlonepair || 
      a3->isdrudlonepair || a4->isdrudlonepair) return 0;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
  if ( ! tuple ) return -10;
  tuple->next = a1->impropers;
  tuple->atom[0] = a2;
  tuple->atom[1] = a3;
  tuple->atom[2] = a4;
  tuple->del = 0;
  a1->impropers = tuple;
#endif

  return 0;
}

static int topo_mol_add_improper(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_improper_t *def) {
  topo_mol_improper_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  topo_mol_ident_t t1, t2, t3, t4;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL, *res4=NULL;
#endif

  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)
  
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  if ( ! a4 ) return -9;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->impropers;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->impropers;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->impropers;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->impropers;
  tuple->atom[3] = a4;
  tuple->del = 0;
  a1->impropers = tuple;
  a2->impropers = tuple;
  a3->impropers = tuple;
  a4->impropers = tuple;
#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return -3;
  }

  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];

  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return -5;
  }
  
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  
  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return -7;
  }

  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  
  if (def->atom4 && !def->atom4->patch) {
    res4 = topo_mol_get_res(mol,&t4,def->rel4);
    if (!res4->reordered) a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if ( !a4 ) {
    t4.aname = def->atomstr4;
    a4 = topo_mol_get_atom(mol,&t4,def->rel4);
    if ( !a4 ) return -9;
  }
  if (!a1 || !a2 || !a3 || !a4) return -1;
  if (a1->isdrudlonepair || a2->isdrudlonepair || 
      a3->isdrudlonepair || a4->isdrudlonepair) return 0;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
  if ( ! tuple ) return -10;
  tuple->next = a1->impropers;
  tuple->atom[0] = a2;
  tuple->atom[1] = a3;
  tuple->atom[2] = a4;
  tuple->del = 0;
  a1->impropers = tuple;
#endif 


  return 0;
}

static void topo_mol_del_improper(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_improper_t *def) {
  topo_mol_improper_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  topo_mol_ident_t t1, t2, t3, t4;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL, *res4=NULL;
#endif

  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)

  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  for ( tuple = a1->impropers; tuple;
		tuple = topo_mol_improper_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 && tuple->atom[3] == a4 ) tuple->del = 1;
    if ( tuple->atom[0] == a4 && tuple->atom[1] == a3
	&& tuple->atom[2] == a2 && tuple->atom[3] == a1 ) tuple->del = 1;
  }
  
#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return;
  }

  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];

  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( ! a2 ) return;
  }


  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  
  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return;
  }

  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  
  if (def->atom4 && !def->atom4->patch) {
    res4 = topo_mol_get_res(mol,&t4,def->rel4);
    if (!res4->reordered) a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if ( !a4) {
    t4.aname = def->atomstr4;
    a4 = topo_mol_get_atom(mol,&t4,def->rel4);
    if ( !a4 ) return;
  }
  for ( tuple = a1->impropers; tuple; tuple = tuple->next ) {
    if (tuple->atom[0] == a2 && tuple->atom[1] == a3 && tuple->atom[2] == a4)
      tuple->del = 1;
  }
#endif
  
}

static int add_cmap_to_residues(topo_mol *mol, 
    const topo_mol_residue_t *resl[8], const char *anamel[8]) {
  int i;
  topo_mol_cmap_t *tuple;
  topo_mol_atom_t *al[8];

  if (! mol) return -1;
  for ( i=0; i<8; ++i ) {
    al[i] = topo_mol_get_atom_from_res(resl[i], anamel[i]);
    if (!al[i]) return -2-2*i;
  }
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_cmap_t));
  if ( ! tuple ) return -20;
  for ( i=0; i<8; ++i ) {
    tuple->next[i] = al[i]->cmaps;
    tuple->atom[i] = al[i];
  }
  for ( i=0; i<8; ++i ) {
    /* This must be in a separate loop because atoms may be repeated. */
    al[i]->cmaps = tuple;
  }
  tuple->del = 0;
  return 0;
}

static int topo_mol_add_cmap(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_cmap_t *def) {
  int i;
  topo_mol_cmap_t *tuple;
  topo_mol_atom_t *al[8];
  topo_mol_ident_t tl[8];
  if (! mol) return -1;
  for ( i=0; i<8; ++i ) {
    if ( def->resl[i] < 0 || def->resl[i] >= ntargets ) return -2-2*i;
    tl[i] = targets[def->resl[i]];
    tl[i].aname = def->atoml[i];
    al[i] = topo_mol_get_atom(mol,&tl[i],def->rell[i]);
    if ( ! al[i] ) return -3-2*i;
    
#if defined(NEWPSFGEN)
    if (al[i]->isdrudlonepair) continue;
#endif
  }
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_cmap_t));
  if ( ! tuple ) return -20;
  for ( i=0; i<8; ++i ) {
    tuple->next[i] = al[i]->cmaps;
    tuple->atom[i] = al[i];
  }
  for ( i=0; i<8; ++i ) {
    /* This must be in a separate loop because atoms may be repeated. */
    al[i]->cmaps = tuple;
  }
  tuple->del = 0;
  return 0;
}

static void topo_mol_del_cmap(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_cmap_t *def) {
  int i;
  topo_mol_cmap_t *tuple;
  topo_mol_atom_t *al[8];
  topo_mol_ident_t tl[8];
  if (! mol) return;
  for ( i=0; i<8; ++i ) {
    if ( def->resl[i] < 0 || def->resl[i] >= ntargets ) return;
    tl[i] = targets[def->resl[i]];
    tl[i].aname = def->atoml[i];
    al[i] = topo_mol_get_atom(mol,&tl[i],def->rell[i]);
    if ( ! al[i] ) return;
  }
  for ( tuple = al[0]->cmaps; tuple;
		tuple = topo_mol_cmap_next(tuple,al[0]) ) {
    int match1, match2;
    match1 = 0;
    for ( i=0; i<4 && (tuple->atom[i] == al[i]); ++i );
    if ( i == 4 ) match1 = 1;
    for ( i=0; i<4 && (tuple->atom[i] == al[3-i]); ++i );
    if ( i == 4 ) match1 = 1;
    match2 = 0;
    for ( i=0; i<4 && (tuple->atom[4+i] == al[4+i]); ++i );
    if ( i == 4 ) match2 = 1;
    for ( i=0; i<4 && (tuple->atom[4+i] == al[7-i]); ++i );
    if ( i == 4 ) match2 = 1;
    if ( match1 && match2 ) tuple->del = 1;
  }
}

#if !defined(NEWPSFGEN)
  
static int add_exclusion_to_residues(topo_mol *mol, 
    const topo_mol_residue_t *res1, const char *aname1,
    const topo_mol_residue_t *res2, const char *aname2) {
      
#else

static int add_exclusion_to_residues(topo_mol *mol, 
    const topo_mol_residue_t *res1, topo_defs_atom_t *atom1,
    const topo_mol_residue_t *res2, topo_defs_atom_t *atom2) {

#endif

  topo_mol_exclusion_t *tuple;
  topo_mol_atom_t *a1 = 0, *a2 = 0 ;

#if !defined(NEWPSFGEN)

  a1 = topo_mol_get_atom_from_res(res1, aname1);
  a2 = topo_mol_get_atom_from_res(res2, aname2);
  
#else
  if (!atom1 || !atom2) {
    return -1;
  }
  if (!res1->reordered) {
    a1 = res1->atomArray[atom1->atomIndex];
  } else {
    a1 = topo_mol_get_atom_from_res(res1, atom1->name);
  }
  
  if (!res2->reordered) {
    a2 = res2->atomArray[atom2->atomIndex];
  } else {
    a2 = topo_mol_get_atom_from_res(res2, atom2->name);
  }
  
#endif

  if (!a1 || !a2) return -1;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_exclusion_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->exclusions;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->exclusions;
  tuple->atom[1] = a2;
  tuple->del = 0;
  a1->exclusions = tuple;
  a2->exclusions = tuple;
  return 0;
}

#if !defined(NEWPSFGEN)

/* seems like these two functions were depricated */
static int topo_mol_add_exclusion(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_exclusion_t *def) {
  topo_mol_exclusion_t *tuple;
  topo_mol_atom_t *a1, *a2;
  topo_mol_ident_t t1, t2;
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_exclusion_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->exclusions;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->exclusions;
  tuple->atom[1] = a2;
  tuple->del = 0;
  a1->exclusions = tuple;
  a2->exclusions = tuple;
  return 0;
}

static void topo_mol_del_exclusion(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_exclusion_t *def) {
  topo_mol_exclusion_t *tuple;
  topo_mol_atom_t *a1, *a2;
  topo_mol_ident_t t1, t2;
  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  for ( tuple = a1->exclusions; tuple;
		tuple = topo_mol_exclusion_next(tuple,a1) ) {
    if ( tuple->atom[0] == a1 && tuple->atom[1] == a2 ) tuple->del = 1;
    if ( tuple->atom[0] == a2 && tuple->atom[1] == a1 ) tuple->del = 1;
  }
}

#endif

#if !defined(NEWPSFGEN)

static int add_conformation_to_residues(topo_mol *mol, 
    const topo_mol_residue_t *res1, const char *aname1,
    const topo_mol_residue_t *res2, const char *aname2,
    const topo_mol_residue_t *res3, const char *aname3,
    const topo_mol_residue_t *res4, const char *aname4, 
    topo_defs_conformation_t *def) {
#else

static int add_conformation_to_residues(topo_mol *mol,
    const topo_mol_ident_t *targets, 
    const topo_mol_residue_t *res1,
    const topo_mol_residue_t *res2,
    const topo_mol_residue_t *res3,
    const topo_mol_residue_t *res4, topo_defs_conformation_t *def) {

#endif

  topo_mol_conformation_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  
#if !defined(NEWPSFGEN)

  a1 = topo_mol_get_atom_from_res(res1, aname1);
  a2 = topo_mol_get_atom_from_res(res2, aname2);
  a3 = topo_mol_get_atom_from_res(res3, aname3);
  a4 = topo_mol_get_atom_from_res(res4, aname4);
  if (!a1 || !a2 || !a3 || !a4) return -1;
#else

  if (res1->reordered || !def->atom1 || def->atom1->patch) {
    a1 = topo_mol_get_atom_from_res(res1, def->atomstr1);
  } else {
    a1 = res1->atomArray[def->atom1->atomIndex];
  }
  
  if (res2->reordered || !def->atom2 || def->atom2->patch) {
    a2 = topo_mol_get_atom_from_res(res2, def->atomstr2);
  } else {
    a2 = res2->atomArray[def->atom2->atomIndex];
  }
  
  if (res3->reordered || !def->atom3 || def->atom3->patch) {
    a3 = topo_mol_get_atom_from_res(res3, def->atomstr3);
  } else {
    a3 = res3->atomArray[def->atom3->atomIndex];
  }
  
  if (res4->reordered ||  !def->atom4) {
    a4 = topo_mol_get_atom_from_res(res4, def->atomstr4);
  } else {
    a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if (!a1 || !a2 || !a3 || !a4) return -1;
  if (a1->isdrudlonepair || a2->isdrudlonepair || 
      a3->isdrudlonepair  || a4->isdrudlonepair) return 0;

#endif

  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_conformation_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->conformations;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->conformations;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->conformations;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->conformations;
  tuple->atom[3] = a4;
  tuple->del = 0;
  tuple->improper = def->improper;
  tuple->dist12 = def->dist12;
  tuple->angle123 = def->angle123;
  tuple->dihedral = def->dihedral;
  tuple->angle234 = def->angle234;
  tuple->dist34 = def->dist34;
  a1->conformations = tuple;
  a2->conformations = tuple;
  a3->conformations = tuple;
  a4->conformations = tuple;
  return 0;
}

static int topo_mol_add_conformation(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_conformation_t *def) {
  topo_mol_conformation_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  topo_mol_ident_t t1, t2, t3, t4;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1=NULL, *res2=NULL, *res3=NULL, *res4=NULL;
#endif
  
  if (! mol) return -1;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return -2;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)
  
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return -3;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( ! a2 ) return -5;
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( ! a3 ) return -7;
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);
  if ( ! a4 ) return -9;

#else

  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( ! a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( ! a1 ) return -3;
  }
  
  if ( def->res2 < 0 || def->res2 >= ntargets ) return -4;
  
  t2 = targets[def->res2];
  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return -5;
  }
  
  if ( def->res3 < 0 || def->res3 >= ntargets ) return -6;
  
  t3 = targets[def->res3];
  if (def->atom3 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( ! a3 ) return -7;
  }
  
  if ( def->res4 < 0 || def->res4 >= ntargets ) return -8;
  
  t4 = targets[def->res4];
  if (def->atom4 && !def->atom4->patch) {
    res4 = topo_mol_get_res(mol,&t4,def->rel4);
    if (!res4->reordered) a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if ( !a4) {
    t4.aname = def->atomstr4;
    a4 = topo_mol_get_atom(mol,&t4,def->rel4);
    if ( !a4 ) return -9;
  }
  if (!a1 || !a2 || !a3 || !a4) return -1;
  if (a1->isdrudlonepair || a2->isdrudlonepair || 
      a3->isdrudlonepair || a4->isdrudlonepair) return 0;

#endif

  tuple = memarena_alloc(mol->arena,sizeof(topo_mol_conformation_t));
  if ( ! tuple ) return -10;
  tuple->next[0] = a1->conformations;
  tuple->atom[0] = a1;
  tuple->next[1] = a2->conformations;
  tuple->atom[1] = a2;
  tuple->next[2] = a3->conformations;
  tuple->atom[2] = a3;
  tuple->next[3] = a4->conformations;
  tuple->atom[3] = a4;
  tuple->del = 0;
  tuple->improper = def->improper;
  tuple->dist12 = def->dist12;
  tuple->angle123 = def->angle123;
  tuple->dihedral = def->dihedral;
  tuple->angle234 = def->angle234;
  tuple->dist34 = def->dist34;
  a1->conformations = tuple;
  a2->conformations = tuple;
  a3->conformations = tuple;
  a4->conformations = tuple;
  return 0;
}

static void topo_mol_del_conformation(topo_mol *mol, const topo_mol_ident_t *targets,
				int ntargets, topo_defs_conformation_t *def) {
  topo_mol_conformation_t *tuple;
  topo_mol_atom_t *a1 =0, *a2 =0, *a3 =0, *a4 =0;
  topo_mol_ident_t t1, t2, t3, t4;
#if defined(NEWPSFGEN)
  topo_mol_residue_t *res1 =0, *res2 =0, *res3 =0, *res4 =0;
#endif

  if (! mol) return;
  if ( def->res1 < 0 || def->res1 >= ntargets ) return;
  t1 = targets[def->res1];
  
#if !defined(NEWPSFGEN)
  
  t1.aname = def->atom1;
  a1 = topo_mol_get_atom(mol,&t1,def->rel1);
  if ( ! a1 ) return;
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  t2 = targets[def->res2];
  t2.aname = def->atom2;
  a2 = topo_mol_get_atom(mol,&t2,def->rel2);
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  t3 = targets[def->res3];
  t3.aname = def->atom3;
  a3 = topo_mol_get_atom(mol,&t3,def->rel3);
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  t4 = targets[def->res4];
  t4.aname = def->atom4;
  a4 = topo_mol_get_atom(mol,&t4,def->rel4);

#else
  
  if (def->atom1 && !def->atom1->patch) {
    res1 = topo_mol_get_res(mol,&t1,def->rel1);
    if (!res1->reordered) a1 = res1->atomArray[def->atom1->atomIndex];
  }
  if ( !a1 ) {
    t1.aname = def->atomstr1;
    a1 = topo_mol_get_atom(mol,&t1,def->rel1);
    if ( !a1 ) return;
  }
  
  if ( def->res2 < 0 || def->res2 >= ntargets ) return;
  
  t2 = targets[def->res2];
  if (def->atom2 && !def->atom2->patch) {
    res2 = topo_mol_get_res(mol,&t2,def->rel2);
    if (!res2->reordered) a2 = res2->atomArray[def->atom2->atomIndex];
  }
  if ( !a2 ) {
    t2.aname = def->atomstr2;
    a2 = topo_mol_get_atom(mol,&t2,def->rel2);
    if ( !a2 ) return ;
  }
  
  if ( def->res3 < 0 || def->res3 >= ntargets ) return;
  
  t3 = targets[def->res3];
  if (def->atom4 && !def->atom3->patch) {
    res3 = topo_mol_get_res(mol,&t3,def->rel3);
    if (!res3->reordered) a3 = res3->atomArray[def->atom3->atomIndex];
  }
  if ( !a3 ) {
    t3.aname = def->atomstr3;
    a3 = topo_mol_get_atom(mol,&t3,def->rel3);
    if ( !a3 ) return;
  }
  
  if ( def->res4 < 0 || def->res4 >= ntargets ) return;
  
  t4 = targets[def->res4];
  if (def->atom4 && !def->atom4->patch) {
    res4 = topo_mol_get_res(mol,&t4,def->rel4);
    if (!res4->reordered) a4 = res4->atomArray[def->atom4->atomIndex];
  }
  if ( !a4) {
    t4.aname = def->atomstr4;
    a4 = topo_mol_get_atom(mol,&t4,def->rel4);
    if ( !a4 ) return;
  }
  
#endif

  for ( tuple = a1->conformations; tuple;
		tuple = topo_mol_conformation_next(tuple,a1) ) {
    if ( tuple->improper == def->improper
	&&  tuple->atom[0] == a1 && tuple->atom[1] == a2
	&& tuple->atom[2] == a3 && tuple->atom[3] == a4 ) tuple->del = 1;
    if ( tuple->improper == def->improper
	&& tuple->atom[0] == a4 && tuple->atom[1] == a3
	&& tuple->atom[2] == a2 && tuple->atom[3] == a1 ) tuple->del = 1;
  }
}

static int topo_mol_auto_angles(topo_mol *mol, topo_mol_segment_t *segp);
static int topo_mol_auto_dihedrals(topo_mol *mol, topo_mol_segment_t *segp);

int topo_mol_end(topo_mol *mol) {
  int i,n;
  int idef;
  topo_defs *defs;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_defs_residue_t *resdef;
  topo_defs_atom_t *atomdef;
  topo_defs_bond_t *bonddef;
  topo_defs_angle_t *angldef;
  topo_defs_dihedral_t *dihedef;
  topo_defs_improper_t *imprdef;
  topo_defs_cmap_t *cmapdef;
  topo_defs_exclusion_t *excldef;
  topo_defs_conformation_t *confdef;
  topo_mol_ident_t target;
  char errmsg[128];
  int firstdefault=0, lastdefault=0;

  if ( ! mol ) return -1;
  if ( ! mol->buildseg ) {
    topo_mol_log_error(mol,"no segment in progress for end");
    return -1;
  }
  seg = mol->buildseg;
  mol->buildseg = 0;
  defs = mol->defs;

  /* add atoms */
  n = hasharray_count(seg->residue_hash);
  for ( i=0; i<n; ++i ) {
    res = &(seg->residue_array[i]);
    idef = hasharray_index(defs->residue_hash,res->name);
    if ( idef == HASHARRAY_FAIL ) {
      sprintf(errmsg,"unknown residue type %s",res->name);
      topo_mol_log_error(mol,errmsg);
      return -1;
    }
    resdef = &(mol->defs->residue_array[idef]);
    if ( resdef->patch ) {
      sprintf(errmsg,"unknown residue type %s",res->name);
      topo_mol_log_error(mol,errmsg);
      return -1;
    }

    /* patches */
    if ( i==0 && ! strlen(seg->pfirst) ) {
      strcpy(seg->pfirst,resdef->pfirst);
      firstdefault = 1;
    }
    if ( i==(n-1) && ! strlen(seg->plast) ) {
      strcpy(seg->plast,resdef->plast);
      lastdefault = 1;
    }
    
#if !defined(NEWPSFGEN)

    for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
      if ( topo_mol_add_atom(mol,&(res->atoms),0,atomdef) ) { 
        sprintf(errmsg,"add atom failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
        return -8;
      }
    }

#else
  
    // add +1 to add a null pointer at the end of the array
    res->atomArray = 
    (topo_mol_atom_t**) calloc(1, (resdef->atomNum + 1)*sizeof(topo_mol_atom_t*));
    res->atomSize = 0;
    res->reordered = 0;
    res->pres[0] = '\0';
    res->lonepairs = 0;
    res->numaniso = 0;
    res->aniso = 0;
    for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
      if ( topo_mol_add_atom(mol,0,atomdef,res) ) { 
        sprintf(errmsg,"add atom failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
        return -9;
      }
    }
    /* the topo_mol_add_atom changes the res->haslonepairs to 1 of the residues
     * contains lonepairs
    */
    if (resdef->lonepairs) {
      res->lonepairs = 1;
    }

#endif

  }  


  for ( i=0; i<n; ++i ) {
    res = &(seg->residue_array[i]);
    idef = hasharray_index(defs->residue_hash,res->name);
    if ( idef == HASHARRAY_FAIL ) {
      sprintf(errmsg,"unknown residue type %s",res->name);
      topo_mol_log_error(mol,errmsg);
      return -1;
    }
    resdef = &(mol->defs->residue_array[idef]);
    target.segid = seg->segid;
    target.resid = res->resid;
    for ( bonddef = resdef->bonds; bonddef; bonddef = bonddef->next ) {
      int ires1, ires2;
      if (bonddef->res1 != 0 || bonddef->res2 != 0) {
        /* 
         * XXX This should be caught much earlier, like when the topology
         * file is initially read in. 
         */
#if !defined(NEWPSFGEN)
  
        sprintf(errmsg, "ERROR: Bad bond definition %s %s-%s; skipping.",
            res->name, bonddef->atom1, bonddef->atom2);
#else

        sprintf(errmsg, "ERROR: Bad bond definition %s %s-%s; skipping.",\
            res->name, bonddef->atomstr1, bonddef->atomstr2);

#endif        
            
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      ires1=bonddef->rel1+i;
      ires2=bonddef->rel2+i;
      if (ires1 < 0 || ires2 < 0 || ires1 >= n || ires2 >= n) {
        
#if !defined(NEWPSFGEN)

        sprintf(errmsg, "Info: skipping bond %s-%s at %s of segment.", 
            bonddef->atom1, bonddef->atom2, i==0 ? "beginning" : "end");
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      if (add_bond_to_residues(mol, 
            &(seg->residue_array[ires1]), bonddef->atom1,
            &(seg->residue_array[ires2]), bonddef->atom2)) {
        sprintf(errmsg, 
            "ERROR: Missing atoms for bond %s(%d) %s(%d) in residue %s:%s",
            bonddef->atom1,bonddef->rel1,bonddef->atom2,bonddef->rel2,
            res->name,res->resid);
        topo_mol_log_error(mol, errmsg);
      }
#else
        sprintf(errmsg, "Info: skipping bond %s-%s at %s of segment.", 
            bonddef->atomstr1, bonddef->atomstr2, i==0 ? "beginning" : "end");
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      if (add_bond_to_residues(mol, &target, 
            &(seg->residue_array[ires1]), 
            &(seg->residue_array[ires2]), bonddef)) {
        sprintf(errmsg, 
            "ERROR: Missing atoms for bond %s(%d) %s(%d) in residue %s:%s",
            bonddef->atomstr1,bonddef->rel1,bonddef->atomstr2,bonddef->rel2,
            res->name,res->resid);
        topo_mol_log_error(mol, errmsg);
      }  
    
#endif

    }
    if ( seg->auto_angles && resdef->angles ) {
      sprintf(errmsg,"Warning: explicit angles in residue %s:%s will be deleted during autogeneration",res->name,res->resid);
      topo_mol_log_error(mol,errmsg);
    }
    for ( angldef = resdef->angles; angldef; angldef = angldef->next ) {
      if ( topo_mol_add_angle(mol,&target,1,angldef) ) {
        sprintf(errmsg,"Warning: add angle failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
    if ( seg->auto_dihedrals && resdef->dihedrals) {
      sprintf(errmsg,"Warning: explicit dihedrals in residue %s:%s will be deleted during autogeneration",res->name,res->resid);
      topo_mol_log_error(mol,errmsg);
    }
    for ( dihedef = resdef->dihedrals; dihedef; dihedef = dihedef->next ) {
      if ( topo_mol_add_dihedral(mol,&target,1,dihedef) ) {
        sprintf(errmsg,"Warning: add dihedral failed in residue %s:%s",res->name,res->resid);
        topo_mol_log_error(mol,errmsg);
      }
    }
    for ( imprdef = resdef->impropers; imprdef; imprdef = imprdef->next ) {
      int ires1, ires2, ires3, ires4;
      if (imprdef->res1 != 0 || imprdef->res2 != 0 || imprdef->res3 != 0 ||
          imprdef->res4 != 0) {
            
#if !defined(NEWPSFGEN)

        sprintf(errmsg, "ERROR: Bad improper definition %s %s-%s-%s-%s; skipping.",
            res->name, imprdef->atom1, imprdef->atom2, imprdef->atom3, 
            imprdef->atom4);

#else

        sprintf(errmsg, "ERROR: Bad improper definition %s %s-%s-%s-%s; skipping.",
            res->name, imprdef->atomstr1, imprdef->atomstr2, imprdef->atomstr3, 
            imprdef->atomstr4);
    
#endif
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      ires1=imprdef->rel1+i;
      ires2=imprdef->rel2+i;
      ires3=imprdef->rel3+i;
      ires4=imprdef->rel4+i;
      if (ires1 < 0 || ires2 < 0 || ires3 < 0 || ires4 < 0 ||
          ires1 >= n || ires2 >= n || ires3 >= n || ires4 >= n) {
            
#if !defined(NEWPSFGEN)

        sprintf(errmsg,"Info: skipping improper %s-%s-%s-%s at %s of segment.", 
            imprdef->atom1, imprdef->atom2, imprdef->atom3, imprdef->atom4,
            i==0 ? "beginning" : "end");

#else

        sprintf(errmsg,"Info: skipping improper %s-%s-%s-%s at %s of segment.", 
            imprdef->atomstr1, imprdef->atomstr2, imprdef->atomstr3, 
            imprdef->atomstr4, i==0 ? "beginning" : "end");
    
#endif            
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      
#if !defined(NEWPSFGEN)

      if (add_improper_to_residues(mol, 
            &(seg->residue_array[ires1]), imprdef->atom1,
            &(seg->residue_array[ires2]), imprdef->atom2,
            &(seg->residue_array[ires3]), imprdef->atom3,
            &(seg->residue_array[ires4]), imprdef->atom4)) {
        sprintf(errmsg, 
            "ERROR: Missing atoms for improper %s(%d) %s(%d) %s(%d) %s(%d)\n\tin residue %s:%s",
            imprdef->atom1,imprdef->rel1,imprdef->atom2,imprdef->rel2,
            imprdef->atom3,imprdef->rel3,imprdef->atom4,imprdef->rel4,
            res->name,res->resid);
        topo_mol_log_error(mol, errmsg);
      }

#else
  
      if (add_improper_to_residues(mol,&target, 
            &(seg->residue_array[ires1]),
            &(seg->residue_array[ires2]),
            &(seg->residue_array[ires3]),
            &(seg->residue_array[ires4]), imprdef)) {
        sprintf(errmsg, 
            "ERROR: Missing atoms for improper %s(%d) %s(%d) %s(%d) %s(%d)\n\tin residue %s:%s",
            imprdef->atomstr1,imprdef->rel1,imprdef->atomstr2,imprdef->rel2,
            imprdef->atomstr3,imprdef->rel3,imprdef->atomstr4,imprdef->rel4,
            res->name,res->resid);
        topo_mol_log_error(mol, errmsg);
      }
      
#endif

    }

    for ( cmapdef = resdef->cmaps; cmapdef; cmapdef = cmapdef->next ) {
      int j, iresl[8];
      const topo_mol_residue_t *resl[8];
      const char *atoml[8];
      for ( j=0; j<8 && (cmapdef->resl[j] == 0); ++j );
      if ( j != 8 ) {
        sprintf(errmsg, "ERROR: Bad cross-term definition %s %s-%s-%s-%s-%s-%s-%s-%s; skipping.",
            res->name, cmapdef->atoml[0], cmapdef->atoml[1],
            cmapdef->atoml[2], cmapdef->atoml[3], cmapdef->atoml[4],
            cmapdef->atoml[5], cmapdef->atoml[6], cmapdef->atoml[7]);
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      for ( j=0; j<8; ++j ) {
        iresl[j] = cmapdef->rell[j]+i;
      }
      for ( j=0; j<8 && (iresl[j] >= 0) && (iresl[j] < n); ++j );
      if ( j != 8 ) {
        sprintf(errmsg,"Info: skipping cross-term %s-%s-%s-%s-%s-%s-%s-%s at %s of segment.", 
            cmapdef->atoml[0], cmapdef->atoml[1],
            cmapdef->atoml[2], cmapdef->atoml[3], cmapdef->atoml[4],
            cmapdef->atoml[5], cmapdef->atoml[6], cmapdef->atoml[7],
            i==0 ? "beginning" : "end");
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      for ( j=0; j<8; ++j ) {
        resl[j] = &seg->residue_array[iresl[j]];
        atoml[j] = cmapdef->atoml[j];
      }
      if (add_cmap_to_residues(mol, resl, atoml) ) {
        sprintf(errmsg, 
            "ERROR: Missing atoms for cross-term  %s(%d) %s(%d) %s(%d) %s(%d) %s(%d) %s(%d) %s(%d) %s(%d)\n\tin residue %s:%s",
            cmapdef->atoml[0],cmapdef->rell[0],
            cmapdef->atoml[1],cmapdef->rell[1],
            cmapdef->atoml[2],cmapdef->rell[2],
            cmapdef->atoml[3],cmapdef->rell[3],
            cmapdef->atoml[4],cmapdef->rell[4],
            cmapdef->atoml[5],cmapdef->rell[5],
            cmapdef->atoml[6],cmapdef->rell[6],
            cmapdef->atoml[7],cmapdef->rell[7],
            res->name,res->resid);
        topo_mol_log_error(mol, errmsg);
      }
    }
    for ( excldef = resdef->exclusions; excldef; excldef = excldef->next ) {
      int ires1, ires2;
      if (excldef->res1 != 0 || excldef->res2 != 0) {
        
#if !defined(NEWPSFGEN)

        sprintf(errmsg, "ERROR: Bad exclusion definition %s %s-%s; skipping.",
            res->name, excldef->atom1, excldef->atom2);
            
#else 

        sprintf(errmsg, "ERROR: Bad exclusion definition %s %s-%s; skipping.",
            res->name, excldef->atom1->name, excldef->atom2->name);
            
#endif

        topo_mol_log_error(mol, errmsg);
        continue;
      }
      ires1=excldef->rel1+i;
      ires2=excldef->rel2+i;
      if (ires1 < 0 || ires2 < 0 || ires1 >= n || ires2 >= n) {

#if !defined(NEWPSFGEN)

        sprintf(errmsg, "Info: skipping exclusion %s-%s at %s of segment.", 
            excldef->atom1, excldef->atom2, i==0 ? "beginning" : "end");
            
#else

        sprintf(errmsg, "Info: skipping exclusion %s-%s at %s of segment.", 
            excldef->atom1->name, excldef->atom2->name, i==0 ? "beginning" : "end");
            
#endif
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      if (add_exclusion_to_residues(mol, 
            &(seg->residue_array[ires1]), excldef->atom1,
            &(seg->residue_array[ires2]), excldef->atom2)) {
              
#if !defined(NEWPSFGEN)
  
        sprintf(errmsg, 
            "ERROR: Missing atoms for exclusion %s(%d) %s(%d) in residue %s:%s",
            excldef->atom1,excldef->rel1,excldef->atom2,excldef->rel2,
            res->name,res->resid);

#else

        sprintf(errmsg, 
            "ERROR: Missing atoms for exclusion %s(%d) %s(%d) in residue %s:%s",
            excldef->atom1->name,excldef->rel1,excldef->atom2->name,excldef->rel2,
            res->name,res->resid);

#endif    
        
        topo_mol_log_error(mol, errmsg);
      }
    }

    for ( confdef = resdef->conformations; confdef; confdef = confdef->next ) {
      int ires1, ires2, ires3, ires4;
      if (confdef->res1 != 0 || confdef->res2 != 0 || confdef->res3 != 0 ||
          confdef->res4 != 0) {
            
#if !defined(NEWPSFGEN)

        sprintf(errmsg, "ERROR: Bad conformation definition %s %s-%s-%s-%s; skipping.",
            res->name, confdef->atom1, confdef->atom2, confdef->atom3, 
            confdef->atom4);

#else

        sprintf(errmsg, "ERROR: Bad conformation definition %s %s-%s-%s-%s; skipping.",
            res->name, confdef->atomstr1, confdef->atomstr2, confdef->atomstr3, 
            confdef->atomstr4);

#endif   
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      ires1=confdef->rel1+i;
      ires2=confdef->rel2+i;
      ires3=confdef->rel3+i;
      ires4=confdef->rel4+i;
    
            
#if !defined(NEWPSFGEN)

      if (ires1 < 0 || ires2 < 0 || ires3 < 0 || ires4 < 0 ||
          ires1 >= n || ires2 >= n || ires3 >= n || ires4 >= n) {
        sprintf(errmsg,"Info: skipping conformation %s-%s-%s-%s at %s of segment.", 
            confdef->atom1, confdef->atom2, confdef->atom3, confdef->atom4,
            i==0 ? "beginning" : "end");
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      if (add_conformation_to_residues(mol, 
            &(seg->residue_array[ires1]), confdef->atom1,
            &(seg->residue_array[ires2]), confdef->atom2,
            &(seg->residue_array[ires3]), confdef->atom3,
            &(seg->residue_array[ires4]), confdef->atom4, confdef)) {
        sprintf(errmsg, "Warning: missing atoms for conformation %s %s-%s-%s-%s; skipping.",
            res->name, confdef->atom1, confdef->atom2, confdef->atom3, 
            confdef->atom4);
        topo_mol_log_error(mol, errmsg);
      }
      
#else

      if (ires1 < 0 || ires2 < 0 || ires3 < 0 || ires4 < 0 ||
          ires1 >= n || ires2 >= n || ires3 >= n || ires4 >= n) {
        sprintf(errmsg,"Info: skipping conformation %s-%s-%s-%s at %s of segment.", 
            confdef->atomstr1, confdef->atomstr2, confdef->atomstr3, 
            confdef->atomstr4, i==0 ? "beginning" : "end");
        topo_mol_log_error(mol, errmsg);
        continue;
      }
      if (add_conformation_to_residues(mol,&target, 
            &(seg->residue_array[ires1]),
            &(seg->residue_array[ires2]),
            &(seg->residue_array[ires3]),
            &(seg->residue_array[ires4]), confdef)) {
        sprintf(errmsg, "Warning: missing atoms for conformation %s %s-%s-%s-%s; skipping.",
            res->name, confdef->atomstr1, confdef->atomstr2, confdef->atomstr3, 
            confdef->atomstr4);
        topo_mol_log_error(mol, errmsg);
      }
      
#endif 
  
    }
  }
  
#if defined(NEWPSFGEN)
  {
    int ires, nres;
    topo_defs_anisotropy_t *aniso;
    nres = hasharray_count(seg->residue_hash);
  
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      /* assign the lonepairs hosts and coordinates*/
      idef = hasharray_index(defs->residue_hash,res->name);
      resdef = &(mol->defs->residue_array[idef]);
    
      if (res->lonepairs) {
        for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
          if (atomdef->islonepair) {
            if (topo_mol_update_lonepair(mol,atomdef,res, 0, 0) ) {
              sprintf(errmsg,"failed to update lonepair in residue %s:%s",
              res->name,res->resid);
              topo_mol_log_error(mol,errmsg);
              return -2;
            }
          }
        }
      }

      /* update the anisotropy fields of the residue */
      for ( aniso = resdef->aniso; aniso; aniso = aniso->next ) {
          if (topo_mol_add_anisotropy(mol, aniso, res, &target) ) {
          sprintf(errmsg,"failed to add anisotropy in residue %s:%s",
          res->name,res->resid);
          topo_mol_log_error(mol,errmsg);
          return -2;
        } 
      }  
    }
  }  
#endif

  /* apply patches, last then first because dipeptide patch ACED depends on CT3 atom NT */

  res = &(seg->residue_array[n-1]);
  if ( ! strlen(seg->plast) ) strcpy(seg->plast,"NONE");

  target.segid = seg->segid;
  target.resid = res->resid;
  if ( topo_mol_patch(mol, &target, 1, seg->plast, 0,
	seg->auto_angles, seg->auto_dihedrals, lastdefault) ) return -10;

  res = &(seg->residue_array[0]);
  if ( ! strlen(seg->pfirst) ) strcpy(seg->pfirst,"NONE");

  target.segid = seg->segid;
  target.resid = res->resid;
  if ( topo_mol_patch(mol, &target, 1, seg->pfirst, 1,
	seg->auto_angles, seg->auto_dihedrals, firstdefault) ) return -11;

#if defined(NEWPSFGEN) && 0
  int ires, nres, atomid, iseg, nseg;
  nres = hasharray_count(seg->residue_hash);
  
  /* update the atomid to ease the guess of the angles and dihedrals */
  atomid = 0;
  for ( ires=0; ires<nres; ++ires ) {
    res = &(seg->residue_array[ires]);
    for ( i = 0; i < res->atomSize; i++ ) {
      res->atomArray[i]->atomid = ++atomid;
    }
  }
  
#endif
  
  if (seg->auto_angles && topo_mol_auto_angles(mol, seg)) return -12;
  if (seg->auto_dihedrals && topo_mol_auto_dihedrals(mol, seg)) return -13;

  return 0;
}


int topo_mol_regenerate_resids(topo_mol *mol) {
  int ires, nres, iseg, nseg, npres;
  int prevresid, resid, npatchresptrs, ipatch;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_patch_t **patchptr, *patch;
  topo_mol_patchres_t *patchres, **patchresptrs;
  char newresid[NAMEMAXLEN+20], (*newpatchresids)[NAMEMAXLEN];

  if (! mol) return -1;

  nseg = hasharray_count(mol->segment_hash);
  npatchresptrs=0;

  /* clean patches so only valid items remain */
  for ( patchptr = &(mol->patches); *patchptr; ) {
    npres=0;
    for ( patchres = (*patchptr)->patchresids; patchres; patchres = patchres->next ) {
      ++npres;
      /* Test the existence of segid:resid for the patch */
      if (!topo_mol_validate_patchres(mol,(*patchptr)->pname,patchres->segid, patchres->resid)) {
        break;
      }
    }
    if ( patchres ) {  /* remove patch from list */
      *patchptr = (*patchptr)->next;
      continue;
    }
    npatchresptrs += npres;
    patchptr = &((*patchptr)->next);  /* continue to next patch */
  }

  patchresptrs = malloc(npatchresptrs * sizeof(topo_mol_patchres_t*));
  if ( ! patchresptrs ) return -5;
  newpatchresids = calloc(npatchresptrs, NAMEMAXLEN);
  if ( ! newpatchresids ) return -6;

  for ( ipatch=0, patch = mol->patches; patch; patch = patch->next ) {
    for ( patchres = patch->patchresids; patchres; patchres = patchres->next ) {
      patchresptrs[ipatch++] = patchres;
    }
  }

  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if ( ! seg ) continue;
    nres = hasharray_count(seg->residue_hash);
    if ( hasharray_clear(seg->residue_hash) == HASHARRAY_FAIL ) return -2;

    prevresid = -100000;
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      resid = atoi(res->resid);
      if ( resid <= prevresid ) resid = prevresid + 1;
      sprintf(newresid, "%d", resid);
      if ( NAMETOOLONG(newresid) ) return -3;
      if ( strcmp(res->resid, newresid) ) { /* changed, need to check patches */
        for ( ipatch=0; ipatch < npatchresptrs; ++ipatch ) {
          if ( ( ! strcmp(seg->segid, patchresptrs[ipatch]->segid) ) &&
               ( ! strcmp(res->resid, patchresptrs[ipatch]->resid) ) ) {
            sprintf(newpatchresids[ipatch], "%d", resid);
          }
        }
      }
      sprintf(res->resid, "%d", resid);
      if ( hasharray_reinsert(seg->residue_hash,res->resid,ires) != ires ) return -4;
      prevresid = resid;
    }
  }

  for ( ipatch=0; ipatch < npatchresptrs; ++ipatch ) {
    if ( newpatchresids[ipatch][0] ) {
      strcpy(patchresptrs[ipatch]->resid,newpatchresids[ipatch]);
    }
  }

  free(patchresptrs);
  free(newpatchresids);
  return 0;
}

int topo_mol_regenerate_angles(topo_mol *mol) {
  int errval;
  if ( mol ) {
    memarena_destroy(mol->angle_arena);
    mol->angle_arena = memarena_create();
  }
  errval = topo_mol_auto_angles(mol,0);
  if ( errval ) {
    char errmsg[128];
    sprintf(errmsg,"Error code %d",errval);
    topo_mol_log_error(mol,errmsg);
  }
  return errval;
}

int topo_mol_regenerate_dihedrals(topo_mol *mol) {
  int errval;
  if ( mol ) {
    memarena_destroy(mol->dihedral_arena);
    mol->dihedral_arena = memarena_create();
  }
  errval = topo_mol_auto_dihedrals(mol,0);
  if ( errval ) {
    char errmsg[128];
    sprintf(errmsg,"Error code %d",errval);
    topo_mol_log_error(mol,errmsg);
  }
  return errval;
}

#if !defined(NEWPSFGEN)

static int is_hydrogen(topo_mol_atom_t *atom) {
  return ( atom->mass < 3.5 && atom->name[0] == 'H' );
}

static int is_oxygen(topo_mol_atom_t *atom) {
  return ( atom->mass > 14.5 && atom->mass < 18.5 && atom->name[0] == 'O' );
}

#else
int is_hydrogen(topo_mol_atom_t *atom) {
  int i = 0;
  while (isspace(atom->name[i]) && i < NAMEMAXLEN ) i++;
  return ( atom->mass < 3.5 && atom->name[i] == 'H' );
}

int is_oxygen(topo_mol_atom_t *atom) {
  int i = 0;
  while (isspace(atom->name[i]) && i < NAMEMAXLEN ) i++;
  return ( atom->mass > 14.5 && atom->mass < 18.5 && atom->name[i] == 'O' );
}

#endif


static int topo_mol_auto_angles(topo_mol *mol, topo_mol_segment_t *segp) {
  int ires, nres, iseg, nseg;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_bond_t *b1, *b2;
  topo_mol_angle_t *tuple;
  topo_mol_atom_t *atom=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  
#if defined(NEWPSFGEN)  
  int i; /* atom counter in the for loop */
  int atomid = 0;
#endif 
  
  if (! mol) return -1;
  nseg = segp ? 1 : hasharray_count(mol->segment_hash);

  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = segp ? segp : mol->segment_array[iseg];
    if ( ! seg ) continue;

    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      
#if !defined(NEWPSFGEN)  
      for ( atom = res->atoms; atom; atom = atom->next ) {
#else
      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        atom->atomid = ++atomid;
#endif

        if ( ! segp ) { atom->angles = NULL; }
        for ( tuple = atom->angles; tuple;
		tuple = topo_mol_angle_next(tuple,atom) ) {
          tuple->del = 1;
        }
      }
    }
  }

  for ( iseg=0; iseg<nseg; ++iseg ) {
  seg = segp ? segp : mol->segment_array[iseg];
  if ( ! seg ) continue;

  nres = hasharray_count(seg->residue_hash);
  for ( ires=0; ires<nres; ++ires ) {
    res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)

    for ( atom = res->atoms; atom; atom = atom->next ) {
      a2 = atom;
      for ( b1 = atom->bonds; b1; b1 = topo_mol_bond_next(b1,atom) ) {
        if ( b1->del ) continue;
        if ( b1->atom[0] == atom ) a1 = b1->atom[1];
        else if ( b1->atom[1] == atom ) a1 = b1->atom[0];
        else return -5;
        b2 = b1;  while ( (b2 = topo_mol_bond_next(b2,atom)) ) {
          if ( b2->del ) continue;
          if ( b2->atom[0] == atom ) a3 = b2->atom[1];
          else if ( b2->atom[1] == atom ) a3 = b2->atom[0];
          else return -6;
          if ( is_hydrogen(a2) && ( ! topo_mol_bond_next(b2,atom) ) &&
               ( ( is_hydrogen(a1) && is_oxygen(a3) ) ||
                 ( is_hydrogen(a3) && is_oxygen(a1) ) ) )
            continue;  /* extra H-H bond on water */
          tuple = memarena_alloc(mol->angle_arena,sizeof(topo_mol_angle_t));
          if ( ! tuple ) return -10;
          tuple->next[0] = a1->angles;
          tuple->atom[0] = a1;
          tuple->next[1] = a2->angles;
          tuple->atom[1] = a2;
          tuple->next[2] = a3->angles;
          tuple->atom[2] = a3;
          tuple->del = 0;
          a1->angles = tuple;
          a2->angles = tuple;
          a3->angles = tuple;
        }
      }
    }
#else

    for ( i = 0; i < res->atomSize; i++ ) {
      atom = res->atomArray[i];
      a2 = atom;
      for ( b1 = atom->bonds; b1; b1 = topo_mol_bond_next(b1,atom) ) {
        if ( b1->del ) continue;
        if ( b1->atom[0] == atom ) a1 = b1->atom[1];
        else if ( b1->atom[1] == atom ) a1 = b1->atom[0];
        else return -5;
        b2 = b1;  while ( (b2 = topo_mol_bond_next(b2,atom)) ) {
          if ( b2->del ) continue;
          if ( b2->atom[0] == atom ) a3 = b2->atom[1];
          else if ( b2->atom[1] == atom ) a3 = b2->atom[0];
          if ( is_hydrogen(a2) && ( ! topo_mol_bond_next(b2,atom) ) &&
               ( ( is_hydrogen(a1) && is_oxygen(a3) ) ||
                 ( is_hydrogen(a3) && is_oxygen(a1) ) ) )
            continue;  /* extra H-H bond on water */
          tuple = memarena_alloc(mol->angle_arena,sizeof(topo_mol_angle_t));
          if ( ! tuple ) return -10;
          tuple->next[0] = a1->angles;
          tuple->atom[0] = a1;
          tuple->next[1] = a2->angles;
          tuple->atom[1] = a2;
          tuple->next[2] = a3->angles;
          tuple->atom[2] = a3;
          tuple->del = 0;
          a1->angles = tuple;
          a2->angles = tuple;
          a3->angles = tuple;
        }
      }
    }
#endif
  }
  }

  return 0;
}

static int topo_mol_auto_dihedrals(topo_mol *mol, topo_mol_segment_t *segp) {
  int ires, nres, iseg, nseg, atomid;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_dihedral_t *tuple;
  topo_mol_atom_t *atom=NULL, *a1=NULL, *a2=NULL, *a3=NULL, *a4=NULL;
  
#if !defined(NEWPSFGEN)
  int found, count1, count2;  
  topo_mol_angle_t *g1, *g2;
#else
  int i; /* atom counter in the for loop */
  topo_mol_bond_t *bond1=NULL, *bond2=NULL, *bond3=NULL; /* bond to find dihedrals */
  atomid = 0;
#endif 
  
  if (! mol) return -1;
  nseg = segp ? 1 : hasharray_count(mol->segment_hash);
  
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = segp ? segp : mol->segment_array[iseg];
    if ( ! seg ) continue;

    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      
#if !defined(NEWPSFGEN)

      for ( atom = res->atoms; atom; atom = atom->next ) {
        if ( ! segp ) { atom->dihedrals = NULL; }
        for ( tuple = atom->dihedrals; tuple;
    tuple = topo_mol_dihedral_next(tuple,atom) ) {
          tuple->del = 1;
        }
      }    
#else

      for ( i = 0; i < res->atomSize; i++ ) {
        /* Although the topo_mol_auto_angles might already had reset the 
         * atomid, it is faster to reset it that check if for every atom
        */
        res->atomArray[i]->atomid = ++atomid;
        if ( ! segp ) { res->atomArray[i]->dihedrals = NULL; }
        for ( tuple = res->atomArray[i]->dihedrals; tuple;
            tuple = tuple->next) {
          tuple->del = 1;
        }
      }
          
#endif
        
    }
  }

#if !defined(NEWPSFGEN)
  /*  number atoms, needed to avoid duplicate dihedrals below  */
  /*  assumes no inter-segment bonds if segp is non-null  */
  atomid = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = segp ? segp : mol->segment_array[iseg];
    if ( ! seg ) continue;

    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);



    for ( atom = res->atoms; atom; atom = atom->next ) {
      atom->atomid = ++atomid;
    }
    }
  }
  count1 = count2 = 0;
#else

  /* in the NEWPSFGEN this step is done in the topo_mol_end().
   * before the topo_mol_auto_angles and topo_auto_dihedrals*/
  
#endif

 
  for ( iseg=0; iseg<nseg; ++iseg ) {
  seg = segp ? segp : mol->segment_array[iseg];
  if ( ! seg ) continue;

  nres = hasharray_count(seg->residue_hash);
  for ( ires=0; ires<nres; ++ires ) {
    res = &(seg->residue_array[ires]);
    
#if !defined(NEWPSFGEN)  

    for ( atom = res->atoms; atom; atom = atom->next ) {
      for ( g1 = atom->angles; g1; g1 = topo_mol_angle_next(g1,atom) ) {
        if ( g1->del ) continue;
        if ( g1->atom[1] != atom ) continue;
        for ( g2 = atom->angles; g2; g2 = topo_mol_angle_next(g2,atom) ) {
          if ( g2->del ) continue;
          if ( g2->atom[1] == atom ) continue;
          found = 0;
          if ( g2->atom[0] == atom ) {  /*  XBX BXX  */
            if ( g2->atom[1] == g1->atom[0] ) {  /*  CBA BCD  */
              a1 = g1->atom[2];
              a2 = g1->atom[1];  /* == g2->atom[0] */
              a3 = g1->atom[0];  /* == g2->atom[1] */
              a4 = g2->atom[2];
              found = ( a1->atomid < a4->atomid );
              if ( a1 != a4 ) ++count2;
            } else if ( g2->atom[1] == g1->atom[2] ) {  /*  ABC BCD  */
              a1 = g1->atom[0];
              a2 = g1->atom[1];  /* == g2->atom[0] */
              a3 = g1->atom[2];  /* == g2->atom[1] */
              a4 = g2->atom[2];
              found = ( a1->atomid < a4->atomid );
              if ( a1 != a4 ) ++count2;
            }
          } else if ( g2->atom[2] == atom ) {  /*  XBX XXB  */
            if ( g2->atom[1] == g1->atom[0] ) {  /*  CBA DCB  */
              a1 = g1->atom[2];
              a2 = g1->atom[1];  /* == g2->atom[2] */
              a3 = g1->atom[0];  /* == g2->atom[1] */
              a4 = g2->atom[0];
              found = ( a1->atomid < a4->atomid );
              if ( a1 != a4 ) ++count2;
            } else if ( g2->atom[1] == g1->atom[2] ) {  /*  ABC DCB  */
              a1 = g1->atom[0];
              a2 = g1->atom[1];  /* == g2->atom[2] */
              a3 = g1->atom[2];  /* == g2->atom[1] */
              a4 = g2->atom[0];
              found = ( a1->atomid < a4->atomid );
              if ( a1 != a4 ) ++count2;
            }
          } else return -6;
          if ( ! found ) continue;
          ++count1;
          tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
          if ( ! tuple ) return -10;
          tuple->next[0] = a1->dihedrals;
          tuple->atom[0] = a1;
          tuple->next[1] = a2->dihedrals;
          tuple->atom[1] = a2;
          tuple->next[2] = a3->dihedrals;
          tuple->atom[2] = a3;
          tuple->next[3] = a4->dihedrals;
          tuple->atom[3] = a4;
          tuple->del = 0;
          a1->dihedrals = tuple;
          a2->dihedrals = tuple;
          a3->dihedrals = tuple;
          a4->dihedrals = tuple;
        }
      }
    }
  }
  }
  if ( count2 != 2 * count1 ) return -15;  /* missing dihedrals */
  
#else

    for ( i = 0; i < res->atomSize; i++ ) {
      a1 = res->atomArray[i];
      if (a1->del) continue; 
      for ( bond1 = a1->bonds; bond1; bond1 = topo_mol_bond_next(bond1,a1) ) {
        
        if (bond1->del) continue;
        
        a2 = bond1->atom[0];
        if (a2->atomid == a1->atomid) a2 = bond1->atom[1];

        if (!a2) continue;
        
        for ( bond2 = a2->bonds; bond2;
      		bond2 = topo_mol_bond_next(bond2,a2) ) {
          
          if (bond2->del) continue;
          
          a3 = bond2->atom[0];
          if (a3->atomid == a2->atomid) a3 = bond2->atom[1];

          if (!a3 || a3->atomid == a1->atomid) continue;

          if ( is_hydrogen(a2) && ( ! topo_mol_bond_next(bond2,atom) ) &&
               ( ( is_hydrogen(a1) && is_oxygen(a3) ) ||
                 ( is_hydrogen(a3) && is_oxygen(a1) ) ) )
            continue;  /* extra H-H bond on water */
            
          for ( bond3 = a3->bonds; bond3;
        		bond3 = topo_mol_bond_next(bond3,a3) ) {
            
            if (bond3->del) continue;
            
            a4 = bond3->atom[0];
            if (a4->atomid == a3->atomid) a4 = bond3->atom[1];
            if (!a4 || a4->atomid == a2->atomid || a4->atomid <= a1->atomid) continue;
            
            tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
            if ( ! tuple ) return -10;
            tuple->next = a1->dihedrals;
            tuple->atom[0] = a2;
            tuple->atom[1] = a3;
            tuple->atom[2] = a4;
            tuple->del = 0;
            a1->dihedrals = tuple;
          }
        }
      }
    }
  }
 }
#endif
  
  return 0;
}

int topo_mol_patch(topo_mol *mol, const topo_mol_ident_t *targets,
                        int ntargets, const char *rname, int prepend,
			int warn_angles, int warn_dihedrals, int deflt) {
  int idef;
  topo_defs_residue_t *resdef;
  topo_defs_atom_t *atomdef;
  topo_defs_bond_t *bonddef;
  topo_defs_angle_t *angldef;
  topo_defs_dihedral_t *dihedef;
  topo_defs_improper_t *imprdef;
  topo_defs_cmap_t *cmapdef;
  topo_defs_conformation_t *confdef;
  topo_mol_residue_t *res=NULL, *oldres=NULL;
  
#if !defined(NEWPSFGEN)  
  topo_mol_atom_t *oldatoms=NULL;
#else 
  topo_defs_anisotropy_t *anisodef;
#endif
  
  char errmsg[128];

  if ( ! mol ) return -1;
  if ( mol->buildseg ) return -2;
  if ( ! mol->defs ) return -3;

  idef = hasharray_index(mol->defs->residue_hash,rname);
  if ( idef == HASHARRAY_FAIL ) {
    sprintf(errmsg,"unknown patch type %s",rname);
    topo_mol_log_error(mol,errmsg);
    return -4;
  }
  resdef = &(mol->defs->residue_array[idef]);
  if ( ! resdef->patch ) {
    sprintf(errmsg,"unknown patch type %s",rname);
    topo_mol_log_error(mol,errmsg);
    return -5;
  }

  oldres = 0;
  for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
    if ( atomdef->res < 0 || atomdef->res >= ntargets ) return -6;
    res = topo_mol_get_res(mol,&targets[atomdef->res],atomdef->rel);
    if ( ! res ) return -7;
    
      
#if !defined(NEWPSFGEN)  
  
    if ( atomdef->del ) {
      topo_mol_del_atom(res,atomdef->name);
      oldres = 0;
      continue;
    }
    if ( res != oldres ) {
      oldres = res;
      oldatoms = res->atoms;
    }
    if ( atomdef->type[0] == '\0' ) {
      topo_mol_find_atom(&(res->atoms),oldatoms,atomdef->name);
    } else if ( topo_mol_add_atom(mol,&(res->atoms),oldatoms,atomdef) ) {
      sprintf(errmsg,"add atom failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
      return -8;
    }

#else
        
    if ( atomdef->del ) {
      topo_mol_del_atom(res,atomdef->name);
      oldres = 0;
      continue;
    }
    if ( res != oldres ) {
      oldres = res;
    }
  
    if (atomdef->type[0] != '\0' && topo_mol_add_atom(mol,resdef,atomdef,res) ) {
      sprintf(errmsg,"add atom failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
      return -8;
    }
    if (res->pres[0] == '\0' || strcmp(res->pres,resdef->name)) {
      strcpy(res->pres,resdef->name);
    }
    

#endif

  }
  
#if defined(NEWPSFGEN) 
    
  if (resdef->lonepairs == 1) {
    
    for ( atomdef = resdef->atoms; atomdef; atomdef = atomdef->next ) {
      if (atomdef->islonepair && !atomdef->del) {
        if ( topo_mol_update_lonepair(mol,atomdef,res,targets,1) ) { 
          sprintf(errmsg,"failed to update lonepair in residue %s:%s",res->name,res->resid);
          topo_mol_log_error(mol,errmsg);
          return -11;
        }
      }
    }
  }
  
  /* define the anisotropy in case of drude ff*/
  for (anisodef = resdef->aniso; anisodef; anisodef = anisodef->next) {
    /* since the anisotropy is defined at the residue level, in the case of patches
     * affecting two residues (e.g. DISU) the anisotropy is defined in the first residue
     * used to define the anisotropy entry in the patch
    */
    res = topo_mol_get_res(mol,&targets[anisodef->res[0]],anisodef->rel[0]); 
    topo_mol_add_anisotropy(mol,anisodef,res, targets);
  }
#endif

  for ( bonddef = resdef->bonds; bonddef; bonddef = bonddef->next ) {
    if ( bonddef->del ) topo_mol_del_bond(mol,targets,ntargets,bonddef);
    else if ( topo_mol_add_bond(mol,targets,ntargets,bonddef) ) {
      sprintf(errmsg,"Warning: add bond failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  if ( warn_angles && resdef->angles ) {
    sprintf(errmsg,"Warning: explicit angles in patch %s will be deleted during autogeneration",rname);
    topo_mol_log_error(mol,errmsg);
  }
  for ( angldef = resdef->angles; angldef; angldef = angldef->next ) {
    if ( angldef->del ) topo_mol_del_angle(mol,targets,ntargets,angldef);
    else if ( topo_mol_add_angle(mol,targets,ntargets,angldef) ) {
      sprintf(errmsg,"Warning: add angle failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  if ( warn_dihedrals && resdef->dihedrals ) {
    sprintf(errmsg,"Warning: explicit dihedrals in patch %s will be deleted during autogeneration",rname);
    topo_mol_log_error(mol,errmsg);
  }
  for ( dihedef = resdef->dihedrals; dihedef; dihedef = dihedef->next ) {
    if ( dihedef->del ) topo_mol_del_dihedral(mol,targets,ntargets,dihedef);
    else if ( topo_mol_add_dihedral(mol,targets,ntargets,dihedef) ) {
      sprintf(errmsg,"Warning: add dihedral failed in patch %s",rname);
        topo_mol_log_error(mol,errmsg);
      }
  }
  for ( imprdef = resdef->impropers; imprdef; imprdef = imprdef->next ) {
    if ( imprdef->del ) topo_mol_del_improper(mol,targets,ntargets,imprdef);
    else if ( topo_mol_add_improper(mol,targets,ntargets,imprdef) ) {
      sprintf(errmsg,"Warning: add improper failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  for ( cmapdef = resdef->cmaps; cmapdef; cmapdef = cmapdef->next ) {
    if ( cmapdef->del ) topo_mol_del_cmap(mol,targets,ntargets,cmapdef);
    else if ( topo_mol_add_cmap(mol,targets,ntargets,cmapdef) ) {
      sprintf(errmsg,"Warning: add cross-term failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }
  for ( confdef = resdef->conformations; confdef; confdef = confdef->next ) {
    if ( confdef->del ) topo_mol_del_conformation(mol,targets,ntargets,confdef);
    else if ( topo_mol_add_conformation(mol,targets,ntargets,confdef) ) {
      sprintf(errmsg,"Warning: add conformation failed in patch %s",rname);
      topo_mol_log_error(mol,errmsg);
    }
  }

  if (strncasecmp(rname,"NONE",4)) {
    int ret;
    ret = topo_mol_add_patch(mol,rname,deflt);
    if (ret<0) {
      sprintf(errmsg,"Warning: Listing patch %s failed!",rname);
      topo_mol_log_error(mol,errmsg);
   }
    for ( idef=0; idef<ntargets; idef++ ) {
      // printf("%s:%s ", targets[idef].segid,targets[idef].resid);
      topo_mol_add_patchres(mol,&targets[idef]);
    }
    printf("\n");
  }
  return 0;
}

#if !defined(NEWPSFGEN)

int topo_mol_multiply_atoms(topo_mol *mol, const topo_mol_ident_t *targets,
						int ntargets, int ncopies) {
  int ipass, natoms, iatom, icopy;
  const topo_mol_ident_t *target;
  int itarget;
  topo_mol_atom_t *atom, **atoms;
  topo_mol_residue_t *res;
  topo_mol_segment_t *seg;
  int nres, ires;
  

  if (!mol) return -1;

  /* Quiet compiler warnings */
  natoms = 0; 
  atoms = NULL;

  /* two passes needed to find atoms */
  /* the fist passage is to just count the atoms to know the dimension of the 
   * array. Also the second passage updates the bonds, angles, and etc.
  */
  for (ipass=0; ipass<2; ++ipass) {
    if ( ipass ) atoms = memarena_alloc(mol->arena,
				natoms*sizeof(topo_mol_atom_t*));
    natoms = 0;
    /* walk all targets */
    for (itarget=0; itarget<ntargets; ++itarget) {
      target = targets + itarget;
      
      if (!target->resid) { /* whole segment */
        seg = topo_mol_get_seg(mol,target);
        if ( ! seg ) return -2;
        nres = hasharray_count(seg->residue_hash);
        for ( ires=0; ires<nres; ++ires ) {
          res = &(seg->residue_array[ires]);

          for ( atom = res->atoms; atom; atom = atom->next ) {
            if ( ipass ) atoms[natoms] = atom;
            ++natoms;
          }
            
        }
        continue;
      }

      if (!target->aname) { /* whole residue */
        res = topo_mol_get_res(mol,target,0);
        if ( ! res ) return -3;

        for ( atom = res->atoms; atom; atom = atom->next ) {
          if ( ipass ) atoms[natoms] = atom;
          ++natoms;
        }
                
        continue;
      }

      /* one atom */
      atom = topo_mol_get_atom(mol,target,0);
      if ( ! atom ) return -4;
      if ( ipass ) atoms[natoms] = atom;
      ++natoms;
    }
  }
  
  /* make one copy on each pass through loop */
  for (icopy=1; icopy<ncopies; ++icopy) {

  /* copy the actual atoms */
  for (iatom=0; iatom<natoms; ++iatom) {
    topo_mol_atom_t *newatom;     

    atom = atoms[iatom];
    
    if ( atom->copy ) {
      topo_mol_log_error(mol,"an atom occurs twice in the selection");
      return -20;
    }
    newatom = memarena_alloc(mol->arena,sizeof(topo_mol_atom_t));
    if ( ! newatom ) return -5;
    memcpy(newatom,atom,sizeof(topo_mol_atom_t));   

    atom->next = newatom;

    atom->copy = newatom;
    newatom->bonds = 0;
    newatom->angles = 0;
    newatom->dihedrals = 0;
    newatom->impropers = 0;
    newatom->cmaps = 0;
    newatom->exclusions = 0;
    newatom->conformations = 0;



  }

  /* copy associated bonds, etc. */
  for (iatom=0; iatom<natoms; ++iatom) {
    topo_mol_atom_t *a1, *a2, *a3, *a4;
    topo_mol_bond_t *bondtmp;
    topo_mol_angle_t *angletmp;
    topo_mol_dihedral_t *dihetmp;
    topo_mol_improper_t *imprtmp;
    topo_mol_cmap_t *cmaptmp;
    topo_mol_exclusion_t *excltmp;
    topo_mol_conformation_t *conftmp;
    atom = atoms[iatom];
    

    for ( bondtmp = atom->bonds; bondtmp;
		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
      topo_mol_bond_t *tuple;
      if ( bondtmp->del ) continue;
      if ( bondtmp->atom[0] == atom || ( ! bondtmp->atom[0]->copy ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
      if ( ! tuple ) return -6;
      a1 = bondtmp->atom[0]->copy; if ( ! a1 ) a1 = bondtmp->atom[0];
      a2 = bondtmp->atom[1]->copy; if ( ! a2 ) a2 = bondtmp->atom[1];
      tuple->next[0] = a1->bonds;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->bonds;
      tuple->atom[1] = a2;
      tuple->del = 0;
      a1->bonds = tuple;
      a2->bonds = tuple;
    }

    for ( angletmp = atom->angles; angletmp;
		angletmp = topo_mol_angle_next(angletmp,atom) ) {
      topo_mol_angle_t *tuple;
      if ( angletmp->del ) continue;
      if ( angletmp->atom[0] == atom || ( ! angletmp->atom[0]->copy
      && ( angletmp->atom[1] == atom || ( ! angletmp->atom[1]->copy ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->angle_arena,sizeof(topo_mol_angle_t));
      if ( ! tuple ) return -7;
      a1 = angletmp->atom[0]->copy; if ( ! a1 ) a1 = angletmp->atom[0];
      a2 = angletmp->atom[1]->copy; if ( ! a2 ) a2 = angletmp->atom[1];
      a3 = angletmp->atom[2]->copy; if ( ! a3 ) a3 = angletmp->atom[2];
      tuple->next[0] = a1->angles;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->angles;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->angles;
      tuple->atom[2] = a3;
      tuple->del = 0;
      a1->angles = tuple;
      a2->angles = tuple;
      a3->angles = tuple;
    }

    for ( dihetmp = atom->dihedrals; dihetmp;
    dihetmp = topo_mol_dihedral_next(dihetmp,atom) ) {
      topo_mol_dihedral_t *tuple;
      if ( dihetmp->del ) continue;
      if ( dihetmp->atom[0] == atom || ( ! dihetmp->atom[0]->copy
      && ( dihetmp->atom[1] == atom || ( ! dihetmp->atom[1]->copy
      && ( dihetmp->atom[2] == atom || ( ! dihetmp->atom[2]->copy ) ) ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
      if ( ! tuple ) return -8;
      a1 = dihetmp->atom[0]->copy; if ( ! a1 ) a1 = dihetmp->atom[0];
      a2 = dihetmp->atom[1]->copy; if ( ! a2 ) a2 = dihetmp->atom[1];
      a3 = dihetmp->atom[2]->copy; if ( ! a3 ) a3 = dihetmp->atom[2];
      a4 = dihetmp->atom[3]->copy; if ( ! a4 ) a4 = dihetmp->atom[3];
      tuple->next[0] = a1->dihedrals;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->dihedrals;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->dihedrals;
      tuple->atom[2] = a3;
      tuple->next[3] = a4->dihedrals;
      tuple->atom[3] = a4;
      tuple->del = 0;
      a1->dihedrals = tuple;
      a2->dihedrals = tuple;
      a3->dihedrals = tuple;
      a4->dihedrals = tuple;

    }
    
    for ( imprtmp = atom->impropers; imprtmp;
		imprtmp = topo_mol_improper_next(imprtmp,atom) ) {
      topo_mol_improper_t *tuple;
      if ( imprtmp->del ) continue;
      if ( imprtmp->atom[0] == atom || ( ! imprtmp->atom[0]->copy
      && ( imprtmp->atom[1] == atom || ( ! imprtmp->atom[1]->copy
      && ( imprtmp->atom[2] == atom || ( ! imprtmp->atom[2]->copy ) ) ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
      if ( ! tuple ) return -9;
      a1 = imprtmp->atom[0]->copy; if ( ! a1 ) a1 = imprtmp->atom[0];
      a2 = imprtmp->atom[1]->copy; if ( ! a2 ) a2 = imprtmp->atom[1];
      a3 = imprtmp->atom[2]->copy; if ( ! a3 ) a3 = imprtmp->atom[2];
      a4 = imprtmp->atom[3]->copy; if ( ! a4 ) a4 = imprtmp->atom[3];
      tuple->next[0] = a1->impropers;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->impropers;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->impropers;
      tuple->atom[2] = a3;
      tuple->next[3] = a4->impropers;
      tuple->atom[3] = a4;
      tuple->del = 0;
      a1->impropers = tuple;
      a2->impropers = tuple;
      a3->impropers = tuple;
      a4->impropers = tuple;
    }

    for ( cmaptmp = atom->cmaps; cmaptmp;
		cmaptmp = topo_mol_cmap_next(cmaptmp,atom) ) {
      topo_mol_atom_t *al[8];
      topo_mol_cmap_t *tuple;
      int ia, skip;
      if ( cmaptmp->del ) continue;
      skip = 0;
      for ( ia = 0; ia < 8; ++ia ) {
        if ( cmaptmp->atom[ia] == atom ) { skip = 0; break; }
        if ( cmaptmp->atom[ia]->copy ) { skip = 1; break; }
      }
      if ( skip ) continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_cmap_t));
      if ( ! tuple ) return -9;
      for ( ia = 0; ia < 8; ++ia ) {
        topo_mol_atom_t *ai;
        ai = cmaptmp->atom[ia]->copy;
        if ( ! ai ) ai = cmaptmp->atom[ia];
        al[ia] = ai;
        tuple->next[ia] = ai->cmaps;
        tuple->atom[ia] = ai;
      }
      for ( ia = 0; ia < 8; ++ia ) {
        /* This must be in a separate loop because atoms may be repeated. */
        al[ia]->cmaps = tuple;
      }
      tuple->del = 0;
    }
    for ( excltmp = atom->exclusions; excltmp;
		excltmp = topo_mol_exclusion_next(excltmp,atom) ) {
      topo_mol_exclusion_t *tuple;
      if ( excltmp->del ) continue;
      if ( excltmp->atom[0] == atom || ( ! excltmp->atom[0]->copy ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_exclusion_t));
      if ( ! tuple ) return -6;
      a1 = excltmp->atom[0]->copy; if ( ! a1 ) a1 = excltmp->atom[0];
      a2 = excltmp->atom[1]->copy; if ( ! a2 ) a2 = excltmp->atom[1];
      tuple->next[0] = a1->exclusions;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->exclusions;
      tuple->atom[1] = a2;
      tuple->del = 0;
      a1->exclusions = tuple;
      a2->exclusions = tuple;
    }
    for ( conftmp = atom->conformations; conftmp;
		conftmp = topo_mol_conformation_next(conftmp,atom) ) {
      topo_mol_conformation_t *tuple;
      if ( conftmp->del ) continue;
      if ( conftmp->atom[0] == atom || ( ! conftmp->atom[0]->copy
      && ( conftmp->atom[1] == atom || ( ! conftmp->atom[1]->copy
      && ( conftmp->atom[2] == atom || ( ! conftmp->atom[2]->copy ) ) ) ) ) ) ;
      else continue;
      tuple = memarena_alloc(mol->arena,sizeof(topo_mol_conformation_t));
      if ( ! tuple ) return -10;
      a1 = conftmp->atom[0]->copy; if ( ! a1 ) a1 = conftmp->atom[0];
      a2 = conftmp->atom[1]->copy; if ( ! a2 ) a2 = conftmp->atom[1];
      a3 = conftmp->atom[2]->copy; if ( ! a3 ) a3 = conftmp->atom[2];
      a4 = conftmp->atom[3]->copy; if ( ! a4 ) a4 = conftmp->atom[3];
      tuple->next[0] = a1->conformations;
      tuple->atom[0] = a1;
      tuple->next[1] = a2->conformations;
      tuple->atom[1] = a2;
      tuple->next[2] = a3->conformations;
      tuple->atom[2] = a3;
      tuple->next[3] = a4->conformations;
      tuple->atom[3] = a4;
      tuple->del = 0;
      tuple->improper = conftmp->improper;
      tuple->dist12 = conftmp->dist12;
      tuple->angle123 = conftmp->angle123;
      tuple->dihedral = conftmp->dihedral;
      tuple->angle234 = conftmp->angle234;
      tuple->dist34 = conftmp->dist34;
      a1->conformations = tuple;
      a2->conformations = tuple;
      a3->conformations = tuple;
      a4->conformations = tuple;
    }
  }

  /* clean up copy pointers */
  for (iatom=0; iatom<natoms; ++iatom) {
    atom = atoms[iatom];
    if ( atom->partition == 0 ) atom->partition = 1;
    atom->copy->partition = atom->partition + 1;
    atoms[iatom] = atom->copy;
    atom->copy = 0;
  }

  } /* icopy */

  return 0;  /* success */
}

#else 
int topo_mol_multiply_atoms_exec(topo_mol *mol, topo_mol_residue_t *res, 
          int ncopies, topo_mol_atom_t *atomToCopy);

int topo_mol_multiply_atom_update_bonds(topo_mol *mol, topo_mol_residue_t *res, 
          int ncopies, topo_mol_atom_t *atomToCopy);

int topo_mol_multiply_atoms(topo_mol *mol, const topo_mol_ident_t *targets,
            int ntargets, int ncopies) {
  int ipass;
  const topo_mol_ident_t *target;
  int itarget;
  topo_mol_atom_t *atom;
  topo_mol_residue_t *res;
  topo_mol_segment_t *seg, *segprev;
  int nres, ires, i, error;

  topo_mol_atom_t *a1, *a2, *a3, *aux;
  topo_mol_dihedral_t *dihetmp;
  topo_mol_improper_t *imprtmp;
  if (!mol) return -1;

  /* two passes needed to find atoms */
  /* the fist passage is to just count the atoms to know the dimension of the 
   * array. Also the second passage updates the bonds, angles, and etc.
  */


  for (ipass=0; ipass<2; ++ipass) {
    /* walk all targets */
   
    for (itarget=0; itarget<ntargets; ++itarget) {
      target = targets + itarget;
      
      if (!target->resid) { /* whole segment */
        seg = topo_mol_get_seg(mol,target);
        if ( ! seg ) return -2;
        nres = hasharray_count(seg->residue_hash);
        for ( ires=0; ires<nres; ++ires ) {
          res = topo_mol_get_res(mol,target,0);
          if ( ! res ) return -3;

          res = &(seg->residue_array[ires]);
          if (!ipass) {
            error = topo_mol_multiply_atoms_exec(mol, res, ncopies,0);
          
          } else {
            error = topo_mol_multiply_atom_update_bonds(mol, res, ncopies,0 );
          }
          
          if (error) return error ;

        }
        continue;
      }
      if (!target->aname) {/* whole residue */

        res = topo_mol_get_res(mol,target,0);
        if ( ! res ) return -3;

        if (!ipass) {
          error = topo_mol_multiply_atoms_exec(mol, res, ncopies,0);
        } else {
          error = topo_mol_multiply_atom_update_bonds(mol, res, ncopies, 0);
        }
        continue;
      }

      /* one atom */
      if (!ipass) {
        res = topo_mol_get_res(mol,target,0);
        if ( ! res ) return -3;

        atom = topo_mol_get_atom(mol,target,0);
        if ( ! atom ) return -4;
        error = topo_mol_multiply_atoms_exec(mol, res, ncopies, atom);
      } else {
        res = topo_mol_get_res(mol,target,0);
        atom = topo_mol_get_atom(mol,target,0);
        error = topo_mol_multiply_atom_update_bonds(mol, res, ncopies, atom);
      }
      
      
    }
    
  }

  segprev = 0;
  for (itarget=0; itarget<ntargets; ++itarget) {

    target = targets + itarget;
    seg = topo_mol_get_seg(mol,target);
    if ( ! seg ) return -2;

    if (seg == segprev) continue;

    nres = hasharray_count(seg->residue_hash);
   
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        
        if (atom->partition == 0) {
          aux = atom;
        } else {
          aux = atom->copy;
        }

        for ( dihetmp = atom->dihedrals; dihetmp;
              dihetmp = dihetmp->next ) {

          topo_mol_dihedral_t *tuple;
          if ( dihetmp->del) continue;

          if (  atom->copy ||
                dihetmp->atom[0]->copy || 
                dihetmp->atom[1]->copy || 
                dihetmp->atom[2]->copy ) {
            tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
            if ( ! tuple ) return -8;
            a1 = dihetmp->atom[0]->copy; if ( ! a1 ) a1 = dihetmp->atom[0];
            a2 = dihetmp->atom[1]->copy; if ( ! a2 ) a2 = dihetmp->atom[1];
            a3 = dihetmp->atom[2]->copy; if ( ! a3 ) a3 = dihetmp->atom[2];
           
            tuple->atom[0] = a1;
            tuple->atom[1] = a2;
            tuple->atom[2] = a3;
            tuple->del = 0;

            /* Add the new dihedral/improper with the copies after the current position
             * so in case the copy also points to another copy, the loop
             * will be able to identify the copy
            */ 
            if (atom->partition == 0) {
              if (dihetmp->next) {
                tuple->next = dihetmp->next;
              } else {
                tuple->next = 0;
              }
              dihetmp->next = tuple;
            } else {
              tuple->next = aux->dihedrals;
              aux->dihedrals = tuple;
            }
          }
        }

        for ( imprtmp = atom->impropers; imprtmp;
            imprtmp = imprtmp->next ) {
          topo_mol_improper_t *tuple;
          if ( imprtmp->del ) continue;

          if (  atom->copy ||
                imprtmp->atom[0]->copy || 
                imprtmp->atom[1]->copy || 
                imprtmp->atom[2]->copy) {
            tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));
            if ( ! tuple ) return -9;
            a1 = imprtmp->atom[0]->copy; if ( ! a1 ) a1 = imprtmp->atom[0];
            a2 = imprtmp->atom[1]->copy; if ( ! a2 ) a2 = imprtmp->atom[1];
            a3 = imprtmp->atom[2]->copy; if ( ! a3 ) a3 = imprtmp->atom[2];
            tuple->atom[0] = a1;
            tuple->atom[1] = a2;
            tuple->atom[2] = a3;
            tuple->del = 0;

            if (atom->partition == 0) {
              if (imprtmp->next) {
                tuple->next = imprtmp->next;
              } else {
                tuple->next = 0;
              }
              imprtmp->next = tuple;
            } else {
              tuple->next = aux->impropers;
              aux->impropers = tuple;
            }
          }
        }
      }
    }
    segprev = seg;
  }
  return 0;
}

int topo_mol_multiply_atoms_exec(topo_mol *mol, topo_mol_residue_t *res, 
          int ncopies, topo_mol_atom_t *atomToCopy) {

  int i, h, j =0, newlimit=0;
  topo_mol_atom_t **tmp;

  // new limit is ncopies times the size of the current array. 
  newlimit = res->atomSize * ncopies;
  // In case of multiplying single atoms
  if (atomToCopy) {
    newlimit = res->atomSize + ncopies -1;
  }
  

  // create a tmp array to hold the multiplied residue
  tmp = (topo_mol_atom_t **) calloc(newlimit, sizeof(topo_mol_atom_t*));
  j = 0;
  // the array has the atoms multiplied next to each other
  for  (i = 0; i < newlimit; i++ ) {
    tmp[i] = res->atomArray[j];

    if (atomToCopy && strcmp(tmp[i]->name, atomToCopy->name) ) {
      j++;
      continue;
    }

    for (h = i + 1 ; h < i + ncopies; h++) {
      tmp[h] = memarena_alloc(mol->arena,sizeof(topo_mol_atom_t));
      if ( ! tmp[h] ) return -5;
      memcpy(tmp[h],tmp[i],sizeof(topo_mol_atom_t));
      tmp[h]->bonds = 0;
      tmp[h]->angles = 0;
      tmp[h]->dihedrals = 0;
      tmp[h]->impropers = 0;
      tmp[h]->cmaps = 0;
      tmp[h]->exclusions = 0;
      tmp[h]->conformations = 0;
      tmp[h - 1]->copy = tmp[h];
      tmp[h]->copy = 0;
    }

    if (!atomToCopy) {
      /* increment the number of copies in the cases of multiplying residues
       * The cases of multiplying only atom, the copy should increment +
      */
      i += ncopies -1;
    } else if (atomToCopy && !strcmp(tmp[i]->name, atomToCopy->name)) {
      i += ncopies -1;
    }

    j++;
  }
  // Assign new array size
  res->atomSize = newlimit;
  res->reordered = 1;

  res->atomArray = (topo_mol_atom_t **) realloc(res->atomArray, (newlimit)*sizeof(topo_mol_atom_t*));

  // Copy the tmp array to the res->atomArray original array
  memcpy(&(res->atomArray), &tmp,sizeof(newlimit));


  // TODO: Copy the information about the lonepairs and drude particles 
/* #if defined(NEWPSFGEN)   
     newatom->lonepair = 0;
     newatom->islonepair = 0;
     newatom->del = 0;
     newatom->dxyz = 0;
 #endif
*/
  return 0;
}


int topo_mol_multiply_atom_update_bonds(topo_mol *mol, topo_mol_residue_t *res, 
          int ncopies, topo_mol_atom_t *atomToCopy) {
  topo_mol_atom_t *a1, *a2, *a3, *a4;
  topo_mol_bond_t *bondtmp;
  topo_mol_angle_t *angletmp;
  topo_mol_cmap_t *cmaptmp;
  topo_mol_exclusion_t *excltmp;
  topo_mol_conformation_t *conftmp;
  int i, h =0;
  topo_mol_atom_t *atom;
  /* Copy the rest of the atoms' filed, like bonds, dihedrals and such.
   * The way the following works is that the first atom points to the copy,
   * and the copy of the first atom, points to the next copy and so on...
  */
  for  (i = 0; i < res->atomSize; i++) {
    
    if (atomToCopy) {
      if (!strcmp(res->atomArray[i]->name, atomToCopy->name) && \
        res->atomArray[i]->partition == 0) res->atomArray[i]->partition = 1;
    } else {
      if (res->atomArray[i]->partition == 0 ) res->atomArray[i]->partition = 1;
    }  

    for (h = i ; h < i + ncopies; h++) {
      if (atomToCopy && strcmp(res->atomArray[i]->name, atomToCopy->name)) {
        break;
      }
      atom = res->atomArray[h];
      if (!atom->copy) continue; 
      atom->copy->partition = atom->partition + 1;

      for ( bondtmp = atom->bonds; bondtmp;
            bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
        topo_mol_bond_t *tuple;
        if ( bondtmp->del ) continue;
        if ( bondtmp->atom[0] == atom || ( ! bondtmp->atom[0]->copy ) ) ;
        else continue;
        tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
        if ( ! tuple ) return -6; /* XXX what is -6? */
        a1 = bondtmp->atom[0]->copy; if ( ! a1 ) {
          a1 = bondtmp->atom[0];
        }
        a2 = bondtmp->atom[1]->copy; if ( ! a2 ) {
          a2 = bondtmp->atom[1];
        }
        tuple->next[0] = a1->bonds;
        tuple->atom[0] = a1;
        tuple->next[1] = a2->bonds;
        tuple->atom[1] = a2;
        tuple->del = 0;
        a1->bonds = tuple;
        a2->bonds = tuple;

        /* check if the non-copy part of the bond (if any) has any atom that as copied
         * declared in the dihedrals and impropers
        */

      }

      for ( angletmp = atom->angles; angletmp;
            angletmp = topo_mol_angle_next(angletmp,atom) ) {
        topo_mol_angle_t *tuple;
        if ( angletmp->del ) continue;
        if ( angletmp->atom[0] == atom || ( ! angletmp->atom[0]->copy
        && ( angletmp->atom[1] == atom || ( ! angletmp->atom[1]->copy ) ) ) ) ;
        else continue;
        tuple = memarena_alloc(mol->angle_arena,sizeof(topo_mol_angle_t));
        if ( ! tuple ) return -7; /* XXX what is -7? */
        a1 = angletmp->atom[0]->copy; if ( ! a1 ) a1 = angletmp->atom[0];
        a2 = angletmp->atom[1]->copy; if ( ! a2 ) a2 = angletmp->atom[1];
        a3 = angletmp->atom[2]->copy; if ( ! a3 ) a3 = angletmp->atom[2];
        tuple->next[0] = a1->angles;
        tuple->atom[0] = a1;
        tuple->next[1] = a2->angles;
        tuple->atom[1] = a2;
        tuple->next[2] = a3->angles;
        tuple->atom[2] = a3;
        tuple->del = 0;
        a1->angles = tuple;
        a2->angles = tuple;
        a3->angles = tuple;
      }



      for ( cmaptmp = atom->cmaps; cmaptmp;
          cmaptmp = topo_mol_cmap_next(cmaptmp,atom) ) {
        topo_mol_atom_t *al[8];
        topo_mol_cmap_t *tuple;
        int ia, skip;
        if ( cmaptmp->del ) continue;
        skip = 0;
        for ( ia = 0; ia < 8; ++ia ) {
          if ( cmaptmp->atom[ia] == atom ) { skip = 0; break; }
          if ( cmaptmp->atom[ia]->copy ) { skip = 1; break; }
        }
        if ( skip ) continue;
        tuple = memarena_alloc(mol->arena,sizeof(topo_mol_cmap_t));
        if ( ! tuple ) return -9; /* XXX what is -9? */
        for ( ia = 0; ia < 8; ++ia ) {
          topo_mol_atom_t *ai;
          ai = cmaptmp->atom[ia]->copy;
          if ( ! ai ) ai = cmaptmp->atom[ia];
          al[ia] = ai;
          tuple->next[ia] = ai->cmaps;
          tuple->atom[ia] = ai;
        }
        for ( ia = 0; ia < 8; ++ia ) {
          /* This must be in a separate loop because atoms may be repeated. */
          al[ia]->cmaps = tuple;
        }
        tuple->del = 0;
      }
      
      
      for ( excltmp = atom->exclusions; excltmp;
          excltmp = topo_mol_exclusion_next(excltmp,atom) ) {
        topo_mol_exclusion_t *tuple;
        if ( excltmp->del ) continue;
        if ( excltmp->atom[0] == atom || ( ! excltmp->atom[0]->copy ) ) ;
        else continue;
        tuple = memarena_alloc(mol->arena,sizeof(topo_mol_exclusion_t));
        if ( ! tuple ) return -6; /* XXX what is -6? */
        a1 = excltmp->atom[0]->copy; if ( ! a1 ) a1 = excltmp->atom[0];
        a2 = excltmp->atom[1]->copy; if ( ! a2 ) a2 = excltmp->atom[1];
        tuple->next[0] = a1->exclusions;
        tuple->atom[0] = a1;
        tuple->next[1] = a2->exclusions;
        tuple->atom[1] = a2;
        tuple->del = 0;
        a1->exclusions = tuple;
        a2->exclusions = tuple;
      }

      for ( conftmp = atom->conformations; conftmp;
            conftmp = topo_mol_conformation_next(conftmp,atom) ) {
        topo_mol_conformation_t *tuple;
        if ( conftmp->del ) continue;
        if ( conftmp->atom[0] == atom || ( ! conftmp->atom[0]->copy
        && ( conftmp->atom[1] == atom || ( ! conftmp->atom[1]->copy
        && ( conftmp->atom[2] == atom || ( ! conftmp->atom[2]->copy ) ) ) ) ) ) ;
        else continue;
        tuple = memarena_alloc(mol->arena,sizeof(topo_mol_conformation_t));
        if ( ! tuple ) return -10; /* XXX what is -10? */
        a1 = conftmp->atom[0]->copy; if ( ! a1 ) a1 = conftmp->atom[0];
        a2 = conftmp->atom[1]->copy; if ( ! a2 ) a2 = conftmp->atom[1];
        a3 = conftmp->atom[2]->copy; if ( ! a3 ) a3 = conftmp->atom[2];
        a4 = conftmp->atom[3]->copy; if ( ! a4 ) a4 = conftmp->atom[3];
        tuple->next[0] = a1->conformations;
        tuple->atom[0] = a1;
        tuple->next[1] = a2->conformations;
        tuple->atom[1] = a2;
        tuple->next[2] = a3->conformations;
        tuple->atom[2] = a3;
        tuple->next[3] = a4->conformations;
        tuple->atom[3] = a4;
        tuple->del = 0;
        tuple->improper = conftmp->improper;
        tuple->dist12 = conftmp->dist12;
        tuple->angle123 = conftmp->angle123;
        tuple->dihedral = conftmp->dihedral;
        tuple->angle234 = conftmp->angle234;
        tuple->dist34 = conftmp->dist34;
        a1->conformations = tuple;
        a2->conformations = tuple;
        a3->conformations = tuple;
        a4->conformations = tuple;
      }
    }
    if (!atomToCopy) {
      /* increment the number of copies in the cases of multiplying residues
       * The cases of multiplying only atom, the copy should increment +
      */
      i += ncopies -1;
    } else if (atomToCopy && !strcmp(res->atomArray[i]->name, atomToCopy->name)) {
      i += ncopies -1;
    }
  }

  /* XXX there was no explicit return value here!!!!! */
  return  0;
}
#endif

  
/* API function */
void topo_mol_delete_atom(topo_mol *mol, const topo_mol_ident_t *target) {
  topo_mol_residue_t *res;
  topo_mol_segment_t *seg;
  int ires, iseg;
#if defined(NEWPSFGEN)
  int i;
#endif

  if (!mol) return;

  iseg = hasharray_index(mol->segment_hash,target->segid);
  if ( iseg == HASHARRAY_FAIL ) {
    char errmsg[50];
    sprintf(errmsg,"no segment %s",target->segid);
    topo_mol_log_error(mol,errmsg);
    return;
  }
  seg = mol->segment_array[iseg];
  
  if (!target->resid) {
    /* Delete this segment */
    int nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
  #if !defined(NEWPSFGEN)  
      topo_mol_atom_t *atom;
  #endif
    
      res = &(seg->residue_array[ires]);
      
  #if !defined(NEWPSFGEN)  
      atom = res->atoms;
      while (atom) {
        topo_mol_destroy_atom(atom);
        atom = atom->next;
      }
      res->atoms = 0;
  #else
      for ( i = 0; i < res->atomSize; i++ ) {
        topo_mol_destroy_atom(res->atomArray[i]);
        res->atomArray[i]->del = 1;
      }
      res->atomSize=0;
  #endif  
    }
  
    hasharray_destroy(seg->residue_hash);
    mol->segment_array[iseg] = 0;
    if (hasharray_delete(mol->segment_hash, target->segid) < 0) {
      topo_mol_log_error(mol, "Unable to delete segment");
    }
    return;
  }

  ires = hasharray_index(seg->residue_hash,target->resid);
  if ( ires == HASHARRAY_FAIL ) {
    char errmsg[50];
    sprintf(errmsg,"no residue %s of segment %s",
                                        target->resid,target->segid);
    topo_mol_log_error(mol,errmsg);
    return;
  }
  res = seg->residue_array+ires;  
  
  if (!target->aname) {
    /* Must destroy all atoms in residue, since there may be bonds between
     * this residue and other atoms 
     */
    
#if !defined(NEWPSFGEN)  

    topo_mol_atom_t *atom = res->atoms;
    while (atom) {
      topo_mol_destroy_atom(atom);
      atom = atom->next;
    }
    res->atoms = 0;
    
#else
  
    for ( i = 0; i < res->atomSize; i++ ) {
      topo_mol_destroy_atom(res->atomArray[i]);
      res->atomArray[i]->del = 1;
    }
    res->atomSize=0;
#endif

    hasharray_delete(seg->residue_hash, target->resid); 
    return;
  }
  
  /* Delete one atom */
#if !defined(NEWPSFGEN)  

  topo_mol_destroy_atom(topo_mol_unlink_atom(&(res->atoms),target->aname));

#else
  
  topo_mol_destroy_atom(topo_mol_unlink_atom(res,target->aname));
  
#endif

}


int topo_mol_set_name(topo_mol *mol, const topo_mol_ident_t *target,
                                     const char *name) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom = 0;
#if defined(NEWPSFGEN)  
  int i;
#endif
  
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;

#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else
  
  for ( i = 0; i < res->atomSize; i++ ) {
    if ( ! strcmp(target->aname,res->atomArray[i]->name) ) {
      atom = res->atomArray[i];
      break;
    }
  }

#endif

  if ( ! atom ) return -3;
  strcpy(atom->name,name);
  return 0;
}

int topo_mol_set_resname(topo_mol *mol, const topo_mol_ident_t *target,
                                        const char *rname) {
  topo_mol_residue_t *res;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;
  strcpy(res->name,rname);
  return 0;
}

int topo_mol_set_segid(topo_mol *mol, const topo_mol_ident_t *target,
                                      const char *segid) {
  int iseg, iseg2;
  topo_mol_segment_t *seg;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  seg = topo_mol_get_seg(mol,target);
  if ( ! seg ) return -3;
  iseg = hasharray_delete(mol->segment_hash, seg->segid);
  if ( iseg < 0) {
    topo_mol_log_error(mol, "Unable to delete segment");
    return -4;
  }
  strcpy(seg->segid,segid);
  iseg2 = hasharray_reinsert(mol->segment_hash, seg->segid, iseg);
  if ( iseg != iseg2 ) {
    topo_mol_log_error(mol, "Unable to insert segment");
    return -5;
  }
  return 0;
}

int topo_mol_set_element(topo_mol *mol, const topo_mol_ident_t *target,
                                        const char *element, int replace) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;
  
#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else

  atom = topo_mol_find_atom_by_name(res, target->aname);
        
#endif

  if ( ! atom ) return -3;

  if ( replace || ! strlen(atom->element) ) {
    strcpy(atom->element,element);
  }
  return 0;
}

int topo_mol_set_chain(topo_mol *mol, const topo_mol_ident_t *target,
                                        const char *chain, int replace) {
  topo_mol_residue_t *res;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;

  if ( replace || ! strlen(res->chain) ) {
    strcpy(res->chain,chain);
  }
  return 0;
}

#if defined(NEWPSFGEN)
/* Set the coordinates of the Colinear Lonepairs
 * The scale value is stores as angle variable in the lonepairs
 * data structure
*/

int topo_mol_set_COLINEARLP(topo_mol *mol, topo_mol_atom_t *atomLP) {
  double r12x,r12y,r12z,r12dist,a;
  topo_mol_atom_t *atom1, *atom2;
  topo_mol_lonepair_t *lonepair;
  char msg[128];
  
  atom1 = atomLP->lonepair->atoms[1];
  atom2 = atomLP->lonepair->atoms[2];
  lonepair = atomLP->lonepair;
  
  r12x = atom1->x - atom2->x;
  r12y = atom1->y - atom2->y;
  r12z = atom1->z - atom2->z;
  r12dist = sqrt(r12x*r12x + r12y*r12y + r12z*r12z); 
  
  if (r12dist < 1e-10 || 100. < r12dist) { // same low tolerance as used in CHARMM
    sprintf(msg,"Large/small distance between lonepair reference atoms: %s and %s",
  						atom1->name, atom2->name);
    topo_mol_log_error(mol,msg);
  }
  
  /* a is the normalized scaling factor used to place the lonepair.  
   * From NAMD HomePatch.C The lonepair is placed at a fixed (possibly negative) 
   * distance along the vector from k to j relative to a reference point. 
   * The reference point is computed by multiplying the vector from k to j 
   * by a (possibly negative) scale factor. For example, a scale value of 0 
   * (the default) sets the reference point as rj, while a value of -1 sets 
   * the reference point as rk. A scale value of -0.5 would use the center of 
   * the two hosts.
  */
  a = (lonepair->angle + lonepair->distance/r12dist);
  atomLP->x = atom1->x + a*r12x;
  atomLP->y = atom1->y + a*r12y;
  atomLP->z = atom1->z + a*r12z;
  
  return 0;
}


/* Extraction of the topo_mol_set_xyz to deal with the assignment
 * of xyz based on the Internal Coordinates, which is the same 
 * math applyied to assign the coordinates of the lonepairs, besides
 * the Colinear
*/
int topo_mol_set_IC(topo_mol_atom_t *atom0, topo_mol_atom_t *atom1,
  topo_mol_atom_t *atom2, topo_mol_atom_t *atom3, double distance,
  double angle, double dihedral, int bisector) {
    
  double r12x,r12y,r12z,r23x,r23y,r23z,ix,iy,iz,jx,jy,jz,kx,ky,kz;
  double tx,ty,tz,a,b,c;
  int gwild = 0;
  double vx, vy, vz;
  
  if (!bisector ) {
    vx = atom2->x;
    vy = atom2->y;
    vz = atom2->z;
  } else {
   /* Bisector only valid for guessing lonepairs as "For the BISEctor option, 
    * the dihedral is based on: I,J,(K+L)/2,L where I,J,K,L are the coordinate 
    * vectors of the specified atoms" (CHARMM doc)
    * Since the lonepairs are defined in oposite order, the K and L are atom1
    * and atom2
    */
    vx = (atom1->x + atom2->x)*0.50;
    vy = (atom1->y + atom2->y)*0.50;
    vz = (atom1->z + atom2->z)*0.50;
    
    /* make the distance positive again */
    distance *= -1;
  }
  
  r12x = vx - atom1->x;
  r12y = vy - atom1->y;
  r12z = vz - atom1->z;

  r23x = atom3->x - vx;
  r23y = atom3->y - vy;
  r23z = atom3->z - vz;
  
  /* Computing the x, y and z components of  ix, iy and iz of 
   * r23 needed to compute  
   * the final corss product 
   */
  
  a = sqrt(r23x*r23x + r23y*r23y + r23z*r23z); 
  if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
  ix = a * r23x;
  iy = a * r23y;
  iz = a * r23z;

  /* "Xyzzy (mnemonic)"way to calculate the cross product between r12
   * and r23
   * k = r12 x r23
   * t = (r12 x r23) x r23; t is used and reused, but at the end 
   * stores the final cross term
   */
   
  tx = r12y*r23z - r12z*r23y;
  ty = r12z*r23x - r12x*r23z;
  tz = r12x*r23y - r12y*r23x;
  
  /* Computing the  x, y and z components of 
   * r12 x r23
  */
  a = sqrt(tx*tx + ty*ty + tz*tz);
  if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
  kx = a * tx;
  ky = a * ty;
  kz = a * tz;

  /* Computing the cross product of the vector (r12 x r23) x r23 */
  tx = ky*iz - kz*iy;
  ty = kz*ix - kx*iz;
  tz = kx*iy - ky*ix;
  
  /* Computing the x, y and z components of the final 
   * (r12 x r23) x r23 needed to compute  
   * the final corss product 
  */
  a = sqrt(tx*tx + ty*ty + tz*tz);
  if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
  jx = a * tx;
  jy = a * ty;
  jz = a * tz;
  
  
  
  if (gwild) return gwild;
  

  a = -1.0 * distance * cos(angle);
  b = distance * sin(angle) * cos(dihedral);
  c = distance * sin(angle) * sin(dihedral);
  
  atom0->x = atom3->x + a * ix + b * jx + c * kx;
  atom0->y = atom3->y + a * iy + b * jy + c * ky;
  atom0->z = atom3->z + a * iz + b * jz + c * kz;

  
  return gwild;
}
#endif

int topo_mol_set_xyz(topo_mol *mol, const topo_mol_ident_t *target,
                                        double x, double y, double z) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;
  
#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else

  atom = topo_mol_find_atom_by_name(res, target->aname);
        
#endif

  if ( ! atom ) return -3;

  atom->x = x;
  atom->y = y;
  atom->z = z;
  atom->xyz_state = TOPO_MOL_XYZ_SET;
  return 0;
}

#if defined(NEWPSFGEN)
  
int topo_mol_set_drude_xyz(topo_mol *mol, const topo_mol_ident_t *target,
                                        double x, double y, double z) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  const char *atomname;
  if ( ! mol ) return -1;
  if ( ! target ) return -1;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -1;

  /* delete the D from the atom's name to get drude host name */
  atomname = target->aname +1;
  atom = topo_mol_find_atom_by_name(res, atomname);

  if ( ! atom ) return -1;

  atom->dxyz = (double *)malloc(3 * sizeof(double));
  atom->dxyz[0] = x;
  atom->dxyz[1] = y;
  atom->dxyz[2] = z;
  return 0;
}

#endif

int topo_mol_set_vel(topo_mol *mol, const topo_mol_ident_t *target,
                                        double vx, double vy, double vz) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;
  
#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else

  atom = topo_mol_find_atom_by_name(res, target->aname);
        
#endif

  if ( ! atom ) return -3;

  atom->vx = vx;
  atom->vy = vy;
  atom->vz = vz;
  return 0;
}

int topo_mol_set_mass(topo_mol *mol, const topo_mol_ident_t *target,
                      double mass) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;

#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else

  atom = topo_mol_find_atom_by_name(res, target->aname);
        
#endif

  if ( ! atom ) return -3;

  atom->mass = mass;
  return 0;
}

int topo_mol_set_charge(topo_mol *mol, const topo_mol_ident_t *target,
                        double charge) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;

#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else

  atom = topo_mol_find_atom_by_name(res, target->aname);
        
#endif

  if ( ! atom ) return -3;

  atom->charge = charge;
  return 0;
}

int topo_mol_set_bfactor(topo_mol *mol, const topo_mol_ident_t *target, 
                         double bfactor) {
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;
  if ( ! target ) return -2;
  res = topo_mol_get_res(mol,target,0);
  if ( ! res ) return -3;

#if !defined(NEWPSFGEN)  

  for ( atom = res->atoms; atom; atom = atom->next ) {
    if ( ! strcmp(target->aname,atom->name) ) break;
  }

#else

  atom = topo_mol_find_atom_by_name(res, target->aname);
        
#endif

  if ( ! atom ) return -3;

  atom->partition = bfactor;
  return 0;
}

/* XXX Unused */
int topo_mol_clear_xyz(topo_mol *mol, const topo_mol_ident_t *target) {
  topo_mol_atom_t *atom;
  if ( ! mol ) return -1;

  atom = topo_mol_get_atom(mol,target,0);
  if ( ! atom ) return -2;

  atom->x = 0;
  atom->y = 0;
  atom->z = 0;
  atom->xyz_state = TOPO_MOL_XYZ_VOID;

  return 0;
}


int topo_mol_guess_xyz(topo_mol *mol) {
  char msg[128];
  int iseg,nseg,ires,nres,ucount,i,nk,nu,gcount,gwild,okwild,wcount,hcount;
  int ipass;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom, *a1, *a2, *a3;
  topo_mol_atom_t *ka[4];
  topo_mol_atom_t *ua[4];
  topo_mol_bond_t *bondtmp;
  topo_mol_angle_t *angletmp;

#if !defined(NEWPSFGEN)
#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#endif

#if 0 && defined(NEWPSFGEN)
  /* Not currently in use, but should replace all conversions bellow*/
  const double DegreetoRad = M_PI/180.0;
#endif
  
#if !defined(NEWPSFGEN)  
  double dihedral, angle234, dist34;
  double tx, ty, tz;
#else
  double dihedral, angle, dist;
  topo_mol_atom_t *tmplp;
  int lpcount, bisector;
#endif

  topo_mol_atom_t **uatoms;
  topo_mol_conformation_t *conf;
  double r12x,r12y,r12z,r12,r23x,r23y,r23z,r23,ix,iy,iz,jx,jy,jz,kx,ky,kz;
  double a,b,c;

  if ( ! mol ) return -1;

  ucount = 0;
  hcount = 0; /* non-hydrogen atom guessed coordinates count  */

#if defined(NEWPSFGEN)
  lpcount = 0;
#endif  
  
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)  

      for ( atom = res->atoms; atom; atom = atom->next ) {
        if ( atom->xyz_state != TOPO_MOL_XYZ_SET ) {
          ++ucount;
          if ( atom->mass > 2.5 ) ++hcount;
        }
      }

#else
      for ( i = 0; i < res->atomSize; i++ ) {
        if (res->atomArray[i]->xyz_state != TOPO_MOL_XYZ_SET && 
            !res->atomArray[i]->isdrudlonepair) {
          ++ucount;
          if ( res->atomArray[i]->mass > 2.5 ) ++hcount;
        } else if (res->atomArray[i]->xyz_state != TOPO_MOL_XYZ_SET && 
                   res->atomArray[i]->isdrudlonepair) {
          ++lpcount;
        }
      }
      
#endif

    }
  }
  sprintf(msg,"Info: guessing coordinates for %d atoms (%d non-hydrogen)",
						ucount, hcount);
  topo_mol_log_error(mol,msg);

  uatoms = (topo_mol_atom_t**) malloc(ucount*sizeof(topo_mol_atom_t*));
  if ( ! uatoms ) return -2;
  ucount = 0;

#if defined(NEWPSFGEN)

  if (lpcount) {
    sprintf(msg,"Info: guessing coordinates for %d lonepairs",
              lpcount);
    topo_mol_log_error(mol,msg);    
    
  }

#endif  

  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)  

      for ( atom = res->atoms; atom; atom = atom->next ) {
        if ( atom->xyz_state != TOPO_MOL_XYZ_SET ) uatoms[ucount++] = atom;
      }

#else

      for ( i = 0; i < res->atomSize; i++) {
        if (res->atomArray[i]->xyz_state != TOPO_MOL_XYZ_SET && 
            !res->atomArray[i]->isdrudlonepair) {
            
            uatoms[ucount++] = res->atomArray[i];
        } 
      }
      
#endif

    }
  }

  for ( i=0; i<ucount; ++i ) uatoms[i]->xyz_state = TOPO_MOL_XYZ_VOID;

  /* everything below based on atom 4 unknown, all others known */

  /* from the CHARMM docs:

    Normal IC table entry:
                I
                 \
                  \
                   J----K
                         \
                          \
                           L
        values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)

    Improper type of IC table entry:
                I        L
                 \     /
                  \   /
                   *K
                   |
                   |
                   J
        values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)

  */

  gcount = 1; /* Total guesses count  */
  okwild = 0; /* Atom not guessed?!  */
  wcount = 0; /* total atom with coordinates guessed count  */
  hcount = 0; /* non-hydrogen atom with coordinates guessed count  */
  while ( gcount || ! okwild ) {
    if ( gcount == 0 ) { 
      if ( okwild ) break; 
      else okwild = 1; 
    }
    gcount = 0;
    for ( i=0; i<ucount; ++i ) { 
      atom = uatoms[i];
      if ( atom->xyz_state != TOPO_MOL_XYZ_VOID ) continue;
      /* Assign coordinates based on the IC table
       * The * in the IC definiton signals an improper angle 
       * instead of a dihedral angle
      */
      for ( conf = atom->conformations; conf;
		      conf = topo_mol_conformation_next(conf,atom) ) {
        if ( conf->del ) continue;
        /* Since the conformation data structure stores the info in 
         * duplicate, it is necessary check the order of ther atoms[0], atoms[1]
         * ,atoms[2] and atoms[3], to match a1,a2 and a3 and to match
         * the distance and the angle to be used
        */
        else if ( conf->atom[0] == atom &&
      		conf->atom[1]->xyz_state != TOPO_MOL_XYZ_VOID &&
      		conf->atom[2]->xyz_state != TOPO_MOL_XYZ_VOID &&
      		conf->atom[3]->xyz_state != TOPO_MOL_XYZ_VOID ) {
          if ( conf->improper ) {
            a1 = conf->atom[3]; a2 = conf->atom[1]; a3 = conf->atom[2];
			
#if !defined(NEWPSFGEN)
            dist34 = conf->dist12;
            angle234 = conf->angle123 * (M_PI/180.0);
#else
            dist = conf->dist12;
            angle = conf->angle123 * (M_PI/180.0);
#endif
            dihedral = -1.0 * conf->dihedral * (M_PI/180.0);
          } else {
            a1 = conf->atom[3]; a2 = conf->atom[2]; a3 = conf->atom[1];
            
#if !defined(NEWPSFGEN) 
          dist34 = conf->dist12;
          angle234 = conf->angle123 * (M_PI/180.0);
	
#else
            dist = conf->dist12;
            angle = conf->angle123 * (M_PI/180.0);
#endif
            dihedral = conf->dihedral * (M_PI/180.0);
          } 
        } 
        else if ( conf->atom[3] == atom &&
      		conf->atom[2]->xyz_state != TOPO_MOL_XYZ_VOID &&
      		conf->atom[1]->xyz_state != TOPO_MOL_XYZ_VOID &&
      		conf->atom[0]->xyz_state != TOPO_MOL_XYZ_VOID ) {
          if ( conf->improper ) {
            a1 = conf->atom[0]; a2 = conf->atom[1]; a3 = conf->atom[2];
#if !defined(NEWPSFGEN) 
            dist34 = conf->dist34;
            angle234 = conf->angle234 * (M_PI/180.0);
#else
            dist = conf->dist34;
            angle = conf->angle234 * (M_PI/180.0);
#endif
            dihedral = conf->dihedral * (M_PI/180.0);
          } else {
            a1 = conf->atom[0]; a2 = conf->atom[1]; a3 = conf->atom[2];

#if !defined(NEWPSFGEN)
            dist34 = conf->dist34;
            angle234 = conf->angle234 * (M_PI/180.0);
#else
            dist = conf->dist34;
            angle = conf->angle234 * (M_PI/180.0);
#endif
            dihedral = conf->dihedral * (M_PI/180.0);
          } 
        } 
        else continue;

        gwild = 0;/* the vector cross product is 0, so it will fail the guess*/

#if !defined(NEWPSFGEN)
        if ( dist34 == 0.0 ) { dist34 = 1.0; gwild = 1; }
        if ( angle234 == 0.0 ) { angle234 = 109.0*M_PI/180.0; gwild = 1; }
#else
        if ( dist == 0.0 ) { dist = 1.0; gwild = 1; }
        if ( angle == 0.0 ) { angle = 109.0*M_PI/180.0; gwild = 1; }
#endif

#if !defined(NEWPSFGEN)        
        /* Compute the cross vector formed by the two vectores 
         * a2->a1 (r12) and a3->a2(r23) as following:
         * (r12 x r23) x r23
         * a is the amplitute of the cross vectores to normalize the 
         * standard basis vectors i j and k components of the vector 
         * V = Vx*i + Vy*j + Vz*K
         * a = 1.0/a as i,j and 0 <= K <= 1
        */
        
        /* Compute all x,y and z components of the vectores */
        r12x = a2->x - a1->x;
        r12y = a2->y - a1->y;
        r12z = a2->z - a1->z;
        
        r23x = a3->x - a2->x;
        r23y = a3->y - a2->y;
        r23z = a3->z - a2->z;
        
        /* Computing the x, y and z components of  ix, iy and iz of 
         * r23 needed to compute  
         * the final corss product 
        */
        
        a = sqrt(r23x*r23x + r23y*r23y + r23z*r23z); 
        if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
        ix = a * r23x;
        iy = a * r23y;
        iz = a * r23z;
      
        /* "Xyzzy (mnemonic)"way to calculate the cross product between r12
         * and r23
         * k = r12 x r23
         * t = (r12 x r23) x r23; t is used and reused, but at the end 
         * stores the final cross term
         */
        tx = r12y*r23z - r12z*r23y;
        ty = r12z*r23x - r12x*r23z;
        tz = r12x*r23y - r12y*r23x;
        
        /* Computing the  x, y and z components of 
         * r12 x r23
        */
        a = sqrt(tx*tx + ty*ty + tz*tz);
        if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
        kx = a * tx;
        ky = a * ty;
        kz = a * tz;
      
        /* Computing the cross product of the vector (r12 x r23) x r23 */
        tx = ky*iz - kz*iy;
        ty = kz*ix - kx*iz;
        tz = kx*iy - ky*ix;
        
        /* Computing the x, y and z components of the final 
         * (r12 x r23) x r23 needed to compute  
         * the final corss product 
        */
        a = sqrt(tx*tx + ty*ty + tz*tz);
        if ( a == 0.0 ) gwild = 1; else a = 1.0 / a;
        jx = a * tx;
        jy = a * ty;
        jz = a * tz;
        
#if !defined(NEWPSFGEN)
        a = -1.0 * dist34 * cos(angle234);
        b = dist34 * sin(angle234) * cos(dihedral);
        c = dist34 * sin(angle234) * sin(dihedral);
#else  
        a = -1.0 * dist * cos(angle);
        b = dist * sin(angle) * cos(dihedral);
        c = dist * sin(angle) * sin(dihedral);
#endif

        if ( gwild && ! okwild ) continue;
        if ( okwild ) {
          ++wcount;
          if ( atom->mass > 2.5 ) ++hcount;
        }
        
        atom->x = a3->x + a * ix + b * jx + c * kx;
        atom->y = a3->y + a * iy + b * jy + c * ky;
        atom->z = a3->z + a * iz + b * jz + c * kz;

#else
      gwild = topo_mol_set_IC(atom, a1, a2, a3, dist, angle, dihedral, 0);
      if ( gwild && ! okwild ) continue;
      if ( okwild ) {
        ++wcount;
        if ( atom->mass > 2.5 ) ++hcount;
      }
#endif
      
        atom->xyz_state = okwild ? TOPO_MOL_XYZ_BADGUESS : TOPO_MOL_XYZ_GUESS; 
        ++gcount;
        break;  /* don't re-guess this atom */
      }
    }
  }
  
  /* look for bad angles due to swapped atom names */
  for ( i=0; i<ucount; ++i ) { 
    atom = uatoms[i];
    /* only look for errors in guessed atoms */
    if ( atom->xyz_state == TOPO_MOL_XYZ_VOID ||
         atom->xyz_state == TOPO_MOL_XYZ_SET ) continue;

    for ( angletmp = atom->angles; angletmp;
		    angletmp = topo_mol_angle_next(angletmp,atom) ) {
      if ( angletmp->del ) continue;
      if ( angletmp->atom[0] == atom ) {
        a1 = angletmp->atom[2];
        a2 = angletmp->atom[1];
      } else if ( angletmp->atom[2] == atom ) {
        a1 = angletmp->atom[0];
        a2 = angletmp->atom[1];
      } else continue;
      /* only use set atoms, don't hid topology file errors */
      if ( a1->xyz_state != TOPO_MOL_XYZ_SET ) continue;
      if ( a2->xyz_state != TOPO_MOL_XYZ_SET ) continue;

      r12x = a2->x - a1->x;
      r12y = a2->y - a1->y;
      r12z = a2->z - a1->z;
      r12 = sqrt(r12x*r12x + r12y*r12y + r12z*r12z);
      r23x = atom->x - a2->x;
      r23y = atom->y - a2->y;
      r23z = atom->z - a2->z;
      r23 = sqrt(r23x*r23x + r23y*r23y + r23z*r23z);
      /* assume wrong if angle is less than 45 degrees 
       * A · B = AB cosθ 
       * dot product of r12 . r23 < |r12|*|r23|* -cos(45) 
       * -cos(45) = -0.70
       */
      if ( r12x*r23x + r12y*r23y + r12z*r23z < r12 * r23 * -0.7 ) {
        sprintf(msg, "Warning: failed to guess coordinate due to bad angle %s %s %s",
            a1->name, a2->name, atom->name);
        topo_mol_log_error(mol, msg);
        if ( atom->xyz_state == TOPO_MOL_XYZ_BADGUESS ) {
          --wcount;
          if ( atom->mass > 2.5 ) --hcount;
        }
        --gcount;
        atom->xyz_state = TOPO_MOL_XYZ_VOID; 
        break;
      }
    }
  }

  /* fallback rules for atoms without conformation records */
  for ( ipass=0; ipass<2; ++ipass ) {  /* don't do entire chain */
    for ( i=0; i<ucount; ++i ) { atom = uatoms[i];
      if ( atom->xyz_state != TOPO_MOL_XYZ_VOID ) continue;

      /* pick heaviest known atom we are bonded to (to deal with water) */
      a1 = 0;
      for ( bondtmp = atom->bonds; bondtmp;
  		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
        if ( bondtmp->atom[0] == atom ) a2 = bondtmp->atom[1];
        else a2 = bondtmp->atom[0];
        if ( a2->xyz_state == TOPO_MOL_XYZ_VOID ) continue;
        if ( a1 == 0 || a2->mass > a1->mass ) a1 = a2;
      }
      if ( a1 == 0 ) continue;
      atom = a1;

      /* find all bonded atoms known and unknown coordinates */
      nk = 0;  nu = 0;
      for ( bondtmp = atom->bonds; bondtmp;
  		bondtmp = topo_mol_bond_next(bondtmp,atom) ) {
        if ( bondtmp->del ) continue;
        if ( bondtmp->atom[0] == atom ) a2 = bondtmp->atom[1];
        else a2 = bondtmp->atom[0];
        if ( a2->xyz_state == TOPO_MOL_XYZ_VOID ) {
          if ( nu < 4 ) ua[nu++] = a2;
        } else {
          if ( nk < 4 ) ka[nk++] = a2;
        }
      }

      if ( ipass ) {  /* hydrogens only on second pass */
        int j;
        for ( j=0; j<nu && ua[j]->mass < 2.5; ++j );
        if ( j != nu ) continue;
      }

      if ( nu + nk > 4 ) continue;  /* no intuition beyond this case */

      if ( nk == 0 ) {  /* not bonded to any known atoms */
        a1 = ua[0];
        a1->x = atom->x + 1.0;
        a1->y = atom->y;
        a1->z = atom->z;
        a1->xyz_state = TOPO_MOL_XYZ_BADGUESS;
        ++gcount;  ++wcount;
        if ( a1->mass > 2.5 ) ++hcount;
        continue;
      }

      if ( nk == 1 ) {  /* bonded to one known atom */
        a1 = ka[0];
        ix = a1->x - atom->x;
        iy = a1->y - atom->y;
        iz = a1->z - atom->z;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        jx = -1.0 * iy;  jy = ix;  jz = 0;
        if ( jx*jx + jy*jy + jz*jz < 0.1 ) {
          jx = 0;  jy = -1.0 * iz;  jz = iy;
        }
        a = sqrt(jx*jx+jy*jy+jz*jz);
        if ( a ) a = 1.0 / a;  else continue;
        jx *= a; jy *= a; jz *= a;
        if ( nu == 1 ) {  /* one unknown atom */
          a = cos(109.0*M_PI/180.0);
          b = sin(109.0*M_PI/180.0);
          a2 = ua[0];
          a2->x = atom->x + a * ix + b * jx;
          a2->y = atom->y + a * iy + b * jy;
          a2->z = atom->z + a * iz + b * jz;
          a2->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a2->mass > 2.5 ) ++hcount;
        } else if ( nu == 2 ) {  /* two unknown atoms */
          a = cos(120.0*M_PI/180.0);
          b = sin(120.0*M_PI/180.0);
          a1 = ua[0];
          a2 = ua[1];
          a1->x = atom->x + a * ix + b * jx;
          a1->y = atom->y + a * iy + b * jy;
          a1->z = atom->z + a * iz + b * jz;
          a2->x = atom->x + a * ix - b * jx;
          a2->y = atom->y + a * iy - b * jy;
          a2->z = atom->z + a * iz - b * jz;
          a1->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a1->mass > 2.5 ) ++hcount;
          a2->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a2->mass > 2.5 ) ++hcount;
        } else { /* three unknown atoms */
          a1 = ua[0];
          a2 = ua[1];
          a3 = ua[2];
          /* only handle this case if at least two are hydrogens */
          if ( a1->mass > 2.5 && a2->mass > 2.5 ) continue;
          if ( a1->mass > 2.5 && a3->mass > 2.5 ) continue;
          if ( a2->mass > 2.5 && a3->mass > 2.5 ) continue;
          kx = iy*jz - iz*jy;
          ky = iz*jx - ix*jz;
          kz = ix*jy - iy*jx;
          a = sqrt(kx*kx+ky*ky+kz*kz);
          if ( a ) a = 1.0 / a;  else continue;
          kx *= a; ky *= a; kz *= a;
          a = cos(109.0*M_PI/180.0);
          b = sin(109.0*M_PI/180.0);
          a1->x = atom->x + a * ix + b * jx;
          a1->y = atom->y + a * iy + b * jy;
          a1->z = atom->z + a * iz + b * jz;
          c = b * sin(120.0*M_PI/180.0);
          b *= cos(120.0*M_PI/180.0);
          a2->x = atom->x + a * ix + b * jx + c * kx;
          a2->y = atom->y + a * iy + b * jy + c * ky;
          a2->z = atom->z + a * iz + b * jz + c * kz;
          a3->x = atom->x + a * ix + b * jx - c * kx;
          a3->y = atom->y + a * iy + b * jy - c * ky;
          a3->z = atom->z + a * iz + b * jz - c * kz;
          a1->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a1->mass > 2.5 ) ++hcount;
          a2->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a2->mass > 2.5 ) ++hcount;
          a3->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a3->mass > 2.5 ) ++hcount;
        }
        continue;
      }

      if ( nk == 2 ) {  /* bonded to two known atoms */
        a1 = ka[0];
        ix = a1->x - atom->x;
        iy = a1->y - atom->y;
        iz = a1->z - atom->z;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        jx = ix;  jy = iy;  jz = iz;
        a1 = ka[1];
        ix = a1->x - atom->x;
        iy = a1->y - atom->y;
        iz = a1->z - atom->z;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        kx = jx - ix;  ky = jy - iy;  kz = jz - iz;
        jx += ix;  jy += iy;  jz += iz;
        a = sqrt(jx*jx+jy*jy+jz*jz);
        if ( a ) a = 1.0 / a;  else continue;
        jx *= a; jy *= a; jz *= a;
        if ( nu == 1 ) {  /* one unknown atom */
          a2 = ua[0];
          a2->x = atom->x - jx;
          a2->y = atom->y - jy;
          a2->z = atom->z - jz;
          a2->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a2->mass > 2.5 ) ++hcount;
        } else {  /* two unknown atoms */
          a1 = ua[0];
          a2 = ua[1];
          /* only handle this case if both are hydrogens */
          if ( a1->mass > 2.5 || a2->mass > 2.5 ) continue;
          a = sqrt(kx*kx+ky*ky+kz*kz);
          if ( a ) a = 1.0 / a;  else continue;
          kx *= a; ky *= a; kz *= a;
          ix = jy*kz - jz*ky;
          iy = jz*kx - jx*kz;
          iz = jx*ky - jy*kx;
          a = sqrt(ix*ix+iy*iy+iz*iz);
          if ( a ) a = 1.0 / a;  else continue;
          ix *= a; iy *= a; iz *= a;

#if !defined(NEWPSFGEN)
          angle234 = (180.0-0.5*109.0)*M_PI/180.0;
          a = sin(angle234);
          b = cos(angle234);
#else
          angle = (180.0-0.5*109.0)*M_PI/180.0;
          a = sin(angle);
          b = cos(angle);
#endif
          
          a1->x = atom->x + a * ix + b * jx;
          a1->y = atom->y + a * iy + b * jy;
          a1->z = atom->z + a * iz + b * jz;
          a2->x = atom->x - a * ix + b * jx;
          a2->y = atom->y - a * iy + b * jy;
          a2->z = atom->z - a * iz + b * jz;
          a1->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a1->mass > 2.5 ) ++hcount;
          a2->xyz_state = TOPO_MOL_XYZ_BADGUESS;
          ++gcount;  ++wcount;
          if ( a2->mass > 2.5 ) ++hcount;
        }
        continue;
      }

      if ( nk == 3 ) {  /* bonded to three known atoms */
        a1 = ka[0];
        ix = a1->x - atom->x;
        iy = a1->y - atom->y;
        iz = a1->z - atom->z;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        jx = ix;  jy = iy;  jz = iz;
        a1 = ka[1];
        ix = a1->x - atom->x;
        iy = a1->y - atom->y;
        iz = a1->z - atom->z;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        jx += ix;  jy += iy;  jz += iz;
        a1 = ka[2];
        ix = a1->x - atom->x;
        iy = a1->y - atom->y;
        iz = a1->z - atom->z;
        a = sqrt(ix*ix+iy*iy+iz*iz);
        if ( a ) a = 1.0 / a;  else continue;
        ix *= a; iy *= a; iz *= a;
        jx += ix;  jy += iy;  jz += iz;
        a = sqrt(jx*jx+jy*jy+jz*jz);
        if ( a ) a = 1.0 / a;  else continue;
        a2 = ua[0];
        a2->x = atom->x - a * jx;
        a2->y = atom->y - a * jy;
        a2->z = atom->z - a * jz;
        a2->xyz_state = TOPO_MOL_XYZ_BADGUESS;
        ++gcount;  ++wcount;
        if ( a2->mass > 2.5 ) ++hcount;
        continue;
      }

    }
  }


  gcount = 0;
  for ( i=0; i<ucount; ++i ) {
    if ( uatoms[i]->xyz_state == TOPO_MOL_XYZ_VOID ) ++gcount;
  }
  if ( wcount ) {
    sprintf(msg,"Warning: poorly guessed coordinates for %d atoms (%d non-hydrogen):", wcount, hcount);
    topo_mol_log_error(mol,msg);
    for ( iseg=0; iseg<nseg; ++iseg ) {
      seg = mol->segment_array[iseg];
      if (! seg) continue;
      nres = hasharray_count(seg->residue_hash);
      for ( ires=0; ires<nres; ++ires ) {
        res = &(seg->residue_array[ires]);


#if !defined(NEWPSFGEN)  

        for ( atom = res->atoms; atom; atom = atom->next ) {
          if ( atom->xyz_state == TOPO_MOL_XYZ_BADGUESS) {
            sprintf(msg, "Warning: poorly guessed coordinate for atom %s\t %s:%s\t  %s",
                atom->name, res->name, res->resid, seg->segid);

#else

        for ( i = 0; i < res->atomSize; i++ ) {
          if ( res->atomArray[i]->xyz_state == TOPO_MOL_XYZ_BADGUESS && 
              !res->atomArray[i]->isdrudlonepair) {
            sprintf(msg, "Warning: poorly guessed coordinate for atom %s\t %s:%s\t  %s",
                res->atomArray[i]->name, res->name, res->resid, seg->segid);
#endif

            topo_mol_log_error(mol, msg);
          }
        }
      }
    }
  }
  if ( gcount ) {
    sprintf(msg,"Warning: failed to guess coordinates for %d atoms",gcount);
    topo_mol_log_error(mol,msg);
  }

  
#if defined(NEWPSFGEN)

  /* Now that the coordinates of all atoms were set, set the lonepairs */
  
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
        
      for ( i = 0; i < res->atomSize; i++ ) {
        if (res->atomArray[i]->xyz_state != TOPO_MOL_XYZ_GUESS && 
            res->atomArray[i]->isdrudlonepair) {
          tmplp = res->atomArray[i];
          if (tmplp->lonepair->lptype == COLINEARLP) {
            if (topo_mol_set_COLINEARLP(mol, tmplp)) {
              sprintf(msg,"Warning: poorly guessed coordinates for %s lonepair", 
              tmplp->name);
              topo_mol_log_error(mol,msg);
            } 
            
            tmplp->xyz_state = TOPO_MOL_XYZ_SET;
            
          } else {
            bisector = 0;
            if (tmplp->lonepair->lptype == BISECTORLP ) {
              bisector = 1;
            }
            if (topo_mol_set_IC(tmplp->lonepair->atoms[0], tmplp->lonepair->atoms[3],
              tmplp->lonepair->atoms[2],tmplp->lonepair->atoms[1], tmplp->lonepair->distance,
              tmplp->lonepair->angle, tmplp->lonepair->dihedral, bisector)) {
                sprintf(msg,"Warning: poorly guessed coordinates for %s lonepair", 
                tmplp->name);
                topo_mol_log_error(mol,msg);
            } 
            
            tmplp->xyz_state = TOPO_MOL_XYZ_SET;
            
          }
        } else if (res->atomArray[i]->alpha) {
          res->atomArray[i]->dxyz = (double *)malloc(3 * sizeof(double));
          res->atomArray[i]->dxyz[0] = res->atomArray[i]->x;
          res->atomArray[i]->dxyz[1] = res->atomArray[i]->y;
          res->atomArray[i]->dxyz[2] = res->atomArray[i]->z;
        }
      }
    }
  }
#endif
    free((void*)uatoms);
  return 0;
}


/* Copied and modified from topo_mol_segment */
int topo_mol_add_patch(topo_mol *mol, const char *pname, int deflt) {
  topo_mol_patch_t **patches;
  topo_mol_patch_t *patchtmp;
  if ( ! mol ) return -1;
  if ( NAMETOOLONG(pname) ) return -2;
  patches = &(mol->patches);
  
  patchtmp = 0;
  patchtmp = memarena_alloc(mol->arena,sizeof(topo_mol_patch_t));
  if ( ! patchtmp ) return -3;
  
  strcpy(patchtmp->pname,pname);
  patchtmp->patchresids = 0;

  patchtmp->npres = 0;
  patchtmp->deflt = deflt;
  patchtmp->next = 0;
/*    printf("add_patch %i %s;\n", mol->npatch, patchtmp->pname);   */

  if (mol->npatch==0) {
    *patches = patchtmp;
  } else {
    mol->curpatch->next = patchtmp;
  }
  mol->curpatch = patchtmp;

  mol->npatch++;
  return 0;
}


/* Copied and modified from topo_mol_residue */
int topo_mol_add_patchres(topo_mol *mol, const topo_mol_ident_t *target) {
  topo_mol_patch_t *patch;
  topo_mol_patchres_t **patchres;
  topo_mol_patchres_t *patchrestmp;
  if ( ! mol ) return -1;
  if ( NAMETOOLONG(target->segid) ) return -2;
  if ( NAMETOOLONG(target->resid) ) return -2;

  patch = mol->curpatch; 
  patchres = &(patch->patchresids);
  patchrestmp = 0;
  patchrestmp = memarena_alloc(mol->arena,sizeof(topo_mol_patchres_t));
  if ( ! patchrestmp ) return -3;

  strcpy(patchrestmp->segid,target->segid);
  strcpy(patchrestmp->resid,target->resid);
/*   printf("add_patchres %i %s:%s;\n", patch->npres, patchrestmp->segid, patchrestmp->resid);  */
  patch->npres++;
  /* patchrestmp->next = *patchres;  old code builds list in reverse order */
  patchrestmp->next = NULL;
  while ( *patchres ) { patchres = &((*patchres)->next); }
  *patchres = patchrestmp;
  return 0;
}


/* Test the existence of segid:resid for the patch */
int topo_mol_validate_patchres(topo_mol *mol, const char *pname, const char *segid, const char *resid) {
  topo_mol_ident_t target;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  target.segid = segid;
  target.resid = resid;
  seg = topo_mol_get_seg(mol,&target);
  if ( ! seg ) {
    char errmsg[50];
    sprintf(errmsg,"Segment %s not exsisting, skipping patch %s.\n",segid,pname);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }
  res = topo_mol_get_res(mol,&target,0);
  if ( ! res ) {
    char errmsg[50];
    sprintf(errmsg,"Residue %s:%s not exsisting, skipping patch %s.\n",segid,resid,pname);
    topo_mol_log_error(mol,errmsg);
    return 0;
  }
  return 1;
}
