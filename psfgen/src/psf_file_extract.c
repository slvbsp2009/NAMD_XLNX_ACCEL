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
 *      $RCSfile: psf_file_extract.c,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.42 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "psf_file.h"
#include "psf_file_extract.h"
#include "pdb_file.h"
#include "topo_mol_struct.h"

#if defined(NEWPSFGEN)
#include "topo_defs.h"
#endif
/* General note: in a few places I read various arrays in reverse order. 
   That's because I want psf files emitted by psfgen to have all the atoms,
   bonds, etc. in the same order as in the original psf file.  We have to 
   reverse it because adding to a linked list reverses the order.  Actually,
   if the original psf file comes from some other program, then we might 
   change the order of the bonds, angles, etc. but we can at least guarantee
   that if we read a psf file written by psfgen, then write it out again,
   the output will match the input exactly.
*/

#if !defined(NEWPSFGEN)
/* Read in all psf atom information using this struct */
struct psfatom {
  char name[10];
  char atype[10];
  char resname[10];
  char segname[10];
  char resid[10];
  char element[3];
  double charge, mass;
};
typedef struct psfatom psfatom;

#else 

static topo_mol_segment_t *get_segment(topo_mol *mol, const char *segname);

static topo_mol_residue_t *get_residue(topo_mol_segment_t *seg, 
        const char *resid);

#endif


#define PSF_RECORD_LENGTH 	200

static int extract_patches(FILE *file, topo_mol *mol) { 
  char inbuf[PSF_RECORD_LENGTH+2];
  int npatch = 0;
  
  /* Read comments; get patch info */
  while (!feof(file)) {
    if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file)) {
      /* EOF with no NATOM */
      return -1;  
    }
    if (strlen(inbuf) > 0) {
      if (strstr(inbuf, "REMARKS")) {
	char *pbuf;
	if (strstr(inbuf, "REMARKS topology ")) {
	  char topofile[256];
	  pbuf = strstr(inbuf, "topology");
	  pbuf=pbuf+strlen("topology");
	  sscanf(pbuf, "%s", topofile);	  
	  topo_defs_add_topofile(mol->defs, topofile);
	}
	if (strstr(inbuf, "REMARKS patch ") || strstr(inbuf, "REMARKS defaultpatch ")) {
	  char pres[NAMEMAXLEN], segres[2*NAMEMAXLEN];
	  char s[NAMEMAXLEN], r[NAMEMAXLEN];
	  topo_mol_ident_t target;
	  pbuf = strstr(inbuf, "patch");
	  pbuf=pbuf+5;
	  sscanf(pbuf, "%s", pres);
	  if (strcmp(pres,"----")) {
	    if (strstr(inbuf, "REMARKS defaultpatch")) {
	      topo_mol_add_patch(mol,pres,1);
	    } else {
	      topo_mol_add_patch(mol,pres,0);
	    }
	  }
	  pbuf = strstr(pbuf, pres)+strlen(pres);
	  while (sscanf(pbuf, "%s", segres)==1) { 
	    int slen;
	    slen = strcspn(segres,":");
	    strncpy(s, segres, slen);
	    s[slen] = '\0';
	    strcpy(r, strchr(segres,':')+1);
	    target.segid = s;
	    target.resid = r;
	    topo_mol_add_patchres(mol,&target);
	    pbuf = strstr(pbuf,segres)+strlen(segres);
	  }
	  npatch++;
	}
      } else {
	if (strstr(inbuf, "NATOM")) {
          rewind(file);
	  return npatch;
        }
      }
    }
  } ;
  return npatch;
}


static int extract_segment_extra_data(FILE *file, topo_mol *mol) {
  char inbuf[PSF_RECORD_LENGTH+2];
  
  /* Read comments; get patch info */
  while (!feof(file)) {
    if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file)) {
      /* EOF with no NATOM */
      return -1;  
    }
    if (strlen(inbuf) > 0) {
      if (strstr(inbuf, "REMARKS")) {
	char *pbuf;
	if (strstr(inbuf, "REMARKS segment ")) {
	  char segid[NAMEMAXLEN], pfirst[NAMEMAXLEN], plast[NAMEMAXLEN];
	  char angles[20], diheds[20], tmp[NAMEMAXLEN];
	  topo_mol_segment_t *seg = NULL;
	  int id;
	  pbuf = strstr(inbuf, "segment");
	  pbuf += strlen("segment");
	  sscanf(pbuf, "%s %s %s %s %s %s %s %s %s", segid, tmp, tmp, pfirst, tmp, plast, tmp, angles, diheds);
	  if ( (id = hasharray_index(mol->segment_hash, segid)) != HASHARRAY_FAIL) {
	    /* Then the segment exists.  Look it up and return it. */
	    seg = mol->segment_array[id];
	    strcpy(strchr(pfirst,';'),"");
	    strcpy(strchr(plast, ';'),"");
	    strcpy(seg->pfirst,pfirst);
	    strcpy(seg->plast, plast);
	    seg->auto_angles = 0; 
	    if (!strcmp(angles,"angles")) {
	      seg->auto_angles = 1; 
	    }
	    seg->auto_dihedrals = 0; 
	    if (!strcmp(diheds,"dihedrals")) {
	      seg->auto_dihedrals = 1; 
	    }
	  } 
	}
      }
    }
  }
  return 0;
}

static int extract_bonds(FILE *file, int fw, topo_mol *mol, int natoms, 
                         topo_mol_atom_t **molatomlist) {

  int *bonds;
  int i, nbonds;

  /* Build bonds */
  nbonds = psf_start_block(file, "NBOND");
  if (nbonds < 0) {
    return -1; 
  }
  bonds = (int *)malloc(2*nbonds*sizeof(int));

  if (psf_get_bonds(file, fw, nbonds, bonds)) {
    free(bonds);
    return -1;
  }
 
  for (i=nbonds-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2;
    topo_mol_bond_t *tuple;
    int ind1, ind2; 
  
    ind1 = bonds[2*i]-1; 
    ind2 = bonds[2*i+1]-1;
    if (ind1 < 0 || ind2 < 0 || ind1 >= natoms || ind2 >= natoms) {
      /* Bad indices, abort now */
      free(bonds);
      return -1;
    }
    
    atom1 = molatomlist[ind1];
    atom2 = molatomlist[ind2];


#if defined(NEWPSFGEN)
    /* Skip the lone pairs. The bonds are restored if the VPBONDS are set
     * to 1 later during the written process
     */
    if (atom1->isdrudlonepair || atom2->isdrudlonepair) {
      continue;
    }
#endif

    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_bond_t));
    tuple->next[0] = atom1->bonds;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->bonds;
    tuple->atom[1] = atom2;
    tuple->del = 0;
 
    atom1->bonds = tuple; 
    atom2->bonds = tuple;

  }
  free(bonds);
  return 0;
}

static int extract_angles(FILE *file, int fw, topo_mol *mol, int natoms, 
                         topo_mol_atom_t **molatomlist) {

  int i, nangles;
  int *angles;
  
  nangles = psf_start_block(file, "NTHETA");
  if (nangles < 0) return -1; 
  angles = (int *)malloc(3*nangles*sizeof(int));

  if (psf_get_angles(file, fw, nangles, angles)) {
    free(angles); 
    return -1; 
  } 
  
  for (i=nangles-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2, *atom3;
    topo_mol_angle_t *tuple;

    atom1 = molatomlist[angles[3*i]-1];
    atom2 = molatomlist[angles[3*i+1]-1];
    atom3 = molatomlist[angles[3*i+2]-1];
   
    tuple = memarena_alloc(mol->angle_arena,sizeof(topo_mol_angle_t));
    tuple->next[0] = atom1->angles;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->angles;
    tuple->atom[1] = atom2;
    tuple->next[2] = atom3->angles;
    tuple->atom[2] = atom3;
    tuple->del = 0;
 
    atom1->angles = tuple; 
    atom2->angles = tuple;
    atom3->angles = tuple;
  }
  free(angles);
  return 0;
}

static int extract_dihedrals(FILE *file, int fw, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist) {

  int i, ndihedrals;
  int *dihedrals;

  ndihedrals = psf_start_block(file, "NPHI");
  if (ndihedrals < 0) return -1; 
  
  dihedrals = (int *)malloc(4*ndihedrals*sizeof(int));

  if (psf_get_dihedrals(file, fw, ndihedrals, dihedrals)) {
    free(dihedrals); 
    return -1;
  }
   
  for (i=ndihedrals-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2, *atom3, *atom4;
    topo_mol_dihedral_t *tuple;


    atom1 = molatomlist[dihedrals[4*i]-1];
    atom2 = molatomlist[dihedrals[4*i+1]-1];
    atom3 = molatomlist[dihedrals[4*i+2]-1];
    atom4 = molatomlist[dihedrals[4*i+3]-1];

    tuple = memarena_alloc(mol->dihedral_arena,sizeof(topo_mol_dihedral_t));
    
#if !defined(NEWPSFGEN)
    tuple->next[0] = atom1->dihedrals;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->dihedrals;
    tuple->atom[1] = atom2;
    tuple->next[2] = atom3->dihedrals;
    tuple->atom[2] = atom3;
    tuple->next[3] = atom4->dihedrals;
    tuple->atom[3] = atom4;
    tuple->del = 0;

    atom1->dihedrals = tuple;
    atom2->dihedrals = tuple;
    atom3->dihedrals = tuple;
    atom4->dihedrals = tuple;
#else
    tuple->next = atom1->dihedrals;
    tuple->atom[0] = atom2;
    tuple->atom[1] = atom3;
    tuple->atom[2] = atom4;
    tuple->del = 0;
    atom1->dihedrals = tuple;
#endif
  }
  free(dihedrals);
  return 0;
}

static int extract_impropers(FILE *file, int fw, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist) {

  int i, nimpropers;
  int *impropers;
  
  nimpropers = psf_start_block(file, "NIMPHI");
  if (nimpropers < 0) return -1; 
  impropers = (int *)malloc(4*nimpropers*sizeof(int));

  if (psf_get_impropers(file, fw, nimpropers, impropers)) {
    free(impropers); 
    return -1;
  } 
    
  for (i=nimpropers-1; i >= 0; i--) {
    topo_mol_atom_t *atom1, *atom2, *atom3, *atom4;
    topo_mol_improper_t *tuple;

    atom1 = molatomlist[impropers[4*i]-1];
    atom2 = molatomlist[impropers[4*i+1]-1];
    atom3 = molatomlist[impropers[4*i+2]-1];
    atom4 = molatomlist[impropers[4*i+3]-1];
    
    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_improper_t));

#if !defined(NEWPSFGEN)   
    
    tuple->next[0] = atom1->impropers;
    tuple->atom[0] = atom1;
    tuple->next[1] = atom2->impropers;
    tuple->atom[1] = atom2;
    tuple->next[2] = atom3->impropers;
    tuple->atom[2] = atom3;
    tuple->next[3] = atom4->impropers;
    tuple->atom[3] = atom4;
    tuple->del = 0;
 
    atom1->impropers = tuple; 
    atom2->impropers = tuple;
    atom3->impropers = tuple;
    atom4->impropers = tuple;
#else
    tuple->next = atom1->impropers;
    tuple->atom[0] = atom2;
    tuple->atom[1] = atom3;
    tuple->atom[2] = atom4;
    tuple->del = 0;

    atom1->impropers = tuple; 

#endif
  }
  free(impropers);
  return 0;
}

static int extract_cmaps(FILE *file, int fw, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist) {

  int i, j, ncmaps;
  int *cmaps;
  
  ncmaps = psf_start_block(file, "NCRTERM");
  if (ncmaps < 0) {
    return 1;
  }
  cmaps = (int *)malloc(8*ncmaps*sizeof(int));

  if (psf_get_cmaps(file, fw, ncmaps, cmaps)) {
    free(cmaps); 
    return -1;
  } 
    
  for (i=ncmaps-1; i >= 0; i--) {
    topo_mol_atom_t *atoml[8];
    topo_mol_cmap_t *tuple;

    tuple = memarena_alloc(mol->arena,sizeof(topo_mol_cmap_t));
    for ( j = 0; j < 8; ++j ) {
      atoml[j] = molatomlist[cmaps[8*i+j]-1];
      tuple->next[j] = atoml[j]->cmaps;
      tuple->atom[j] = atoml[j];
    }
    tuple->del = 0;
    for ( j = 0; j < 8; ++j ) {
      atoml[j]->cmaps = tuple; 
    }
  }

  free(cmaps);
  return 0;
}

static int extract_exclusions(FILE *file, int fw, topo_mol *mol, int natoms, 
                         topo_mol_atom_t **molatomlist) {
  int *exclusions, *exclusion_indices;
  int i, j, nexclusions;
  int exclusion_index = 0; 

  /* Read explicit exclusion list */
  nexclusions = psf_start_block(file, "NNB");
  if (nexclusions < 0) {
    return -1; 
  }
  exclusions = (int *)malloc(nexclusions*sizeof(int));
  exclusion_indices = (int *)malloc(natoms*sizeof(int));

  if (psf_get_exclusions(file, fw, nexclusions, exclusions, natoms, exclusion_indices)) {
    free(exclusions);
    free(exclusion_indices);
    return -1;
  }
 
  for (i=0; i < natoms; i++) {
    topo_mol_atom_t *atom1, *atom2;
    topo_mol_exclusion_t *excl;
    atom1 = molatomlist[i];

    if (exclusion_indices[i] > nexclusions || exclusion_indices[i] < 0) {
      free(exclusions);
      free(exclusion_indices);
      printf("lol3[%d]\n", i);
      return -1;
    }
    for (j = exclusion_indices[i]-1; j >= exclusion_index; j--) {
      int ind2 = exclusions[j] - 1;
      if (ind2 < 0 || ind2 >= natoms) {
        free(exclusions);
        free(exclusion_indices);
        return -1;
      }
      atom2 = molatomlist[ind2];
      excl = memarena_alloc(mol->arena,sizeof(topo_mol_exclusion_t));
      excl->next[0] = atom1->exclusions;
      excl->atom[0] = atom1;
      excl->next[1] = atom2->exclusions;
      excl->atom[1] = atom2;
      excl->del = 0;
      atom1->exclusions = excl; 
      atom2->exclusions = excl;
    }
    exclusion_index = exclusion_indices[i];
  }

  free(exclusions);
  free(exclusion_indices);
  return 0;
}

#if defined(NEWPSFGEN)
static int extract_lonepairs(FILE *file, int fw, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist) {

  lonepair *lpair = NULL;
  int i, j, nlpairs = 0, auxlpnum = 0, numlphosts = 0;
  int *lphosts = NULL;
  char inbuf[PSF_RECORD_LENGTH+2];
  char *hostptr=NULL;
  int chr=0;
  
  /* Build lone pairs */
  nlpairs = psf_start_block(file, "NUMLP");
  if (nlpairs < 0) {
    return -1; 
  }

  lpair = malloc(sizeof(lonepair)); /* temp variable to store the info read from 
                                     * the file
                                     */
  for (i = 0; i < natoms; i++) {
      
    if (molatomlist[i]->isdrudlonepair == ISLONEPAIR) {

      if (psf_get_lonepair_info(file, fw, lpair)) {
        return -2;
      }

      molatomlist[i]->lonepair = 
                      memarena_alloc(mol->arena,sizeof(topo_mol_lonepair_t));

      if (lpair->numhost == 2) {
        molatomlist[i]->lonepair->lptype = COLINEARLP;
        
      } else if (lpair->numhost == 3) {

        if (lpair->distance < 0) {
          molatomlist[i]->lonepair->lptype = BISECTORLP;
        } else {
          molatomlist[i]->lonepair->lptype = RELATIVELP;
        }
      }
      
      molatomlist[i]->lonepair->distance = lpair->distance;
      
      /* Covert the angles and dihedral units */
      if (molatomlist[i]->lonepair->lptype != COLINEARLP) {
        lpair->angle *= (M_PI/180.0);
      }
      molatomlist[i]->lonepair->angle = lpair->angle;
      molatomlist[i]->lonepair->dihedral = lpair->dihedral * (M_PI/180.0);
      auxlpnum++;
      if (auxlpnum == nlpairs) {
        break;
      }
    } 
  }

  /* Assign the host atoms to the lone pairs */
  auxlpnum = 0;
  lphosts = (int *)malloc(8*sizeof(int)); /* maximum number of host is 8 */

  /* Instead of having a function in the psf_file.c to extract the information
   * of the host of the lonepairs, it is better to have in here, as more 
   * control of loop is necessary. The number of hosts is not fixed, meaning that
   * the lonepairs can have 2 or 3 (or more) hosts. so we need to preserve the 
   * buffer (line) between different lonepair definition. 
   */

  for (i = 0; i < nlpairs; i++) {

    
    for (j = 0; j < 8; j++) {
      lphosts[j] = -1;
    }

    j = 0;

    do {

      if ((chr % 8) == 0) {
        /* must read next line */
        if(!fgets(inbuf,PSF_RECORD_LENGTH+2,file)) {
          /* early EOF encountered */
          break;
        }
        hostptr = inbuf;

      }

      if((lphosts[j] = psf_get_lonepair_hosts(fw, &hostptr)) < 1) {
        fprintf(stderr,"BAD LONE PAIR LINE IN PSF FILE:\n: %s\n", inbuf);
        return -2;
      }

      /* extract 1 to get the index in the array*/
      lphosts[j]--;
      if (lphosts[j] < 0 || lphosts[j] > natoms) {
        free(lphosts);
        return -4;
      }
      
      if (j == 0) {
        switch (molatomlist[lphosts[0]]->lonepair->lptype) {
          case COLINEARLP: 
            numlphosts = 3;
            break;
          case RELATIVELP:
          case BISECTORLP:
            numlphosts = 4;
            break;
          default:
            numlphosts = 4;
            break;
        }

        molatomlist[lphosts[0]]->lonepair->atoms = 
          (topo_mol_atom_t **)malloc(numlphosts*sizeof(topo_mol_atom_t*));
      } 

      molatomlist[lphosts[0]]->lonepair->atoms[j] = molatomlist[lphosts[j]];

      j++;
      chr++;
    } while (j < numlphosts);

  }
  free(lpair);
  return 0;
}

/* Get the Anisotropy information from psf
 * The anisotropy tensors are defined before the definition of the 
 * atoms involved in the anisotropy. Psfaniso stores the tensors 
 * information to then be associated with the anisotropy at the 
 * the residue level and not atom level like bonds, angles and etc.
 */
static int extract_anisotropy(FILE *file, int fw, topo_mol *mol, int natoms,
                         topo_mol_atom_t **molatomlist, psfatom *atomlist) {

  char inbuf[PSF_RECORD_LENGTH+2];
  int i, j, id=0, naniso=0;
  int *anisohost=NULL;
  psfaniso **aniso=NULL;
  topo_mol_segment_t *seg=NULL;
  topo_mol_residue_t *res=NULL;
  topo_mol_anisotropy_t *newitem=NULL, *tmpaniso=NULL;
  const char *resid=NULL, *segname=NULL;
  char *anisoptr=NULL;


    /* Build lone pairs */
  naniso = psf_start_block(file, "NUMANISO");
  if (naniso < 0) {
    return -1; 
  }

  aniso = (psfaniso **)malloc(naniso*sizeof(psfaniso*));

  anisohost = (int *)malloc(4*sizeof(int));

  if (psf_get_aniso_tensors(file, fw, aniso, naniso) == -1) {
    fprintf(stderr,"Fail to parse the anisotropy section:\n: %s\n", inbuf);
    return -1;
  }


  for (i = 0; i < naniso; i++) {

      if (i % 2 == 0) {
        if(!fgets(inbuf,PSF_RECORD_LENGTH+2,file)) {
          return -1;
        }
        anisoptr = inbuf;
      }
      if((psf_get_aniso_hosts(file, fw, &anisoptr, anisohost)) < 0) {
        fprintf(stderr,"BAD ANISOTROPY LINE IN PSF FILE:\n: %s\n", inbuf);
        return -1;
      }

      resid = atomlist[anisohost[0] -1].resid;
      segname = atomlist[anisohost[0] -1].segname;
      id = hasharray_index(mol->segment_hash, segname);
      if (id == HASHARRAY_FAIL) {
        return -1;
      }
      seg = mol->segment_array[id];
      if (!seg) {
        fprintf(stderr,"ERROR: ANISO CANNOT FIND SEGMENT:\n: %s\n", inbuf); 
        return -1;
      }

      id = hasharray_insert(seg->residue_hash, resid);
      if (id == HASHARRAY_FAIL) {
        return -1;
      }
      res = &(seg->residue_array[id]);
      if (!res) { 
        fprintf(stderr,"ERROR: ANISO CANNOT FIND RESIDUE:\n: %s\n", inbuf); 
        return -1;
      }

      newitem = memarena_alloc(mol->arena,sizeof(topo_mol_anisotropy_t));
      newitem->atoms = (topo_mol_atom_t **)malloc(4*sizeof(topo_mol_atom_t*));

      for (j = 0; j < 4; j++) {
        /* have to deduct 1 from the serial number to have the atom index */
        newitem->atoms[j] = molatomlist[anisohost[j] -1];
      }

      newitem->k11 = aniso[i]->k11;
      newitem->k22 = aniso[i]->k22;
      newitem->k33 = aniso[i]->k33;
      newitem->del = 0;
      newitem->next = 0;

      /* ensure that we keep the order of the anisotropy. Since the original
       * topology files are read from bottom to top, the raw data has to be 
       * set on the inverse order too.
       */
      if (!res->aniso) {
        res->aniso = newitem;
      } else {
        tmpaniso = res->aniso;

        do {
          if (!tmpaniso->next) {
            tmpaniso->next = newitem;
            break;
          } 
          tmpaniso = tmpaniso->next;
        } while (tmpaniso);
      }

      ++res->numaniso;
  }
  

  free(aniso);
  return 0;

}
#endif

/* Return the segment corresponding to the given segname.  If the segname
   doesn't exist, add it.  Return NULL on error.
*/
static topo_mol_segment_t *get_segment(topo_mol *mol, const char *segname) {
  int id;
  topo_mol_segment_t *seg = NULL;
  
  if ( (id = hasharray_index(mol->segment_hash, segname)) != HASHARRAY_FAIL) {
    /* Then the segment exists.  Look it up and return it. */
    seg = mol->segment_array[id];
  } else {
    /* Must create new segment */
    id = hasharray_insert(mol->segment_hash, segname);
    if (id != HASHARRAY_FAIL) {
      seg = mol->segment_array[id] =
            (topo_mol_segment_t *) malloc(sizeof(topo_mol_segment_t)); 
      strcpy(seg->segid, segname);
      seg->residue_hash = hasharray_create(
        (void**) &(seg->residue_array), sizeof(topo_mol_residue_t));
      strcpy(seg->pfirst,"");
      strcpy(seg->plast,"");
      seg->auto_angles = 0; 
      seg->auto_dihedrals = 0; 
    }
  }
  return seg;
}

/* Return a new residue with the given resid.  Add it to the given segment.
   If the resid already exists, return NULL.  Return NULL if there's a problem.
*/

static topo_mol_residue_t *get_residue(topo_mol_segment_t *seg, 
        const char *resid) {
  
  int id;
  topo_mol_residue_t *res;
  
  /* Check that the residue doesn't already exist */
  if ( hasharray_index(seg->residue_hash,resid) != HASHARRAY_FAIL ) {
    return NULL; 
  }
  id = hasharray_insert(seg->residue_hash, resid);
  if (id == HASHARRAY_FAIL) {
    return NULL;
  }
  res = &(seg->residue_array[id]);
  strcpy(res->resid, resid);
  
  return res;
}

/* UPDATE ALL Functions in here to use the res->atomArray instead of the molatomlist  */
int psf_file_extract(topo_mol *mol, FILE *file, FILE *pdbfile, FILE *namdbinfile, FILE *velnamdbinfile, 
                                void *vdata, void *v, void (*print_msg)(void *, void *, const char *)) {
  int i, natoms, charmmext;
  psfatom *atomlist;
  double *atomcoords, *atomvels;
  topo_mol_atom_t **molatomlist;

  long filepos;
  char inbuf[PSF_RECORD_LENGTH+2];

#if defined(NEWPSFGEN)
  int drude, lonepair, ndrude = 0, resnatoms = 0, h = 0;
#endif
  
  /* Read header flags */
  if (feof(file) || (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file))) {
    print_msg(vdata, v,"ERROR: Unable to read psf file");
    return -1;
  }
  if ( strncmp(inbuf, "PSF", 3) ) {
    print_msg(vdata, v,"ERROR: File does not begin with PSF - wrong format?");
    return -1;
  }
  charmmext = ( strstr(inbuf, "EXT") ? 1 : 0 ); 

#if defined(NEWPSFGEN)
  /* Is this a drude force-field structure*/
  drude = ( strstr(inbuf, "DRUDE") ? 1 : 0 ); 
#endif

  /* Read patch info from REMARKS */
  extract_patches(file, mol);

  natoms = psf_start_atoms(file);
  if (natoms < 0) {
    print_msg(vdata, v,"ERROR: Unable to read psf file");
    return -1;
  }
 
  atomlist = (psfatom *)malloc(natoms * sizeof(psfatom));

  /* Read in all atoms */
  for (i=0; i<natoms; i++) {
    psfatom *atom = atomlist + i;
    strcpy(atom->element,"");
    
#if !defined(NEWPSFGEN)
    if (psf_get_atom(file, atom->name,atom->atype,atom->resname, atom->segname,
                     atom->resid, &atom->charge, &atom->mass)
        < 0) {
#else
    if (psf_get_atom(file, atom, drude) < 0) {
#endif

      print_msg(vdata, v,"error reading atoms from psf file");
      free(atomlist);
      return -1;
    }
  }

  /* Optionally read coordinates, insertion code, and element symbol from PDB file */
  atomcoords = 0;
  if ( pdbfile ) {
    char record[PDB_RECORD_LENGTH+2];
    int indx, insertions;
    float x,y,z,o,b;
    char name[8], resname[8], chain[8];
    char segname[8], element[8], resid[8], insertion[8];

    atomcoords = (double *)malloc(natoms * 3 * sizeof(double));

    insertions = 0;
    i=0; 
    do {
      if((indx = read_pdb_record(pdbfile, record)) == PDB_ATOM) {
        psfatom *atom = atomlist + i;
        if ( i >= natoms ) {
          print_msg(vdata, v,"too many atoms in pdb file");
          free(atomlist);
          free(atomcoords);
          return -1;
        }      

        get_pdb_fields(record, name, resname, chain,
                   segname, element, resid, insertion, &x, &y, &z, &o, &b);
        if ( strncmp(atom->name,name,4) ||
             strncmp(atom->resname,resname,4) ||
             strncmp(atom->segname,segname,4) ) {
          print_msg(vdata, v,"atom mismatch in pdb file");
          print_msg(vdata, v,record);
          free(atomlist);
          free(atomcoords);
          return -1;
        }
        if ( insertion[0] != ' ' && insertion[0] != '\0' ) {
          if ( ( strlen(atom->resid ) <= 4 ) && strncmp(atom->resid,resid,4) ) {
            strncpy(atom->resid,resid,7);  atom->resid[7] = 0;
            ++insertions;
          }
          if ( atom->resid[strlen(atom->resid)-1] != insertion[0] ) {
            strncat(atom->resid,insertion,1);
          }
        }

#if defined(NEWPSFGEN)

        strncpy(atom->chain,chain,4);

#endif
        strncpy(atom->element,element,3);  atom->element[2] = 0;
        atomcoords[i*3    ] = x;
        atomcoords[i*3 + 1] = y;
        atomcoords[i*3 + 2] = z;
        ++i;
      }
    } while (indx != PDB_END && indx != PDB_EOF);
    if ( insertions ) {
      char buf[80];
      sprintf(buf, "Found %d mismatched resids with insertion codes in pdb file", insertions);
      print_msg(vdata, v,buf);
    }
    if ( i < natoms ) {
      print_msg(vdata, v,"too few atoms in pdb file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
  }

  if ( namdbinfile ) {
    int numatoms;
    int filen;
    int wrongendian;
    char lenbuf[4];
    char tmpc;

    if ( ! atomcoords ) {
      atomcoords = (double *)malloc(natoms * 3L * sizeof(double));
    }
    fseek(namdbinfile,0,SEEK_END);
    numatoms = (ftell(namdbinfile)-4)/24;
    if (numatoms < 1) {
      print_msg(vdata, v,"namdbin file is too short");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    fseek(namdbinfile,0,SEEK_SET);
    fread(&filen, sizeof(int), 1, namdbinfile);
    wrongendian = 0;
    if (filen != numatoms) {
      wrongendian = 1;
      memcpy(lenbuf, (const char *)&filen, 4);
      tmpc = lenbuf[0]; lenbuf[0] = lenbuf[3]; lenbuf[3] = tmpc;
      tmpc = lenbuf[1]; lenbuf[1] = lenbuf[2]; lenbuf[2] = tmpc;
      memcpy((char *)&filen, lenbuf, 4);
    }
    if (filen != numatoms) {
      print_msg(vdata, v,"inconsistent atom count in namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (numatoms < natoms) {
      print_msg(vdata, v,"too few atoms in namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (numatoms > natoms) {
      print_msg(vdata, v,"too many atoms in namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (wrongendian) {
      print_msg(vdata, v,"namdbin file appears to be other-endian");
    }
    if (fread(atomcoords, sizeof(double), 3L * natoms, namdbinfile)
                                 != (size_t)(3L * natoms)) {
      print_msg(vdata, v,"error reading data from namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (wrongendian) {
      int i;
      char tmp0, tmp1, tmp2, tmp3;
      char *cdata = (char *) atomcoords;
      print_msg(vdata, v,"converting other-endian data from namdbin file");
      for ( i=0; i<3*natoms; ++i, cdata+=8 ) {
        tmp0 = cdata[0]; tmp1 = cdata[1];
        tmp2 = cdata[2]; tmp3 = cdata[3];
        cdata[0] = cdata[7]; cdata[1] = cdata[6];
        cdata[2] = cdata[5]; cdata[3] = cdata[4];
        cdata[7] = tmp0; cdata[6] = tmp1;
        cdata[5] = tmp2; cdata[4] = tmp3;
      }
    }
  }

  atomvels = 0;
  if ( velnamdbinfile ) {
    int numatoms;
    int filen;
    int wrongendian;
    char lenbuf[4];
    char tmpc;
    fseek(velnamdbinfile,0,SEEK_END);
    numatoms = (ftell(velnamdbinfile)-4)/24;
    if (numatoms < 1) {
      print_msg(vdata, v,"velnamdbin file is too short");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    fseek(velnamdbinfile,0,SEEK_SET);
    fread(&filen, sizeof(int), 1, velnamdbinfile);
    wrongendian = 0;
    if (filen != numatoms) {
      wrongendian = 1;
      memcpy(lenbuf, (const char *)&filen, 4);
      tmpc = lenbuf[0]; lenbuf[0] = lenbuf[3]; lenbuf[3] = tmpc;
      tmpc = lenbuf[1]; lenbuf[1] = lenbuf[2]; lenbuf[2] = tmpc;
      memcpy((char *)&filen, lenbuf, 4);
    }
    if (filen != numatoms) {
      print_msg(vdata, v,"inconsistent atom count in namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (numatoms < natoms) {
      print_msg(vdata, v,"too few atoms in namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (numatoms > natoms) {
      print_msg(vdata, v,"too many atoms in namdbin file");
      free(atomlist);
      free(atomcoords);
      return -1;
    }
    if (wrongendian) {
      print_msg(vdata, v,"namdbin file appears to be other-endian");
    }
    atomvels = (double *)malloc(natoms * 3L * sizeof(double));
    if (fread(atomvels, sizeof(double), 3L * natoms, velnamdbinfile)
                                 != (size_t)(3L * natoms)) {
      print_msg(vdata, v,"error reading data from namdbin file");
      free(atomlist);
      free(atomcoords);
      free(atomvels);
      return -1;
    }
    if (wrongendian) {
      long i;
      char tmp0, tmp1, tmp2, tmp3;
      char *cdata = (char *) atomcoords;
      print_msg(vdata, v,"converting other-endian data from namdbin file");
      for ( i=0; i<3L*natoms; ++i, cdata+=8 ) {
        tmp0 = cdata[0]; tmp1 = cdata[1];
        tmp2 = cdata[2]; tmp3 = cdata[3];
        cdata[0] = cdata[7]; cdata[1] = cdata[6];
        cdata[2] = cdata[5]; cdata[3] = cdata[4];
        cdata[7] = tmp0; cdata[6] = tmp1;
        cdata[5] = tmp2; cdata[4] = tmp3;
      }
    }
  }
  
  molatomlist = (topo_mol_atom_t **)malloc(natoms * sizeof(topo_mol_atom_t *));

  i=0; 

#if defined(NEWPSFGEN)
  
  lonepair = 0;

#endif

  while (i < natoms) {
    topo_mol_segment_t *seg;
    topo_mol_residue_t *res;
    topo_mol_atom_t *atomtmp;
    int firstatom, j;
    const char *resid, *segname;

    resid = atomlist[i].resid;
    segname = atomlist[i].segname;
    seg = get_segment(mol, segname);
    if (!seg) { 
      print_msg(vdata, v,"ERROR: unable to get segment!");
      break;
    }
    res = get_residue(seg, resid);
    if (!res) {
      char *buf;
      int len = strlen(resid) + strlen(segname);
      buf = (char *)malloc((50 + len)*sizeof(char));
      sprintf(buf, "Unable to add (duplicate?) residue %s:%s", segname, resid);
      print_msg(vdata, v,buf);
      free(buf);
      break;
    }
    strcpy(res->name, atomlist[i].resname);
    strcpy(res->chain, "");

#if !defined(NEWPSFGEN)

    res->atoms = 0;

#else

    ndrude = 0;
      
#endif

    firstatom = i;
    while (i<natoms && !strcmp(resid, atomlist[i].resid) &&
                       !strcmp(segname, atomlist[i].segname)) {
      /* Add atoms to residue */
      atomtmp = memarena_alloc(mol->arena, sizeof(topo_mol_atom_t));
      atomtmp->bonds = 0;
      atomtmp->angles = 0;
      atomtmp->dihedrals = 0;
      atomtmp->impropers = 0;
      atomtmp->cmaps = 0;
      atomtmp->exclusions = 0;
      atomtmp->conformations = 0;
      strcpy(atomtmp->name, atomlist[i].name);
      strcpy(atomtmp->type, atomlist[i].atype);
      strcpy(atomtmp->element, atomlist[i].element);
      atomtmp->mass = atomlist[i].mass; 
      atomtmp->charge = atomlist[i].charge;
      if (atomcoords) {
        atomtmp->x = atomcoords[i*3    ];
        atomtmp->y = atomcoords[i*3 + 1];
        atomtmp->z = atomcoords[i*3 + 2];
        atomtmp->xyz_state = TOPO_MOL_XYZ_SET;
      } else {
        atomtmp->x = 0;       
        atomtmp->y = 0;       
        atomtmp->z = 0;       
        atomtmp->xyz_state = TOPO_MOL_XYZ_VOID;
      }
      if (atomvels) {
        atomtmp->vx = atomvels[i*3    ];
        atomtmp->vy = atomvels[i*3 + 1];
        atomtmp->vz = atomvels[i*3 + 2];
      } else {
        atomtmp->vx = 0;       
        atomtmp->vy = 0;       
        atomtmp->vz = 0;       
      }
      atomtmp->partition = 0;
      atomtmp->copy = 0;
      atomtmp->atomid = 0;
      
#if defined(NEWPSFGEN)
      /*initialize lonepairs and drude values as 0 */
      atomtmp->del = 0;

      /* The lone pairs have the column 9 of the psf file flagged as -1 
       * and the drude particles -2
       */

      if (atomlist[i].lpd == -1) {
        atomtmp->isdrudlonepair = ISLONEPAIR;
        lonepair= 1;
        /* flag the molecule as containing lonepairs */
        mol->lonepairs = 1;
      } else if (atomlist[i].lpd == -2) {
        atomtmp->isdrudlonepair = ISDRUDE;
      } else {
        atomtmp->isdrudlonepair = 0;
      }
      atomtmp->lonepair = NULL;
      if (drude) {
        atomtmp->alpha = atomlist[i].alpha;
        atomtmp->thole = atomlist[i].thole;
       
      } else {
        atomtmp->alpha = 0;
        atomtmp->thole = 0;
      }
      /*
       * the drude particle is always defined right after the host, so
       * once the drude particle is read, we can update the information in 
       * the host by addressing the molatomlist at the position -1
       */
      if (atomlist[i].lpd == -2) {
        molatomlist[i -1]->dcharge = atomtmp->charge;
        molatomlist[i -1]->dxyz = (double *)malloc(3 * sizeof(double));
        molatomlist[i -1]->dxyz[0] = atomtmp->x;
        molatomlist[i -1]->dxyz[1] = atomtmp->y;
        molatomlist[i -1]->dxyz[2] = atomtmp->z;
        strcpy(molatomlist[i -1]->dname, atomtmp->name);
        ndrude++;
        /* flag the molecule as containing drude particles */
        mol->drude = 1;
      } else {
        atomtmp->dcharge = 0;
      }
      
      atomtmp->dxyz = NULL;
#endif

      /* Save pointer to atom in my table so I can put in the bond 
       * information without having find the atom.
       */
      molatomlist[i] = atomtmp;
      i++;
    }

#if defined(NEWPSFGEN)

    resnatoms = i - firstatom - ndrude;
    res->atomArray = 
    (topo_mol_atom_t**) malloc((resnatoms +1)*sizeof(topo_mol_atom_t*));
    memset(res->atomArray, 0, (resnatoms +1)*sizeof(topo_mol_atom_t*));
    res->atomSize = resnatoms;
    res->reordered = 0;
    res->aniso = 0;
    res->numaniso = 0;
    res->pres[0] = '\0';
    res->lonepairs = 0;
	/* Still need to be tested in case of not having chain information */
    strcpy(res->chain, atomlist[i-1].chain);
#endif

#if !defined(NEWPSFGEN)

    for (j=i-1; j >= firstatom; j--) {
      /* Add new atoms to head of linked list in reverse order, so that
       * the linked list is in the order they appear in the psf file. 
       */
      atomtmp = molatomlist[j];

      atomtmp->next = res->atoms;
      res->atoms = atomtmp;
    }  


#else
    h = 0;
    for (j = firstatom; j < i; j++) {

      res->atomArray[h] = molatomlist[j];

      /* Skip the drude that comes after the drude host */

      if (drude && molatomlist[j]->alpha) {
        j++;
      }
      h++;
    }

#endif

  }  

  if (atomcoords) free(atomcoords);
  atomcoords = 0;
  if (atomvels) free(atomvels);
  atomvels = 0;

  /* Check to see if we broke out of the loop prematurely */
  if (i != natoms) {
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  /* Get the segment patch first,last and auto angles,dihedrals info from psf */
  /* We have to rewind the file and read the info now since it has to be added to */
  /* the existing segments which have just been read. */
  filepos = ftell(file);
  rewind(file);
  extract_segment_extra_data(file, mol);
  fseek(file, filepos, SEEK_SET);

  if (extract_bonds(file, (charmmext ? 10 : 8), mol, natoms, molatomlist)) {
    print_msg(vdata, v,"Error processing bonds");
    free(atomlist);
    free(molatomlist);
    return -1;
  }
 
  if (extract_angles(file, (charmmext ? 10 : 8), mol, natoms, molatomlist)) {
    print_msg(vdata, v,"Error processing angles");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  if (extract_dihedrals(file, (charmmext ? 10 : 8), mol, natoms, molatomlist)) {
    print_msg(vdata, v,"Error processing dihedrals");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  if (extract_impropers(file, (charmmext ? 10 : 8), mol, natoms, molatomlist)) {
    print_msg(vdata, v,"Error processing impropers");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  if (extract_exclusions(file, (charmmext ? 10 : 8), mol, natoms, molatomlist)) {
    print_msg(vdata, v,"Error processing explicit exclusions");
    free(atomlist);
    free(molatomlist);
    return -1;
  }
 

#if defined(NEWPSFGEN)
  /* Update the information for the lone pairs and drude particles */

  if ( lonepair && extract_lonepairs(file, (charmmext ? 10 : 8), mol, natoms, 
                                     molatomlist)) {
    print_msg(vdata, v,"Error processing lone pair section");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  if (drude && extract_anisotropy(file, (charmmext ? 10 : 8), mol, natoms, 
                                  molatomlist, atomlist)) {
    print_msg(vdata, v,"Error processing anisotropy section");
    free(atomlist);
    free(molatomlist);
    return -1;
  }


#endif

  switch (extract_cmaps(file, (charmmext ? 10 : 8), mol, natoms, molatomlist)) {
  case 0:
    break;
  case 1:
    print_msg(vdata, v,"psf file does not contain cross-terms");
    break;
  default:
    print_msg(vdata, v,"Error processing cross-terms");
    free(atomlist);
    free(molatomlist);
    return -1;
  }

  free(atomlist);
  free(molatomlist);
  return 0;
}


