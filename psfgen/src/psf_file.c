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
 *      $RCSfile: psf_file.c,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.16 $      $Date: 2019/07/30 21:53:57 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "psf_file.h"

#if defined(NEWPSFGEN)
#include "psf_file_extract.h"
#endif

#define PSF_RECORD_LENGTH 	160

/*
 * Read in the beginning of the bond/angle/dihed/etc information,
 * but don't read in the data itself.  Returns the number of the record type
 * for the molecule.  If error, returns (-1).
 */
int psf_start_block(FILE *file, const char *blockname) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int nrec = -1;

  /* keep reading the next line until a line with blockname appears */
  do {
    if(inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file)) {
      /* EOF encountered with no blockname line found ==> error, return (-1) */
      return (-1);
    }
    if(strlen(inbuf) > 0 && strstr(inbuf, blockname))
      nrec = atoi(inbuf);
  } while (nrec == -1);

  return nrec;
}


/* return # of atoms, or negative if error */
int psf_start_atoms(FILE *file) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int natom = 0;
  
  /* skip comments; get number of atoms */
  /* Taken from VMD's ReadPSF */
  do {
    if (inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, file)) {
      /* EOF with no NATOM */
      return -1;  
    }
    if (strlen(inbuf) > 0) {
      if (!strstr(inbuf, "REMARKS")) {
        if (strstr(inbuf, "NATOM")) {
          natom = atoi(inbuf);
        }
      }
    }
  } while (!natom);
  return natom;
}

#if !defined(NEWPSFGEN)
int psf_get_atom(FILE *f, char *name, char *atype, char *resname,
                 char *segname, char *resid, double *q, double *m) {

  char inbuf[PSF_RECORD_LENGTH+2];
  int i,num, read_count;

  if(inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, f)) {
    return(-1);
  }
  read_count = sscanf(inbuf, "%d %8s %8s %8s %8s %8s %lf %lf",
    &num, segname, resid, resname, name, atype, q, m);

  if (read_count != 8) {
    fprintf(stderr,"BAD ATOM LINE IN PSF FILE:\n: %s\n", inbuf);
    return -1;
  }
  if (sscanf(atype, "%d", &i) > 0) {
    fprintf(stderr, "PSF file is in CHARMM format; XPLOR format required.\n");
    return -1;
  }
  return num;
}

#else

int psf_get_atom(FILE *f, psfatom *atom, int drude) {

  char inbuf[PSF_RECORD_LENGTH+2];
  int i,num, read_count, tmpcount;

  if(inbuf != fgets(inbuf, PSF_RECORD_LENGTH+1, f)) {
    return(-1);
  }
  if (!drude) {
    read_count = sscanf(inbuf, "%d %8s %8s %8s %8s %8s %lf %lf %d",
      &num, atom->segname, atom->resid, atom->resname, atom->name, atom->atype, 
      &atom->charge, &atom->mass, &atom->lpd);
      tmpcount = 9;
  } else {
    read_count = sscanf(inbuf, "%d %8s %8s %8s %8s %8s %lf %lf %d %lf %lf",
      &num, atom->segname, atom->resid, atom->resname, atom->name, atom->atype, 
      &atom->charge, &atom->mass, &atom->lpd, &atom->alpha, &atom->thole);
      tmpcount = 11;
  }
  /* Initialize chain array*/
  atom->chain[0] = '\0';
  if (read_count != tmpcount) {
    fprintf(stderr,"BAD ATOM LINE IN PSF FILE:\n: %s\n", inbuf);
    return -1;
  }
  if (sscanf(atom->atype, "%d", &i) > 0) {
    fprintf(stderr, "PSF file is in CHARMM format; XPLOR format required.\n");
    return -1;
  }
  return num;
}

#endif

static int atoifw(char **ptr, int fw) {
  char *op = *ptr;
  int ival = 0; /* integer read from the file*/
  int iws = 0; /* size of the string */
  char tmpc;

  sscanf(op, "%d%n", &ival, &iws);
  if ( iws == fw ) { /* "12345678 123..." or " 1234567 123..." */
    *ptr += iws;
  } else if ( iws < fw ) { /* left justified? */
    while ( iws < fw && op[iws] == ' ' ) ++iws;
    *ptr += iws;
  } else if ( iws < 2*fw ) { /* " 12345678 123..." */
    *ptr += iws;
  } else { /* " 123456712345678" or "1234567812345678" */
    tmpc = op[fw];  op[fw] = '\0';
    ival = atoi(op);
    op[fw] = tmpc;
    *ptr += fw; /* go to the next set of character in the line (op)*/
  }
  return ival;
}


int psf_get_bonds(FILE *f, int fw, int n, int *bonds) {
  char inbuf[PSF_RECORD_LENGTH+2];
  char *bondptr = NULL;
  int i=0;
  while (i<n) {
    if((i % 4) == 0) {
      /* must read next line */
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        /* early EOF encountered */
        break;
      }
      bondptr = inbuf;
    }
    if((bonds[2*i] = atoifw(&bondptr,fw)) < 1)
      break;
    if((bonds[2*i+1] = atoifw(&bondptr,fw)) < 1)
      break;
    i++;
  }

  return (i != n);
}


int psf_get_angles(FILE *f, int fw, int n, int *angles) {
  char inbuf[PSF_RECORD_LENGTH+2];
  char *bondptr = NULL;
  int i=0;
  while (i<n) {
    if((i % 3) == 0) {
      /* must read next line */
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        /* early EOF encountered */
        break;
      }
      bondptr = inbuf;
    }
    if((angles[3*i] = atoifw(&bondptr,fw)) < 1)
      break;
    if((angles[3*i+1] = atoifw(&bondptr,fw)) < 1)
      break;
    if((angles[3*i+2] = atoifw(&bondptr,fw)) < 1)
      break;
    i++;
  }

  return (i != n);
}


int psf_get_dihedrals(FILE *f, int fw, int n, int *dihedrals) {
  char inbuf[PSF_RECORD_LENGTH+2];
  char *bondptr = NULL;
  int i=0;
  while (i<n) {
    if((i % 2) == 0) {
      /* must read next line */
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        /* early EOF encountered */
        break;
      }
      bondptr = inbuf;
    }
    if((dihedrals[4*i] = atoifw(&bondptr,fw)) < 1)
      break;
    if((dihedrals[4*i+1] = atoifw(&bondptr,fw)) < 1)
      break;
    if((dihedrals[4*i+2] = atoifw(&bondptr,fw)) < 1)
      break;
    if((dihedrals[4*i+3] = atoifw(&bondptr,fw)) < 1)
      break;
    i++;
  }

  return (i != n);
}


int psf_get_impropers(FILE *f, int fw, int n, int *impropers) {
  
  /* Same format */
  return psf_get_dihedrals(f, fw, n, impropers);
}


int psf_get_cmaps(FILE *f, int fw, int n, int *cmaps) {
  
  /* Same format */
  return psf_get_dihedrals(f, fw, 2*n, cmaps);
}


int psf_get_exclusions(FILE *f, int fw, int nexcl, int *exclusions, int natom, int *exclusion_indices) {
  char inbuf[PSF_RECORD_LENGTH+2];
  char *atomptr = NULL;
  // read list of excluded atoms
  int i=0;
  while (i<nexcl) {
    if((i % 8) == 0) {
      // must read next line
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        // early EOF encountered
        break;
      }
      atomptr = inbuf;
    }
    if((exclusions[i] = atoifw(&atomptr, fw)) < 1)
      break;
    i++;
  }
  if (i != nexcl)
    return 1;

  // read list of exclusion indices for each atom
  i=0;
  while (i<natom) {
    if((i % 8) == 0) {
      // must read next line
      if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
        // early EOF encountered
        break;
      }
      atomptr = inbuf;
    }
    if((exclusion_indices[i] = atoifw(&atomptr, fw)) < 0) {
      break;
    }
    i++;
  }
  return (i != natom);
}

#if defined(NEWPSFGEN)

int psf_get_lonepair_info(FILE *f, int fw, lonepair *lpair) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int read_count=0, tmpcount=0;

  if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
    return -1;
  }

 
  read_count = sscanf(inbuf, "%d %d %s %f %f %f",
    &lpair->numhost, &lpair->lpindex, lpair->type, 
    &lpair->distance, &lpair->angle, &lpair->dihedral);
  
  tmpcount = 6;

  /* the NUMLP section contain a character describing the type of lone pair (lpairs->type)
   * , but  so far this character is always "F", so no need to store this value.
   */
  if (read_count != tmpcount) {
    fprintf(stderr,"BAD LONE PAIR LINE IN PSF FILE:\n: %s\n", inbuf);
    return -2;
  }
  return 0;
}

/* Parse the index section of the lone pairs. */
int psf_get_lonepair_hosts(int fw, char **hostptr) {

  return atoifw(hostptr,fw);
}


int psf_get_aniso_tensors(FILE *f, int fw, psfaniso **aniso, int naniso) {
  char inbuf[PSF_RECORD_LENGTH+2];
  int i, read_count=0, tmpcount=0;
  psfaniso *tmp=NULL;

  i = 0;

  do {

    if(!fgets(inbuf,PSF_RECORD_LENGTH+2,f)) {
      return -1;
    }
    tmp = (psfaniso*)malloc(sizeof(psfaniso));

    read_count = sscanf(inbuf, "%f %f %f", &tmp->k11, &tmp->k22, &tmp->k33);
    
    tmpcount = 3;
    if (read_count != tmpcount) {
      return -1;
    }
    aniso[i] = tmp;

    i++;
  } while (i < naniso);
  
  return 0;
}

/* Parse the index section of the anisotropy. */
int psf_get_aniso_hosts(FILE *f, int fw, char **anisoptr, int *anisohost) {
  int i=0;

  while (i < 4) {

    if((anisohost[i] = atoifw(anisoptr,fw)) < 0) {
      return -1;
    }
    i++;
  }

  return 0;
}

#endif
