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
 *      $RCSfile: psf_file.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $      $Date: 2019/08/01 18:48:37 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef PSF_FILE_H 
#define PSF_FILE_H 

#include <stdio.h>

#if defined(NEWPSFGEN)
#include "psf_file_extract.h"
#endif

int psf_start_atoms(FILE *);
int psf_start_block(FILE *, const char *blockstr);

#if !defined(NEWPSFGEN)
int psf_get_atom(FILE *f, char *name, char *atype, char *resname,
                 char *segname, char *resid, double *q, double *m);
#else                 
int psf_get_atom(FILE *f, psfatom *atom, int drude);
#endif

int psf_get_bonds(FILE *f, int fw, int n, int *bonds);
int psf_get_angles(FILE *f, int fw, int n, int *angles);
int psf_get_dihedrals(FILE *f, int fw, int n, int *dihedrals);
int psf_get_impropers(FILE *f, int fw, int n, int *impropers);
int psf_get_cmaps(FILE *f, int fw, int n, int *cmaps);
int psf_get_exclusions(FILE *f, int fw, int nexcl, int *exclusions, int natom, int *exclusion_indices);


#if defined(NEWPSFGEN)
int psf_get_lonepair_info(FILE *f, int fw, lonepair *lpair);
int psf_get_lonepair_hosts(int fw, char **hostptr);
int psf_get_aniso_tensors(FILE *f, int fw, psfaniso **aniso, int naniso);
int psf_get_aniso_hosts(FILE *f, int fw, char **anisoptr, int *anisohost);
#endif

#endif

