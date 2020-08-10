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
 *      $RCSfile: psf_file_extract.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.10 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef PSF_FILE_READ_H 
#define PSF_FILE_READ_H 

#include <stdio.h>
#include "topo_mol.h"

#if defined(NEWPSFGEN)

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

/* Read in all psf atom information using this struct */
typedef struct psfatom {
  char name[10];
  char atype[10];
  char resname[10];
  char segname[10];
  char chain[10];
  char resid[10];
  char element[3];
  double charge, mass;
  /* data present in the atoms section that can be used to identify the atoms
   * as lone pairs or drudes, alpha and thole are also stored in the atoms 
   * sections in the drude force field (column 10 and 11)
   */
  int lpd; /* this int identifies the atom being read as lone pair (-1)
            * or drude (-2) by reading the 9 column of the psf
            */
  double alpha, thole;

} psfatom;


/* Read in lone pair information using this struct */
typedef struct lonepair {
  char type[10];
  int numhost, lpindex;
  float distance;
  float angle; /* scale in case of the colinear*/
  float dihedral; /* ignored in case of the colinear*/

} lonepair;


/* Read in Anisotropy information using this struct */
typedef struct psfaniso {
  float k11; /* definition of the anisotropic tensor 11 22 defined in the RTF
              * and the 33 is dedeucted from the previous 2
              */
  float k22;
  float k33;

} psfaniso;

#endif

int psf_file_extract(topo_mol *mol, FILE *file, FILE *pdbfile, FILE *namdbinfile, FILE *velnamdbinfile,
                                void *, void *, void (*print_msg)(void *, void *, const char *));

#endif

