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
 *      $RCSfile: psfgen.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.9 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Defines set of data structures used in creation of molecule structures
 *   Exported here so that new modules can be written to interface with psfgen
 *  
 *   To compile the psfgen with new additions, like new atom data structure and
 *   faster dihedral detection, add the NEWPSFGEN flag to the compiler command
 *   The compiler flag NOIO is used to profile psfgen as it prevents the 
 *   writing the psf file
 *
 ***************************************************************************/

#ifndef PSFGEN_H
#define PSFGEN_H

#include "topo_defs.h"
#include "topo_mol.h"
#include "stringhash.h"

/* psfgen-specific data */
struct psfgen_data {
  int id, in_use, all_caps;
  topo_defs *defs;
  topo_mol *mol;
  stringhash *aliases;

  FILE *PSFGENLOGFILE;/* psfgen log file */
  int VPBONDS; /* flag to signal that the bonds between virtual particles 
                * (lone pairs and drude particles) and the their host must 
                * be explicitly printed to the psf
                */

};
typedef struct psfgen_data psfgen_data;

#if defined(NEWPSFGEN)
  


  #ifndef M_PI
  #define M_PI            3.14159265358979323846
  #endif

  #ifndef K_DRUDE
  #define K_DRUDE            500.0
  #endif

#endif

#endif /* PSFGEN_H */
