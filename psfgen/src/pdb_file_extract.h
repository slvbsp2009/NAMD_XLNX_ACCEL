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
 *      $RCSfile: pdb_file_extract.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.5 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef PDB_FILE_EXTRACT_H
#define PDB_FILE_EXTRACT_H

#include <stdio.h>
#include "stringhash.h"
#include "topo_mol.h"

int pdb_file_extract_residues(topo_mol *mol, FILE *file, stringhash *h, int all_caps,
                                void *, void *, void (*print_msg)(void *, void *,const char *));

int pdb_file_extract_coordinates(topo_mol *mol, FILE *file, FILE *namdbinfile,
                                const char *segid, stringhash *h, int all_caps,
                                void *, void *,void (*print_msg)(void *, void *,const char *));

#endif

