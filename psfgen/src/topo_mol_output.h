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
 *      $RCSfile: topo_mol_output.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.7 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef TOPO_MOL_OUTPUT_H
#define TOPO_MOL_OUTPUT_H

#include <stdio.h>
#include "topo_mol.h"

int topo_mol_write_pdb(topo_mol *mol, FILE *file, void *, void *, 
                                void (*print_msg)(void *, void *, const char *));

int topo_mol_write_namdbin(topo_mol *mol, FILE *file, FILE *velfile, void *, void *, 
                                void (*print_msg)(void*, void *, const char *));

int topo_mol_write_psf(topo_mol *mol, FILE *file, int charmmfmt, int nocmap, 
                       int nopatches, void *, void *, 
                       void (*print_msg)(void *, void *, const char *));

#endif

