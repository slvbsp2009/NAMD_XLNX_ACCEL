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
 *      $RCSfile: topo_mol_pluginio.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.10 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef TOPO_MOL_PLUGINIO_H
#define TOPO_MOL_PLUGINIO_H

#include <stdio.h>
#include "topo_mol.h"
#include "stringhash.h"

int topo_mol_read_plugin(topo_mol *mol, const char *pluginname,
                         const char *filename, 
                         const char *coorpluginname, const char *coorfilename,
                         const char *segid, stringhash *h, int all_caps,
                         int coordinatesonly, int residuesonly,
                         void *, void *, void (*print_msg)(void*, void *, const char *));

struct image_spec {
  int na, nb, nc;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
};

int topo_mol_write_plugin(topo_mol *mol, const char *pluginname,
                          const char *filename, struct image_spec *images,
                          void *, void *, void (*print_msg)(void *, void *, const char *));

#endif

