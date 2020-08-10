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
 *      $RCSfile: charmm_parse_topo_defs.h,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef CHARMM_PARSE_TOPO_DEFS_H
#define CHARMM_PARSE_TOPO_DEFS_H

#include <stdio.h>
#include "topo_defs.h"

int charmm_parse_topo_defs(topo_defs *defs, FILE *file, int all_caps,void *vdata, void *v,
                                void (*print_msg)(void *, void *,const char *));

#endif

