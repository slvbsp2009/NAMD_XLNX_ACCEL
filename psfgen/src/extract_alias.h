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
 *      $RCSfile: extract_alias.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $      $Date: 2019/07/24 03:36:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef EXTRACT_ALIAS_H
#define EXTRACT_ALIAS_H

#include "stringhash.h"

#define EXTRACT_ALIAS_FAIL -1

int extract_alias_residue_define(stringhash *h,
			const char *altres, const char *realres);

int extract_alias_atom_define(stringhash *h, const char *resname,
			const char *altatom, const char *realatom);

const char * extract_alias_residue_check(stringhash *h,
						const char *resname);

const char * extract_alias_atom_check(stringhash *h,
			const char *resname, const char *atomname);

#endif

