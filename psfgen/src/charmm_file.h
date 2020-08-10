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
 *      $RCSfile: charmm_file.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $      $Date: 2019/07/24 03:36:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef CHARMM_FILE_H
#define CHARMM_FILE_H

#include <stdio.h>

int charmm_get_tokens(char **tok, int toklen,
			char *sbuf, int sbuflen,
			char *lbuf, int *lineno,
			FILE *stream, int all_caps);

#endif

