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
 *      $RCSfile: memarena.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $      $Date: 2019/07/24 03:36:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef MEMARENA_H
#define MEMARENA_H

struct memarena;
typedef struct memarena memarena;

memarena * memarena_create(void);
void memarena_destroy(memarena *a);

void memarena_blocksize(memarena *a, int blocksize);
void * memarena_alloc(memarena *a, int size);
void * memarena_alloc_aligned(memarena *a, int size, int alignment);

#endif

