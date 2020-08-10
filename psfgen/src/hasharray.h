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
 *      $RCSfile: hasharray.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $      $Date: 2019/07/24 03:36:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef HASHARRAY_H
#define HASHARRAY_H

struct hasharray;
typedef struct hasharray hasharray;

hasharray * hasharray_create(void **itemarray, int itemsize);
int hasharray_clear(hasharray *a);
void hasharray_destroy(hasharray *a);

int hasharray_reinsert(hasharray *a, const char *key, int pos);
int hasharray_insert(hasharray *a, const char *key);

int hasharray_delete(hasharray *a, const char *key);

#define HASHARRAY_FAIL -1

int hasharray_index(hasharray *a, const char *key);

int hasharray_count(hasharray *a);

#endif

