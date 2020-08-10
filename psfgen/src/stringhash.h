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
 *      $RCSfile: stringhash.h,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.3 $      $Date: 2019/07/24 03:36:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#ifndef STRINGHASH_H
#define STRINGHASH_H

struct stringhash;
typedef struct stringhash stringhash;

stringhash * stringhash_create(void);
void stringhash_destroy(stringhash *h);

const char* stringhash_insert(stringhash *h, const char *key, const char *data);

#define STRINGHASH_FAIL 0

const char* stringhash_lookup(stringhash *h, const char *key);

const char* stringhash_delete(stringhash *h, const char *key);

#endif

