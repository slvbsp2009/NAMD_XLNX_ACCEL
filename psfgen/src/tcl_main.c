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
 *      $RCSfile: tcl_main.c,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.10 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#if defined(NAMD_TCL) || ! defined(NAMD_VERSION)

#include <tcl.h>

#if defined(NEWPSFGEN)
#include "psfgen.h"
#endif

extern int Psfgen_Init(Tcl_Interp *);

int main(int argc, char *argv[]) {
  Tcl_Main(argc, argv, Psfgen_Init);
  return 0;
}

#ifdef NAMD_VERSION
/* 
 * Provide user feedback and warnings beyond result values.
 * If we are running interactively, Tcl_Main will take care of echoing results
 * to the console.  If we run a script, we need to output the results
 * ourselves.
 */
#if defined(NEWPSFGEN) 
void newhandle_msg_text(psfgen_data *psfcontext, Tcl_Interp *interp, const char *msg);
#endif

void newhandle_msg(void *vdata, void *v, const char *msg) {
  Tcl_Interp *interp = (Tcl_Interp *)v;

  const char *words[3] = {"puts", "-nonewline", "psfgen) "};
  char *script = NULL;
  
#if defined(NEWPSFGEN)
  ClientData *data = (ClientData *)vdata;  
  // Get the psfge_data structure from data
  psfgen_data *psfcontext = *(psfgen_data **)data;
  /* If the log file was defined, redirect all messages there */
  if (psfcontext->PSFGENLOGFILE != NULL) {
    newhandle_msg_text(psfcontext, interp, msg);
    return;
  }
  
#endif
  // prepend "psfgen) " to all output 
  script = Tcl_Merge(3, words);
  Tcl_Eval(interp,script); 
  Tcl_Free(script);

  // emit the output
  words[1] = msg;
  script = Tcl_Merge(2, words);
  Tcl_Eval(interp,script);
  Tcl_Free(script);
}

/*
 * Same as above but allow user control over prepending of "psfgen) "
 * and newlines.
 */
void newhandle_msg_ex(void *vdata, void *v, const char *msg, int prepend, int newline) {
  Tcl_Interp *interp = (Tcl_Interp *)v;
  const char *words[3] = {"puts", "-nonewline", "psfgen) "};
  char *script = NULL;
  
#if defined(NEWPSFGEN)
  ClientData *data = (ClientData *)vdata; 
  // Get the psfge_data structure from data
  psfgen_data *psfcontext = *(psfgen_data **)data; 
  /* If the log file was defined, redirect all messages there */
  if (psfcontext->PSFGENLOGFILE != NULL) {
    newhandle_msg_text(psfcontext, (void *)v, msg);
    return;
  }
  
#endif
  
  if (prepend) { 
    // prepend "psfgen) " to all output
    script = Tcl_Merge(3, words);
    Tcl_Eval(interp,script);
    Tcl_Free(script);  
  } 
  
  // emit the output
  if (newline) {
    words[1] = msg;
    script = Tcl_Merge(2, words);
  } else {
    words[2] = msg;
    script = Tcl_Merge(3, words);
  }
  Tcl_Eval(interp,script);
  Tcl_Free(script);
} 
#endif

#else

#include <stdio.h>

int main(int argc, char **argv) {
  fprintf(stderr,"%s unavailable on this platform (no Tcl)\n",argv[0]);
  exit(-1);
}

#endif

