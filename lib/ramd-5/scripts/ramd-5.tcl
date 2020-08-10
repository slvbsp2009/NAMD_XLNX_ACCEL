#***********************************************************************************
#                                                                                  *
# Random Acceleration Molecular Dynamics (RAMD)                                    *
# Implementation for NAMD v2.13                                                    *
# Version 5.0.4 July 2018                                                     *
#                                                                                  *
# Copyright (c) 2018, HITS gGmbH, Heidelberg, Germany                              * 
# Authors: Vlad Cojocaru, Stefan Richter, Daria Kokh, Rebecca Wade                 *
# Email: mcmsoft@h-its.org                                                         *
#                                                                                  *
# The first Tcl script to run RAMD in NAMD (v2.5+) was written by Harish Vashisth  *
# Ref: Vashisth H et al, Biophys J. 2008 Nov 1;95(9):4193-204. Epub 2008 Aug 1     *
#                                                                                  *
# The Vashisth script inspired some lines of this script, but mainly this script   *
#    is based on the fortran code written for AMBER 8                              *
#                                                                                  *
# The structure of this script is inspired by the Adaptive Biasing Force module    *
#    distributed with NAMD 2.6                                                     *
#                                                                                  *
# The original RAMD method is described in:                                        *
#    Ref1: Luedemann,S.K.,Lounnas,V.and R.C.Wade.,                                 *
#          J Mol Biol, 303:797-811 (2000)                                          *
#    Ref2: Schleinkofer,K.,Sudarko,Winn,P.,Luedemann,S.K.and R.C.Wade,             *
#          EMBO Reports, 6, 584-589 (2005)                                         *
#                                                                                  *
#                                                                                  *
# Disclaimer:  This script is for research purposes only. EML Research does not    *
#              assume any responsibility for the software or its use.              * 
#                                                                                  *
#   The script along with usage examples is available at                           *
#   https://www.h-its.org/en/research/mcm/software/                                *
#***********************************************************************************

#*******************************************************
# Startup                                              *
#*******************************************************
package provide ramd 5.0

#*******************************************************
# Parameter definitions
#*******************************************************

namespace eval ::RAMD {
 set version "5.0"
 if {! [info exists RAMDdir]} { set RAMDdir [file dirname [info script]] }
 # If it fails, try the local directory
 if { $RAMDdir == "" } { set RAMDdir "." }
 
 TclForces		on
 TclForcesScript	$RAMDdir/ramd-5_script.tcl
 array set defaults {
  ramdSteps               50
  namdVersion             2.13
  forceRAMD                16.0
  rMinRamd                 0.01
  forceOutFreq             10
  firstProtAtom            1
  ramdSeed             14253
  mdSteps                  0
  mdStart                 no
  maxDist                 50
  debugLevel               0
  ramdfilename         "ramd.log"
 }

 set mandatory "firstRamdAtom lastRamdAtom lastProtAtom"
 set silent "rMinMd"

 array set capitals {}
 foreach param [concat $mandatory $silent [array names defaults]] {
  set capitals([string tolower $param]) $param
  # not set yet
  set alreadySet($param) 0
 }
} ;# namespace


proc ramd { keyword value } {
 set ::RAMD::keyword $keyword
 set ::RAMD::value $value

 namespace eval ::RAMD {
  # Build list of all allowed parameter names
  set list [array names capitals]
  set lowercase [string tolower $keyword]
  # Process parameters
  if {[lsearch $list $lowercase] != -1} {
   set keyword $capitals($lowercase)
   if { $alreadySet($keyword) } {
    print "RAMD> WARNING - multiple definitions of parameter $keyword"	
   }
   set $keyword $value
   set alreadySet($keyword) 1
   return
  } else {
   error [format "Unknown RAMD keyword: %s" $keyword]
  }
 } ;# namespace

} ;# proc ramd

# define upper-case synonyms to proc abf 
proc RAMD { keyword value } {
    ramd $keyword $value
}
proc Ramd { keyword value } {
    ramd $keyword $value
}

