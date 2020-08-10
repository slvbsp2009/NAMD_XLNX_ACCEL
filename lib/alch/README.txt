####################################################
#
# Scripts for Alchemical transformations in NAMD
#
####################################################


This directory contains scripts that are helpful for either running alchemical
transformations, or analyzing their results.


# fep.tcl
Driver script to be sourced in a NAMD config file. Provides a concise syntax
to run alchemical transformations in sequential windows.


# deinterleave_idws.py
A tool for post-processing fepout files obtained with NAMD's
Interleaved Double-Wide Sampling (option alchLambdaIDWS)

deinterleave_idws.py depends on SciPy.interpolate; it takes a fepout file containing
interleaved backward and forward energy differences, and splits the data it into distinct
fepout files that can be processed by other tools e.g. the ParseFEP plugin of VMD.
