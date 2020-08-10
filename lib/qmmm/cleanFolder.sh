#!/bin/bash

# If more than one independent QM region exists in the system, but only one QM
# simulation per node is requested in the configuration file, all independent QM
# regions will be ran in the same folder (the idea behind limiting the number of QM
# regions being calculated per node is exactly saving RAM space). Some softwares,
# however, may leave behind temporary files which are QM region-specific, which 
# will lead to geometry errors when running multiple QM regions with one simulation
# per node.
#
# This script erases temporary files left behing by ORCA so that two *different* and 
# independent QM regions can be calculated in the same node (or single computer) 
# in *sequence*, and NOT in *parallel*. Without the script, temporary files left 
# behind after the first QM region was processed would be read by ORCA when 
# calculating the second QM region, but since they are different, ORCA would crash
# correctly claiming that the geometry in temporary files does not match the 
# geometry in the input file.
# 
# This will NOT be more efficient, as independent QM regions could be calculated in
# parallel, speeding up the overall time per simulation step, but the script serves 
# as an example of what could be done with a Secondary Process.

targDir=`dirname "$1"`

echo "Input file: $1"
echo "Target Dir: $targDir"

rm $targDir/*.gbw
rm $targDir/*.prop