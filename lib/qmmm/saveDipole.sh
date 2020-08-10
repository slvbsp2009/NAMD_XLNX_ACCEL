#!/bin/bash

targDir=`dirname "$1"`
ID=`basename "$targDir"`

echo "Input file: $1"
echo "Target Dir: $targDir"
echo "Target ID: $ID"

dp=`cat $1.TmpOut | grep "Total Dipole Moment" | sed -e 's/ \{1,\}/,/g' -e 's/^.*:,\(.*\)/\1/'`
echo "Dipole moment: $dp"

if [ ! -f "/home/melomcr/Research/NAMD_QMMM/TestSystem/dipoleMoments_$ID.csv" ] ; then
echo "x,y,z" > /home/melomcr/Research/NAMD_QMMM/TestSystem/dipoleMoments_$ID.csv
fi

echo "$dp" >> /home/melomcr/Research/NAMD_QMMM/TestSystem/dipoleMoments_$ID.csv