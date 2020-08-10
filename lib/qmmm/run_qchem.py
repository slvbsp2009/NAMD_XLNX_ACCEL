#!/usr/bin/python

# by Maximilian Scheurer (mscheurer@ks.uiuc.edu), June 2017

##################  Imports

from sys import argv as sargv
from sys import exit
from os.path import dirname
import subprocess as sp
import re
import numpy as np
import math


################## Command file parameters
qchemNumberOfCores = 2

qchemConfigLine1 = """\
$rem
JOBTYPE force
METHOD b3lyp
BASIS 6-31G*
SCF_CONVERGENCE 7
$end

$molecule
0 1
"""





# if you specified an excited state ( CIS, ADC or whatever ) in the input file and want to give its gradient to NAMD,
# specify the number of the excited state here,
# so that the python script can read the corresponding excited state energy
# if no exc. state calculation is requested, set excitedState to 0
excitedState=0

# add the following lines to qchemConfigLine1 to calculate e.g. 4 excited states with TDDFT,
# where the gradient will be calculated on the 1st excited state
# CIS_N_ROOTS 4
# CIS_STATE_DERIV 1
# CIS_SINGLETS true
# CIS_TRIPLETS false



################## Processing of file names

inputFilename = sargv[1]

#print(inputFilename)

directory = dirname(inputFilename)

# Prepares the name of the configuration file based on full path to data file provided by NAMD
qchemInFileName = directory + "/"
qchemInFileName += "qmmm.in"

# Name of the qchem log file
qchemOutFileName = qchemInFileName +".out"

# Prepares the file name for the file which will be read by NAMD
finalResFileName = inputFilename
finalResFileName += ".result"

################## Reading and parsing NAMD's data ; Writing gaussian's Input

# Reads NAMD data
infile = open(inputFilename,"r")

line = infile.readline()

# Gets number of atoms in the Quantum Chemistry region (= QM atoms + Link atoms)
numQMatms = int(line.split()[0])
# Gets number of point charges
numPntChr = int(line.split()[1].replace("\n",""))

# print("numQMatms:",numQMatms,"; numPntChr",numPntChr)

# stores all lines written to gaussian's input file
outLinesQM = []

# Identation
ident = "  "

lineIndx = 1
pointChargeList = []
pointChargeDict = {}
charges = []
for line in infile:

    posx = line.split()[0]
    posy = line.split()[1]
    posz = line.split()[2]

    if lineIndx <= numQMatms:

        # qchem's format requires the fileds to be ordered begining with the
        # atom's element symbol, and followed by the XYZ coordinates.

        element = line.split()[3].replace("\n","")

        outLinesQM.append(" ".join([element,posx,posy,posz]) + "\n")

    else:
        #output linebreak to separate atoms from charges
        if lineIndx == numQMatms+1:
            outLinesQM.append("\n")

        # qchems's format requires the fields to be ordered begining with the
        # XYZ, and followed by the charge .
        pos = " ".join(line.split()[0:3])
        charge = line.split()[3]
        charges.append(charge)
        if pos in pointChargeDict:
        	print "double occurence: ", pos, pointChargeDict[pos] , " with new charge ", charge
        	pointChargeDict[pos] += float(charge)
        	print "new merged charge: ", pointChargeDict[pos]
        else:
        	pointChargeDict[pos] = float(charge)


    lineIndx += 1

cnp = np.array(charges,dtype=float)
print "Sum of the charges: ", np.sum(cnp)

outLinesQM.append("$end \n \n $external_charges \n")
for k in pointChargeDict:
	c = pointChargeDict[k]
	p = k.split()
	pointChargeList.append([p[0],p[1],p[2],str(c)])
	outLinesQM.append(" ".join([p[0],p[1],p[2],'{0:.16f}'.format(c)]) + "\n")

outLinesQM.append("$end")
# print len(pointChargeList)
infile.close()

###

with open(qchemInFileName,"w") as outQMFile:

    outQMFile.write(qchemConfigLine1)

    for line in outLinesQM:
        outQMFile.write(line)

    # outQMFile.write(gaussianWhitespace)

################## Run qchem

# We first move the shell to the target directory where calculations are to be
# performed
cmdline = "cd " + directory + "; "
# Then we run qchem with our output file receiving all standard output.

# we probably need to set some environment variables:
import subprocess, os
current_env = os.environ.copy()
current_env["QCSCRATCH"] = directory
# subprocess.Popen(my_command, env=my_env)

#load qchem from modules
cmdline += "qchem -nt %d " % qchemNumberOfCores
cmdline += qchemInFileName + " " + qchemOutFileName

print "command:", cmdline
proc = sp.Popen(args=cmdline, shell=True, env=current_env)
proc.wait()

########### READING QCHEM output
excitedStateEnergy=0

tmpOutFile = open(qchemOutFileName,"r")

qmCharges = []
gradient = []

# Bohr radius for conversion
a0 = 0.52917721067
conversionHB = 627.509469/a0
selfEnergy = 0.0

# Iterates ultil we find the section of output that contains atomic partial
# charges for QM atoms
chargeSection = False

gradientSection = False

scfenergyFound = False

excitedFound = False

iterate = True

selfEnergyFound = False

while iterate:

    line = tmpOutFile.readline()

    if line.find("Ground-State Mulliken Net Atomic Charges") != -1:
        chargeSection = True
        # Skips 3 dividing lines
        for x in range(3):
            line = tmpOutFile.readline()
        continue

    if not scfenergyFound:
        if line.find("SCF   energy in the final basis set") != -1:
            scfenergy = 627.509469 * float(line.strip().split()[-1])
            print "SCF energy: ", scfenergy
            scfenergyFound = True

    if not excitedFound and excitedState != 0:
        if line.find("Excited State") != -1:
            line = line.strip().replace(":","").split()
            if int(line[2]) != excitedState:
                continue
            line =tmpOutFile.readline()
            # check if we really requested the excited state that the gradient will be calculated for!
            while line.find("This state for optimization and/or second-order correction.") == -1:
                line = tmpOutFile.readline()
            line = tmpOutFile.readline()
            line = line.strip()
            reg = re.search("([-+]?\d*\.\d+|\d+)(.+)",line)
            excitedStateEnergy = 627.509469 * float(reg.group(1))
            print "Excited State energy: " , excitedStateEnergy
            excitedFound = True

    if line.find("Gradient of") != -1:
        gradientSection = True
        # line = tmpOutFile.readline()


    if gradientSection:
        gradient = np.zeros(shape=(numQMatms,3))

        ln = 0
        blockColumnNumber = -1
        atomLineIndex = 0
        gradientRead = 0
        while (gradientRead == 0):
            line = tmpOutFile.readline()
            line = line.strip().split()
            if ln % 4 == 0:
                blockColumnNumber = len(line)
                # print(line, blockColumnNumber)
                atomLineIndex += blockColumnNumber
            else:
                xyzIndex = (ln % 4)-1
                # print("current coordinate: "+ str(xyzIndex))
                currentColumnNumber = blockColumnNumber
                columnCounter = 1
                for atomIndex in range(atomLineIndex-currentColumnNumber, atomLineIndex):
                    gradient[atomIndex,xyzIndex] = float(line[columnCounter])
                    #print(atomIndex,columnCounter)
                    if atomIndex == (numQMatms-1) and xyzIndex == 2:
                        gradientRead = 1
                    columnCounter += 1
            ln+=1

        iterate = 0

    if chargeSection:
        qmCharges.append(line.split()[2].replace("\n","").strip())

    if len(qmCharges) == numQMatms:
        chargeSection = False

    if selfEnergyFound != True and line.find("Charge-charge energy") != -1:
	line = line.strip().split()
    	selfEnergy = float(line[-2])
    	print "Self energy of the charges: ", selfEnergy
        selfEnergyFound = True

tmpOutFile.close()

finFile = open(finalResFileName,"w")

if excitedState == 0:
	finFile.write(str( scfenergy - selfEnergy*627.509469 ) + "\n")
	print "Corrected SCF energy: ", scfenergy-selfEnergy*627.509469
else:
	finFile.write(str( excitedStateEnergy - selfEnergy*627.509469 ) + "\n")

forces=np.multiply(-1.0*conversionHB,gradient)
forces=forces.tolist()
#print forces
for i in range(numQMatms):
    finFile.write(' '.join('{0:.7f}'.format(c) for c in forces[i]) + " " + qmCharges[i] + "\n")

finFile.close()


##########

exit(0)
