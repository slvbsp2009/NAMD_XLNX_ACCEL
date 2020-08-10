#!/usr/bin/python

# by Maximilian Scheurer (mscheurer@ks.uiuc.edu), June 2017

##################  Imports

from sys import argv as sargv
from sys import exit
from os.path import dirname
import subprocess as sp
import numpy as np
from scipy import spatial

# TeraChem does not print the point-charge self energy, so we need to
# calculate it ourselves. This may take a lot of time if many point charges
# are included due to O(N^2) scaling!
def computePointChargeSelfEnergy(pclist):
      energy = 0.0
      length = len(pclist)
      pc = np.array(pclist,dtype=np.float64)
      a = np.arange(0,length,1)
      coords = pc[:,:3]
      r = spatial.distance.cdist(coords,coords)/0.52917721067
      idc = np.triu_indices(length,1)
      energy = np.sum(pc[:,3][idc[0]]*pc[:,3][idc[1]]/(r[idc]))
      return energy



################## Hardcoded parameters

# Here we set up a template for the header and options of the input file.

excitedState=0

tcConfigLines1 = """\
basis 6-31G*
coordinates qmmm.xyz
pointcharges point_charges
charge 0
spinmult 1
poptype mulliken
method b3lyp
dftd no
run gradient
end
"""

inputFilename = sargv[1]

directory = dirname(inputFilename)

# Prepares the name of the configuration file based on full path to data file provided by NAMD
tcInFileName = directory + "/"
tcInFileName += "qmmm_tc.in"

coordFile = directory + "/qmmm.xyz"

# Prepares the name of the point charge file based on full path to data file provided by NAMD
#%pointcharges "/dev/shm/NAMD_test/0/qmmm_0.input.pntchrg"
pcFileName = directory + "/point_charges"


# Prepares the name of the output file based on full path to data file provided by NAMD.
# This is where we will direct all tc output, so that we can grab the atom charge
# information to pass back to NAMD.
tcOutFileName = tcInFileName
tcOutFileName += ".out"

# Prepares the file name for the file which will be read by NAMD
finalResFileName = inputFilename
finalResFileName += ".result"

# print("tcInFileName:",tcInFileName)
# print("pcFileName:",pcFileName)
# print("tcOutFileName:",tcOutFileName)
#print("tcGradFileName:",tcGradFileName)
# print("finalResFileName:",finalResFileName)

################## Reading and parsing NAMD's data ; Writing tc's Input
# Reads NAMD data
infile = open(inputFilename,"r")

line = infile.readline()

# Gets number of atoms in the Quantum Chemistry region (= QM atoms + Link atoms)
numQMatms = int(line.split()[0])
# Gets number of point charges
numPntChr = int(line.split()[1].replace("\n",""))

#print("numQMatms:",numQMatms,"; numPntChr",numPntChr)

# stores all lines written to tc's input file
outLinesQM = []

# stores all lines written to tc's point charge file
outLinesPC = []

# The first line in the point charge file is composed only of the total number
# of point charges in the file.

# Identation
ident = "  "

lineIndx = 1
lineIndx = 1
pointChargeList = []
pointChargeDict = {}
charges = []
for line in infile:

    posx = line.split()[0]
    posy = line.split()[1]
    posz = line.split()[2]

    if lineIndx <= numQMatms:


        element = line.split()[3].replace("\n","").capitalize()

        outLinesQM.append(ident + " ".join([element,posx,posy,posz]) + "\n")

    else:


        charge = line.split()[3]
        pos = " ".join(line.split()[0:3])
        charges.append(charge)
        if pos in pointChargeDict:
                print("double occurence: ", pos, pointChargeDict[pos] , " with new charge ", charge)
                pointChargeDict[pos] += float(charge)
                print("new merged charge: ", pointChargeDict[pos])
        else:
                pointChargeDict[pos] = float(charge)


    lineIndx += 1


cnp = np.array(charges,dtype=float)
print("Sum of the charges: ", np.sum(cnp))

outLinesPC.append(str(len(pointChargeDict)) + "\n\n")
for k in pointChargeDict:
        c = pointChargeDict[k]
        p = k.split()
        pointChargeList.append([p[0],p[1],p[2],str(c)])
        outLinesPC.append(" ".join(['{0:.4f}'.format(c),p[0],p[1],p[2]]) + "\n")

pcSelfEnergy = computePointChargeSelfEnergy(pointChargeList)
print("self-energy: " , str(pcSelfEnergy) , " a.u.")

infile.close()

###

with open(tcInFileName,"w") as outQMFile:

    outQMFile.write(tcConfigLines1)


with open(pcFileName,"w") as outPCFile:

    for line in outLinesPC:
        outPCFile.write(line)

with open(coordFile,"w") as cFile:
	cFile.write(str(numQMatms))
	cFile.write("\n\n")
	for line in outLinesQM:
		cFile.write(line)

################## Run tc

# We first move the shell to the target directory where calculations are to be
# performed
import os
current_env = os.environ.copy()

cmdline = "cd " + directory + "; "
# Then we run tc with our output file receiving all standard output.

# here, scipy is not regularly in the python path :(
cmdline += "$TeraChem/bin/terachem "
cmdline += tcInFileName + " > " + tcOutFileName

#print("command:",cmdline)

proc = sp.Popen(args=cmdline, shell=True,env=current_env)
proc.wait()


print("terachem finished running.")
################## Process tc results

scfEnergy = 0
iterate = True

gradientSection = False
excSection = False
grads = []
a0 = 0.52917721067
conversionHB = 627.509469/a0
eEnergy = 0
resultFile = open("scr/results.dat","r")
mainFile = open(tcOutFileName,"r")

if not excitedState:
	while iterate:
		line = mainFile.readline()

		if line.find("dE/dX") != -1:
			gradientSection = True
			line = mainFile.readline()

		if gradientSection == True:
			g = line.strip().split()[0:3]
			grads.append([-1.0*float(x)*conversionHB for x in g])
			if len(grads) == numQMatms:
				gradientSection = False
				break


iterate = True
while iterate:
	line = resultFile.readline()
	if line == None:
		break
	if line.find("Ground state energy (a.u.):") != -1:
		line = resultFile.readline()
		scfEnergy = (float(line.strip()) - pcSelfEnergy)*627.509469
		if excitedState == 0:
			break

	if line.find("dE/dX") != -1:
		gradientSection = True
		line = resultFile.readline()

	if gradientSection == True:
		g = line.strip().split()[1:4]
		grads.append([-1.0*float(x)*conversionHB for x in g])
		if len(grads) == numQMatms:
			gradientSection = False
			if excitedState == 0:
				break

	if line.find("Root	Energy") != -1:
		line = resultFile.readline()
		excSection = True

	if excSection == True:
		state, energy = line.strip().split()
		if int(state) == excitedState:
			# eV to kcal/mol
			eEnergy = float(energy) * 23.06054819
			break


resultFile.close()
finalEnergy = scfEnergy + eEnergy
print(scfEnergy,eEnergy,finalEnergy)
###

# Gets atom charges from the charge_mull.xls output file.
mullikenFile = open("scr/charge_mull.xls","r")

qmCharges = []

iterate = True

while iterate:
    line = mullikenFile.readline()
    qmCharges.append(line.split()[2].replace("\n","").strip())

    if len(qmCharges) == numQMatms:
        break

mullikenFile.close()
# Writes the final output file that will be read by NAMD.
finFile = open(finalResFileName,"w")

# The first line constains the enegy
finFile.write(str(finalEnergy) + "\n")

# And all follwoing lines contain gradients and partial charges
for i in range(numQMatms):

    finFile.write(" ".join(str(x) for x in grads[i]) + " " + qmCharges[i] + "\n")

finFile.close()


##########

exit(0)
