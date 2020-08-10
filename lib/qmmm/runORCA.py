#!/usr/bin/python3

# by Marcelo Melo (melomcr@gmail.com)

# If you would like to run a Quantum Chemistry software which is not supported by NAMD,
# an interface script can be created in order to convert NAMD's input into the format
# expected by the QC software, and subsequently convert the software's output inro the
# format expected by NAMD.
# 
# This script creates an interface for ORCA, so that it can be ran using the custom 
# QM software option provided by NAMD. It is an example script, since ORCA is supported
# natively by NAMD's QM/MM interface, and is intended to be used by those wishing to 
# test theyr own software in hybrid MD simulations through NAMD.

#######  NAMD Files Format
# INPUT: This input file will contain on the first line the number of QM atoms (X)
# and the number of point charges in the file (Y, which may be 0 or more), separated 
# by a space. The following X+Y lines will have four (4) fields: X, Y and Z coordinates,
# and a fourth field which will depend on the type of entry. For QM atoms, the 
# field will contain the element of the QM atom. For point charge lines, the field
# will contain the charge of the point charge.
# 
# OUTPUT: The expected output file whould be placed in the same directory as the 
# input file, and should be named "*inputfile*.result" (meaning it will have the
# same path and name of the input file, plus the suffix '.result'). This file
# should have, on its first line, the energy of the system, and on the following 
# X lines (where X is the number of QM atoms passed in the input file), four (4) 
# fields: the x, y and z components of the TOTAL FORCE applied on that atom, and 
# on the fourth field, the charge of the atom. If the user indicates that charges
# from the QM software should not be used (see "QMChargeMode"), the fourth field
# should have zeroes, but should not be empty.


#######  ORCA Files Format

# ORCA input file is composed of a variable-sized header, where all info pertaining
# to the simulation is specified, and a corrdinates section, where all atom positions
# are specified, along with some informations regarding the overall system (namely
# the charge and multiplicity of the system).
#
# We will get the header and the initial part of the coordinates section from 
# lines specified in this script. The data would be gathered and writen by NAMD,
# but in this case we need to previously prepare it from hardcoded values or by 
# reading from another file.
# 
# ORCA expects a separate file with point charge information. This info is provided
# by NAMD in the same data file, so we will extract it and paste it in a separate 
# file in the same folder. The name and path to the point charge file needs to be 
# included in the configuration portion of ORCA's input file.

##################  Imports

from sys import argv as sargv
from sys import exit 
from os.path import dirname
import subprocess as sp

################## Hardcoded parameters

# Here we set up a template for the header and options of the input file.

# The choice of method and basis-set, as well as parallelization and charge
# output are all made here.
# For ORCA, the keyword "ENGRAD" is essential for the proper output of gradients
# over QM atoms.

orcaConfigLines1 = """\
!  B3LYP 6-31G* Grid4 PAL4 EnGrad TightSCF
%output 
  PrintLevel Mini 
  Print [ P_Mulliken ] 1
  Print [ P_AtCharges_M ] 1
end

"""

# We set up the template for the point charge file indicator

orcaConfigLines2 = "%pointcharges \""

# And set up the header for the coordinates section of the input file.

orcaConfigLines3 = """\
%coords
  CTyp xyz
  Charge 0.000000
  Mult 1.000000
  Units Angs
  coords

"""

################## Processing of file names

inputFilename = sargv[1]

#print(inputFilename)

directory = dirname(inputFilename)

# Prepares the name of the configuration file based on full path to data file provided by NAMD
orcaInFileName = directory + "/"
orcaInFileName += "qmmm.input"

# Prepares the name of the point charge file based on full path to data file provided by NAMD
#%pointcharges "/dev/shm/NAMD_test/0/qmmm_0.input.pntchrg"
pcFileName = orcaInFileName
pcFileName += ".pntchrg"

orcaConfigLines2 += pcFileName + "\"\n"


# Prepares the name of the output file based on full path to data file provided by NAMD.
# This is where we will direct all ORCA output, so that we can grab the atom charge
# information to pass back to NAMD.
orcaOutFileName = orcaInFileName
orcaOutFileName += ".TmpOut"

# Prepares the name of the gradient file based on full path to data file provided by NAMD.
# This is where ORCA will write the gradient information on QM atoms.
orcaGradFileName = orcaInFileName
orcaGradFileName += ".engrad"

# Prepares the file name for the file which will be read by NAMD
finalResFileName = inputFilename
finalResFileName += ".result"

#print("orcaInFileName:",orcaInFileName)
#print("pcFileName:",pcFileName)
#print("orcaOutFileName:",orcaOutFileName)
#print("orcaGradFileName:",orcaGradFileName)
#print("finalResFileName:",finalResFileName)

################## Reading and parsing NAMD's data ; Writing ORCA's Input

# Reads NAMD data
infile = open(inputFilename,"r") 

line = infile.readline()

# Gets number of atoms in the Quantum Chemistry region (= QM atoms + Link atoms)
numQMatms = int(line.split()[0])
# Gets number of point charges
numPntChr = int(line.split()[1].replace("\n",""))

#print("numQMatms:",numQMatms,"; numPntChr",numPntChr)

# stores all lines written to ORCA's input file
outLinesQM = []

# stores all lines written to ORCA's point charge file
outLinesPC = []

# The first line in the point charge file is composed only of the total number
# of point charges in the file.
outLinesPC.append(str(numPntChr) + "\n")

# Identation
ident = "  "

lineIndx = 1
for line in infile:
    
    posx = line.split()[0]
    posy = line.split()[1]
    posz = line.split()[2]
    
    if lineIndx <= numQMatms:
        
        # ORCA's format requires the fileds to be ordered begining with the
        # atom's element symbol, and followed by the XYZ coordinates.
                
        element = line.split()[3].replace("\n","")
        
        outLinesQM.append(ident + " ".join([element,posx,posy,posz]) + "\n")
        
    else:
        
        # ORCA's format requires the fileds to be ordered begining with the
        # charge, and followed by the XYZ coordinates.
        
        charge = line.split()[3]
        
        outLinesPC.append(" ".join([charge,posx,posy,posz]) + "\n")
    
    lineIndx += 1
    

# Finalizes the formating for ORCA's input file, one "end" to terminate the
# atomic coordinates, and nother to terminate the section.
outLinesQM.append(ident + "end" + "\n")
outLinesQM.append("end" + "\n")


infile.close()

###

with open(orcaInFileName,"w") as outQMFile:
    
    outQMFile.write(orcaConfigLines1)
    outQMFile.write(orcaConfigLines2)
    outQMFile.write(orcaConfigLines3)
    
    for line in outLinesQM:
        outQMFile.write(line)

with open(pcFileName,"w") as outPCFile:
    
    for line in outLinesPC:
        outPCFile.write(line)



################## Run ORCA

# We first move the shell to the target directory where calculations are to be
# performed
cmdline = "cd " + directory + "; "
# Then we run orca with our output file receiving all standard output.
cmdline += "/data/Programas/orca_3_0_3_linux_x86-64/orca "
cmdline += orcaInFileName + " > " + orcaOutFileName

#print("command:",cmdline)

proc = sp.Popen(args=cmdline, shell=True)
proc.wait()


# Runs a secondary process
cmdline2 = "/home/melomcr/Research/NAMD_QMMM/Scripts/saveDipole.sh "
cmdline2 += orcaInFileName

proc2 = sp.Popen(args=cmdline2, shell=True)

################## Process ORCA results

# Gets system energy and gradients from ORCA's output file.

gradFile = open(orcaGradFileName,"r")
    
# Skips to the line with number of atoms
for i in range(3):
    gradFile.readline()

orcaNumQMAtms = int(gradFile.readline().replace("\n",""))

# Runs basic sanity check
if orcaNumQMAtms != numQMatms:
    print("ERROR: Expected",numQMatms,"but found",orcaNumQMAtms,"atoms in engrad file!")
    exit(1)

# Skips to the line with final system energy
for i in range(3):
    gradFile.readline()

finalEnergy = gradFile.readline().replace("\n","").strip()

print("ORCA energy: ", finalEnergy,"Eh")

# All energies are given in Eh (Hartree)
# NAMD needs energies in kcal/mol
# The conversion factor is 627.509469
finalEnergy = str( float(finalEnergy) * 627.509469 )

print("ORCA energy: ", finalEnergy,"kcal/mol")

# Skips to the lines with gradients
for i in range(3):
    gradFile.readline()

# Stores all gradients (yes, it would be faster with numpy, but this is just an example)
grads = []
for i in range(orcaNumQMAtms):
    
    grads.append( list() )
    
    # All *GRADIENTS* are given in Eh/a0 (Hartree over Bohr radius)
    # NAMD needs *FORCES* in kcal/mol/angstrons
    # The conversion factor is -1*627.509469/0.529177 = -1185.82151
    for j in range(3):
        gradComp = gradFile.readline().replace("\n","").strip()
        gradComp = float(gradComp) * -1185.82151
        grads[i].append( str(gradComp) )

gradFile.close()


###

# Gets atom charges from the temporary output file.
tmpOutFile = open(orcaOutFileName,"r")

qmCharges = []

# Iterates ultil we find the section of output that contains atomic partial 
# charges for QM atoms
chargeSection = False

iterate = True

while iterate:
    
    line = tmpOutFile.readline()
    
    if line.find("MULLIKEN ATOMIC CHARGES") != -1:
        chargeSection = True
        
        # Skips a dividing line
        line = tmpOutFile.readline()
        
        continue
    
    if chargeSection:
        qmCharges.append(line.split()[3].replace("\n","").strip())
        pass
    
    if len(qmCharges) == numQMatms:
        break

tmpOutFile.close()

# Writes the final output file that will be read by NAMD.

finFile = open(finalResFileName,"w")

# The first line constains the enegy
finFile.write(finalEnergy + "\n")

# And all follwoing lines contain gradients and partial charges
for i in range(numQMatms):
    
    finFile.write(" ".join(grads[i]) + " " + qmCharges[i] + "\n")

finFile.close()


##########

# Makes sure the secondary process has ended befire we return to NAMD.
proc2.wait()


exit(0)
