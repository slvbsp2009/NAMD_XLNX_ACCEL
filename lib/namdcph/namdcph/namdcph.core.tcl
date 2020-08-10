# 
# Core utilities for setting up and running constant pH simulations in NAMD.
#
package require Tcl 8.5

source [file join [file dirname [info script]] "namdtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]
source [file join [file dirname [info script]] "cphtitrator.tcl"]
source [file join [file dirname [info script]] "json.tcl"]

namespace eval ::namdcph {
    namespace import ::cphTitrator::*

    variable TITLE "namdcph)"
    variable SystempH
    variable configFilenames [list]
    variable restartFilename ""
    variable restartFreq 0
    variable outFile ""
    variable numMinSteps 0
    variable excludeList [list]
    variable alchFrcCons 100.0
    variable useHeat 0
    variable MDBasename namdcph.md
    variable SWBasename namdcph.sw
    variable stateInfo [dict create] ;# cphSystem data and defaults
    variable moveInfo [dict create] ;# cphTitrator data and defaults

    namespace export *
}

# =============================================================================
# Main NAMD Routines
# =============================================================================
# ::namdcph::cphRun
#
# Run a set number of neMD/MC cycles.
#
proc ::namdcph::cphRun {numsteps {numcycles 1}} {
    variable ::namdcph::SystempH

    # Initialize NAMD and build a constant pH enabled PSF.
    initialize
    finalize
    if {$::namdcph::numMinSteps > 0} {
        set storedFirstTimeStep [firstTimeStep]
        minimize $::namdcph::numMinSteps
        if {[isset temperature]} {
          reinitvels [temperature]
        } else {
          reinitvels [$::thermostatTempCmd]
        }
        firstTimeStep $storedFirstTimeStep
    }
    set cphlog [openCpHLog]
    set firstCycle 1
    set lastCycle [expr {$firstCycle + $numcycles - 1}] 
    set runArgs [list $numsteps]
    for {set cycle $firstCycle} {$cycle <= $lastCycle} {incr cycle} { 
        # (1) Perform whatever equilibrium sampling is desired.
        #
        runMD {*}$runArgs
        # (2) Propose a move from the full move set.
        #
        lassign [cphTitrator propose $SystempH] accept swNumsteps\
                segresidnameList
        # (3) At this point we have either selected a switch or rejected a
        # whole bunch of moves. If the former, perform the switch.
        #
        if {!$accept} {
            cphPrint "Proposal rejected."
            set runArgs [list norepeat $numsteps]
        } else {
            cphPrint "Proposal accepted, attemping a switch."
            set accept [runSwitch $swNumsteps $segresidnameList]
            set runArgs [list $numsteps]
            # Only accumulate statistics for attempted switches.
            set moveLabel [join $segresidnameList "/"]
            cphTitrator accumulateStats $accept $moveLabel
        }
        cphSystem update $accept $segresidnameList
        # (4) Output cycle information and a checkpoint file if necessary.
        #
        puts $cphlog "[format "%6d" $cycle] [cphSystem get occupancy]"
        flush $cphlog
        writeRestart "[outputname].cphrst" $cycle
    }
    writeRestart force "[outputname].cphrst" $cycle
    # Cleanup temporary files
    file delete {*}[glob [getSWBasename].*]
    file delete {*}[glob [getMDBasename].*]
    close $cphlog
    printTitratorReport
    return
}

#
# FOR ADVANCED USE ONLY!!!
#
# ::namdcph::checkResidueDefinitions
#
# Check that appropriate residue definitions are defined for this residue.
#
# This essentially just scans all possible states for each titratable residue,
# of a given type in the system. THIS OFFERS NO VALIDATION AS TO WHETHER OR NOT
# THE DEFINITIONS MAKE SENSE, it simply checks that the exist.
#
proc ::namdcph::checkResidueDefinitions {resnames} {
    cphWarnExpt
    pH 7.0
    cphForceConstant 0.0
    # Initialize NAMD and build a constant pH enabled PSF.
    initialize
    finalize
    $::thermostatTempCmd 0.0
    outputEnergies 1
    alchLambdaFreq 0
    foreach segresidname [cphSystem get segresidnames] {
        lassign [split $segresidname ":"] segid resid resname
        if {$resname ni $resnames} {
            continue
        } 
        set states [cphSystem get stateList $segresidname]
        foreach state0 $states {
            foreach state1 $states {
                if {[string match $state0 $state1]} {
                    continue
                }
                # Force the initial state as state0
                output [getMDBasename]
                psfgenRead [getMDBasename]
                cphSystem initializeState $segresidname $state0
                alchemify $segresidname
                dealchemify $segresidname
                cphSystem update 1 $segresidname
                # Run an instantaneous switch
                runTestSwitch $segresidname $state1
            }
        }
    }
    return 
}

#
# FOR ADVANCED USE ONLY!!!
#
# ::namdcph::testResidueSwitch
#
# Run an instantenous switch for the given residue between the given states.
#
# This is meant as an interface function for checking energies against
# non-constant pH calculations, but cannot perform that comparison on its own.
#
proc ::namdcph::testResidueSwitch {segresidname state0 state1} {
    cphWarnExpt
    pH 7.0
    cphNumstepsPerSwitch 20
    cphForceConstant 0.0
    # Initialize NAMD and build a constant pH enabled PSF.
    cphSetResidueState $segresidname $state0
    initialize
    finalize
    $::thermostatTempCmd 0.0
    outputEnergies 1
    alchLambdaFreq 0
    return [runTestSwitch $segresidname $state1]
}

#
# FOR ADVANCED USE ONLY!!!
#
# ::namdcph::runTestSwitch
#
# This is just an internal convenience function for running instantaneous test
# switches.
#
# See checkResidueDefinitions and testResidueSwitch
#
proc ::namdcph::runTestSwitch {segresidname state1} {
    # You can't make lists of arrays or arrays of arrays, so the return type
    # has to be a dict of arrays (callback only returns arrays).
    set retVals [dict create]
    run 0
    storeEnergies
    dict set retVals MDEnergies0 [array get ::energyArray]
    cphSystem set trialState $segresidname $state1
    alchemify $segresidname
    alchLambda 0.0
    run 0
    storeEnergies
    dict set retVals SWEnergies0 [array get ::energyArray]
    alchLambda 1.0
    run 0
    storeEnergies
    dict set retVals SWEnergies1 [array get ::energyArray]
    dealchemify $segresidname
    cphSystem update 1 $segresidname
    run 0
    storeEnergies
    dict set retVals MDEnergies1 [array get ::energyArray]
    # Cleanup temporary files
    file delete {*}[glob [getSWBasename].*]
    file delete {*}[glob [getMDBasename].*]
    return $retVals
}

#
# FOR ADVANCED USE ONLY!!!
#
# ::namdcph::cphAnalyzeForce
#
# Test one or more residue definitions.
#
proc ::namdcph::cphAnalyzeForce {dcdfilename segresidname state0 state1} {
    cphWarnExpt
    pH 7.0

    cphSetResidueState $segresidname $state0
    initialize
    set initialPSF [structure]
    set initialPDB [coordinates]
    finalize

    set numsteps 500
    outputEnergies [expr {$numsteps == 0 ? 1 : int($numsteps)}]
    alchLambdaFreq [expr {$numsteps == 0 ? 0 : 1}]

    set nframes -1
    coorfile open dcd $dcdfilename
    while {![coorfile read]} {
        incr nframes

        # We have to do this so that inputs can be correctly loaded...
        set basename [format "%s.%d" [getMDBasename] $nframes]
        output $basename
        file copy -force $initialPSF "$basename.psf"
        file copy -force $initialPDB "$basename.pdb"
        reloadAndReinit $basename false
        # Assign the correct state and build protons or dummy atoms
        cphSystem initializeState $segresidname $state0
        alchemify $segresidname
        dealchemify $segresidname
        cphSystem update 1 $segresidname
        # Now build the alchemical atoms
        cphSystem set trialState $segresidname $state1
        alchemify $segresidname
        firsttimestep 0
        run $numsteps
        if {$numsteps} {
            storeEnergies
            set DeltaE [cphSystem compute switch $segresidname]
            set Work [expr {$::energyArray(CUMALCHWORK) + $DeltaE}]
            set ReducedWork [expr {$Work / [::kBT]}]
            set tmp [printProposalSummary $segresidname]
            cphPrint [format "%s WorkCorr % 10.4f CorrWork % 10.4f"\
                    [join $tmp "/"] $DeltaE $Work]
        }
        dealchemify $segresidname
        file delete {*}[glob $basename.*]
    }
    file delete {*}[glob [getMDBasename].*]
    file delete {*}[glob [getSWBasename].*]
    coorfile close
    return
}

# =============================================================================
# "Setter" Routines - used as new keywords in NAMD
#
#   All of these procedures take the desired value as an argument and return 
# that value. For default values see the decalaration in the namespace 
# header.
# =============================================================================
# -----------------
# Required Keywords
# -----------------
# ::namdcph::pH
#
# pH value for the simulation
#
proc ::namdcph::pH {pHValue} {
    checkIsNumeric "pH" $pHValue
    variable ::namdcph::SystempH $pHValue 
    return
}

# ::namdcph::cphConfigFile
#
# Configuration file for constant pH residues.
#
proc ::namdcph::cphConfigFile {filename} {
    variable ::namdcph::configFilenames
    lappend configFilenames [string trim $filename]
    return
}

# ::namdcph::cphNumstepsPerSwitch
#
# For odd arguments, the first argument is assumed to be a default switch time.
# All remaining arguments are presumed to be label/numsteps pairs for specific
# moves.
#
proc ::namdcph::cphNumstepsPerSwitch {args} {
    variable ::namdcph::moveInfo
    if {[expr {[llength $args] % 2}]} {
        set numsteps [lindex $args 0]
        checkIsNotNegative numSwitchSteps $numsteps
        dict set moveInfo default numsteps [expr {int($numsteps)}]
        set args [lrange $args 1 end]
    }
    checkArglistIsMultiple $args 2
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        lassign [lrange $args $i [expr {$i+1}]] moveLabel numsteps
        checkIsNotNegative numSwitchSteps $numsteps
        dict set moveInfo $moveLabel numsteps [expr {int($numsteps)}]
    }
    return
}

# ---------------------
# Commonly Used Options
# ---------------------
# ::namdcph::cphSetResidueState
#
# Set the state of one or more residues using psfgen syntax.
#
proc ::namdcph::cphSetResidueState {args} {
    checkArglistIsMultiple $args 2
    variable ::namdcph::stateInfo
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        lassign [lrange $args $i [expr {$i+1}]] segresidname state
        dict set stateInfo $segresidname state $state
    }
    return
}

# ::namdcph::cphSetResiduepKai
#
# Set the inherent pKa of one or more residues using psfgen syntax.
#
proc ::namdcph::cphSetResiduepKai {args} {
    checkArglistIsMultiple $args 2
    variable ::namdcph::stateInfo
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        #NB pKai may be a list of values - check individually.
        lassign [lrange $args $i [expr {$i+1}]] segresidname pKai
        foreach pKa $pKai {
            checkIsNumeric pKai $pKa
        }
        dict set stateInfo $segresidname pKai $pKai
    }
    return
}

# ::namdcph::cphRestartFile 
#
# Restart a constant pH run from a restart file.
#
proc ::namdcph::cphRestartFile {filename} {
    variable ::namdcph::restartFilename $filename
    return
}

# ::namdcph::cphRestartFreq
#
# Frequency (in neMD/MC cycles) at which to save constant pH restart files
#
proc ::namdcph::cphRestartFreq {frequency} {
    checkIsNotNegative cphRestartFreq $frequency
    variable ::namdcph::restartFreq $frequency
    return
}

# ::namdcph::cphOutFile
#
# Name for constant pH output file - default is [outputname].cphlog
#
proc ::namdcph::cphOutFile {filename} {
    variable ::namdcph::outFile $filename
    return
}

# ::namdcph::cphProposalWeight
#
# The (unnormalized) proposal weight assigned to each move. The default is for
# all such weights to be equal (uniformally distributed moves).
#
proc ::namdcph::cphProposalWeight {args} {
    variable ::namdcph::moveInfo
    checkArglistIsMultiple $args 2
    for {set i 0} {$i < [llength $args]} {incr i 2} {
        lassign [lrange $args $i [expr {$i+1}]] moveLabel weight
        checkIsNotNegative proposalWeight $weight
        dict set moveInfo $moveLabel weight [expr {1.*$weight}]
    }
    return
}

# ::namdcph::cphProposeProtonTransfer
#
# Add additional "proton transfer" moves to the move set. These new moves are
# only attempted if a proton can be swapped between the given residues. The
# two residues are given by their normal <segid>:<resid>:<resname> designation
# except that they are further separated by a "/".
#
proc ::namdcph::cphProposeProtonTransfer {args} {
    variable ::namdcph::moveInfo
    foreach moveLabel $args {
        dict set moveInfo $moveLabel protonTransfer 1
    }
    return
}

# ::namdcph::cphProposeCotitration
#
# Add additional "co-titration" moves to the move set. These new moves are only
# attempted if two residues are both protonated or deprotonated. They will then
# attempt to change states concurrently. The two residues are given by their
# normal <segid>:<resid>:<resname> designation except that they are further
# separated by a "/".
#
proc ::namdcph::cphProposeCotitration {args} {
    variable ::namdcph::moveInfo
    foreach moveLabel $args {
        dict set moveInfo $moveLabel cotitration 1
    }
    return
}

# ::namdcph::cphMaxProposalAttempts
#
# Number of attempted MC proposals from each move set before giving up.
#
# Values less than 1 default to the number of residues in the system.
# NB: This is not the necessarily the same thing as attempting a move once for
#     each residue.
#
proc ::namdcph::cphMaxProposalAttempts {maxAttempts} {
    checkIsNumeric cphMaxProposalAttempts $maxAttempts
    variable ::namdcph::moveInfo
    dict set moveInfo maxProposalAttempts [expr {int($maxAttempts)}]
    cphWarn "cphMaxProposalAttempts is now deprecated and set to 1!"
    return
}

# ::namdcph::cphNumMinSteps
#
# Number of minimization steps to perform before dynamics.
#
# This is especially useful when the initial states are randomized according to
# the pH.
#
proc ::namdcph::cphNumMinSteps {numsteps} {
    checkIsNotNegative cphNumMinSteps $numsteps
    variable ::namdcph::numMinSteps $numsteps
    return
}

proc ::namdcph::cphExcludeResidue {args} {
    variable ::namdcph::excludeList
    foreach segresidname $args {
         lappend excludeList $segresidname
    }
    return
}

# -------------------
# Specialized Options
# -------------------
# ::namdcph::cphForceConstant
#
# Force constant for zero-length bonds between alchemical atoms in kcal/mol-A^2
#
proc ::namdcph::cphForceConstant {forceConstant} {
    checkIsNotNegative cphForceConstant $forceConstant
    variable ::namdcph::alchFrcCons $forceConstant
    return
}

# ::namdcph::cphMDBasename
#
# Basename for (temporary) constant pH MD input files - the only output of
# interest to the user should be in the usual places
#
proc ::namdcph::cphMDBasename {basename} {
    variable ::namdcph::MDBasename $basename
    return
}

# ::namdcph::cphSwitchBasename
#
# Basename for (temporary) constant pH switch trajectory output files
#
proc ::namdcph::cphSwitchBasename {basename} {
    variable ::namdcph::SWBasename $basename
    return
}

# =============================================================================
# Nonequilibrium Switching Routines
# =============================================================================
# ::namdcph::runMD
proc ::namdcph::runMD {args} {
   run {*}$args 
}

# ::namdcph::runSwitch
#
# Run a neMD/MC switch on the specified residue(s).
#
# This includes multiple steps:
# (1) Modify I/O settings to differentiate from regular MD
# (2) Build alchemically enabled side chains for the desired residues
# (3) Perform a switching trajectory with appropriate momentum reversals
# (4) Perform Metropolis MC using an appropriate work quantity
# (5) Reinitialize MD based on the MC result
#
# Arguments:
# ----------
# numsteps : int
#   Number of steps in the switch
# segresidnameList : list 
#   One or more "<segid>:<resid>" specifications - this is the same syntax as 
#   for the regular psfgen patch command.
#
# Returns:
# --------
# accept : boolean
#   Result of MC accept/reject test 
#
proc ::namdcph::runSwitch {numsteps segresidnameList} {
    variable ::namdcph::useHeat

    # (1) Checkpoint and modify output parameters. 
    #
    storeEnergies
    set storedOutputEnergies [outputEnergies]
    set storedDCDFreq [dcdFreq]
    set storedTimestep $::energyArray(TS)
    outputEnergies $numsteps
    dcdFreq 0
    firstTimestep 0
    if {$::barostatIsSet && [$::barostatCmd]} {
        $::barostatCmd off
    }
    # (2) Build the alchemical switch inputs.
    #
    alchemify $segresidnameList
    # (3) Run the switch trajectory with momentum reversal.
    #
    runprpswitch $numsteps
    outputEnergies $storedOutputEnergies
    dcdFreq $storedDCDFreq
    # (4) Compute the work with state dependent energy shifts.
    #
    storeEnergies
    set DeltaE [cphSystem compute switch $segresidnameList]
    if {$useHeat} {
        # Computed as totalEnergy - totalEnergy0 - heat.
        set Work $::energyArray(WORK)
    } else {
        # Computed as alchemical protocol work.
        set Work $::energyArray(CUMALCHWORK)
    }
    set ReducedWork [expr {($Work + $DeltaE) / [::kBT]}]
    set accept [metropolisAcceptance $ReducedWork]
    set tmp [printProposalSummary $segresidnameList]
    cphPrint [format "%s WorkCorr % 10.4f CorrWork % 10.4f"\
            [join $tmp "/"] $DeltaE $Work]
    # (5) Reinitialize for the next MD step based on accept/reject.
    #
    outputEnergies $storedOutputEnergies
    dcdFreq $storedDCDFreq
    firstTimestep $storedTimestep
    if {$::barostatIsSet && ![$::barostatCmd]} {
        $::barostatCmd on
    }

    if {$accept} {
        cphPrint "Switch accepted!"
        dealchemify $segresidnameList
    } else {
        cphPrint "Switch rejected!"
        alch off
        reloadAndReinit [getMDBasename] false
    }
    return $accept
}

# ::namdcph::alchemify
#
# Reinitialize the system with alchemical sidechains for the given residues.
#
# This includes multiple steps:
# (1) The NAMD state is written to disk and re-read by PSFGEN.
# (2) Appropriate alchemical patches are applied and alchemical atom 
#     coordinates and velocities are sampled.
# (3) New NAMD inputs are written and then re-read. This includes new 
#     extraBonds for the alchemical sidechain.
#
# Arguments:
# ----------
# segresidnameList : list of strings
#   One or more "<segid>:<resid>" specifications - this is the same syntax as 
#   for the regular psfgen patch command.
#
# Returns:
# --------
# None
#
proc ::namdcph::alchemify {segresidnameList} {  
    variable ::namdcph::alchFrcCons
    set T [$::thermostatTempCmd]

    alch on
    # (1) Read in the nonalchemical PSF and apply patches.
    #
    output [getMDBasename]
    psfgenRead [getMDBasename]
    # (2) Apply patches and build coordinates and velocities.
    #
    foreach segresidname $segresidnameList {
        cphSystem alchemifypsf $segresidname $alchFrcCons $T
    }
    regenerate angles dihedrals
    # (3) Write a new set of inputs and reinitialize.
    #
    psfgenWrite [getSWBasename] [mdFilename xsc]
    # Dummy atoms have been built and a PDB has been written. We can now query
    # atom indices and build extraBonds. If psfgenWrite were _not_ called, then
    # the atomid queries would all return zero (that would be bad).
    #
    set ExtraBondsFile [open [swFilename extrabonds] "w"]
    foreach segresidname $segresidnameList {
        puts $ExtraBondsFile\
                [cphSystem get alchBonds $segresidname $alchFrcCons]
    }
    close $ExtraBondsFile
    reloadAndReinit [getSWBasename] true
    return
}

# ::namdcph::dealchemify
#
# Remove alchemical side
#
# Arguments:
# ----------
# segresidnameList : list of strings
#   One or more "<segid>:<resid>" specifications - this is the same syntax as 
#   for the regular psfgen patch command.
#
# Returns:
# --------
# None
#
proc ::namdcph::dealchemify {segresidnameList} {
    output [getSWBasename]
    psfgenRead [getSWBasename]
    foreach segresidname $segresidnameList {
        cphSystem dealchemifypsf $segresidname
    }
    psfgenWrite [getMDBasename] [swFilename xsc]
    alch off
    reloadAndReinit [getMDBasename] false
    return
}

# =============================================================================
# Constant pH specific I/O
# =============================================================================
# ::namdcph::readConfig
#
# Read one or more constant pH config files and return template definitions.
#
# Arguments:
# ----------
# None
#
# Returns:
# --------
# templateDefs : dict
#   A dict conversion of the config file contents from JSON
#
proc ::namdcph::readConfig {} {
    variable ::namdcph::configFilenames

    set templateDefs [dict create]
    foreach configFilename $configFilenames {
        set ConfigFile [open $configFilename "r"]
        set tmpDict [json::json2dict [read $ConfigFile]]
        # Warn the user about re-definitions, but proceed anyway. Note that
        # dict merge overwrites values right to left.
        set oldKeys [dict keys $templateDefs]
        foreach newKey [dict keys $tmpDict] {
            if {$newKey in $oldKeys} {
                cphWarn "Reading and using new definition for residue $newKey"
            }
        } 
        set templateDefs [dict merge $templateDefs $tmpDict] 
        close $ConfigFile
    }
    return $templateDefs
}

# ::namdcph::readRestart
#
# Read in a constant pH state from file and _merge_ the data into the namespace
# variables. That is, give precedence to keyword specifications.
#
# NB! This also requires that cphSystem has been built so that the system and
#     state info can be compared.
#
proc ::namdcph::readRestart {} {
    variable ::namdcph::restartFilename
    variable ::namdcph::SystempH
    variable ::namdcph::moveInfo

    set RestartFile [open $restartFilename "r"]
    set Restart [json::json2dict [read $RestartFile]]
    close $RestartFile

    if {![info exists SystempH] && [dict exists $Restart pH]} {
        pH [dict get $Restart pH]
    }

    if {[dict exists $Restart exclude]} {
        cphExcludeResidue {*}[dict get $Restart exclude]
    }

    # Read state parameters, if present. 
    #
    set stateList {}
    if {[dict exists $Restart states]} {
        set stateList [dict get $Restart states]
    }
    set pKaiList {}
    if {[dict exists $Restart pKais]} {
        set pKaiList [dict get $Restart pKais]
    }

    # Read MC move parameters, if present. Only reset those values which have
    # not been set in via keywords.
    #
    if {[dict exists $Restart MCmoves]} {
        set rstrtMoveInfo [dict get $Restart MCmoves]    
    }
    if {[info exists rstrtMoveInfo]} {
        # This is a bit kludgy - the global data is not a nested dict, so it
        # will crash in the dict for loop. Treat these specially and then unset
        # them...
        if {[dict exists $rstrtMoveInfo maxProposalAttempts]} {
            if {![dict exists $moveInfo maxProposalAttempts]} {
                dict set moveInfo maxProposalAttempts\
                    [dict get $rstrtMoveInfo maxProposalAttempts]
            }
            dict unset rstrtMoveInfo maxProposalAttempts
        }
        dict for {moveLabel data} $rstrtMoveInfo {
            dict for {key value} $data {
                 if {![dict exists $moveInfo $moveLabel $key]} {
                     dict set moveInfo $moveLabel $key $value
                 }
            }
        }
    }

    return [list $stateList $pKaiList]
}

# ::namdcph::writeRestart
#
# Write the current constant pH state information to restartFilename based on 
# the restart frequency. If the "force" keyword precedes the arguments, ignore 
# the restart frequency and write no matter what.
#
# Arguments:
# ----------
# restartFilename : string
#   Name of (JSON format) restart file to write
# cycle : int
#   Last neMD/MC cycle index
#
# Returns:
# --------
# None
#
proc ::namdcph::writeRestart {args} { 
    variable ::namdcph::restartFreq
    variable ::namdcph::SystempH
    variable ::namdcph::excludeList

    if {[string match [lindex $args 0] force]} {
        set restartFilename [lindex $args 1]
        set cycle [expr {int([lindex $args 2])}]
    } else {
        set restartFilename [lindex $args 0]
        set cycle [expr {int([lindex $args 1])}]
        if {!($restartFreq) || !($cycle) || [expr {$cycle % $restartFreq}]} {
            return
        }
    }

    # Here we write the restart as a dict in JSON format, however, this is NOT
    # the same as converting a Tcl dict to JSON, because knowledge of types is
    # lost in the conversion. The easiest (and maybe faster?) way is to build a
    # list that contains all of the dict key/value pairs and then emulate the
    # JSON format with commas and curly braces.
    #
    set rstrtList [list]
    lappend rstrtList "\"cycle\":$cycle"
    lappend rstrtList "\"pH\":$SystempH"
    if {[llength $excludeList]} {
        lappend rstrtList "\"exclude\":[json::list2json $excludeList 1]"
    }
    lappend rstrtList "\"states\":[json::list2json [cphSystem get state] 1]"
    lappend rstrtList "\"pKais\":[json::list2json [cphSystem get pKai]]"
    lappend rstrtList "\"MCmoves\":[json::dict2json [cphTitrator get archive]]"

    # Write everything to file.
    namdFileBackup $restartFilename
    set RestartFile [open $restartFilename "w"]
    puts $RestartFile "\{[join $rstrtList ,]\}"
    close $RestartFile

    # Always checkpoint the PSF and PDB when a restart is written.
    file copy -force [mdFilename psf] "[outputName].psf"
    file copy -force [mdFilename pdb] "[outputName].pdb"
    return
}

# ::namdcph::openCpHLog
#
# Open a new constant pH log for proton occupancies and return the file object.
#
proc ::namdcph::openCpHLog {} {
    variable ::namdcph::outFile
    variable ::namdcph::SystempH

    if {[string length $outFile] > 0} {
        set logFilename $outFile
    } else {
        set logFilename "[outputname].cphlog"
    }
    namdFileBackup $logFilename
    set cphlog [open $logFilename "w"]
    puts $cphlog "#pH $SystempH"
    puts $cphlog "#[join [cphSystem get segresidnames] " "]"
    return $cphlog 
}

# =============================================================================
# Convenience Routines 
# =============================================================================
# ::namdcph::swFilename
#
# Get an appropriate filename for a temporary file used during switches
#
proc ::namdcph::swFilename {ext} {
    return "[getSWBasename].$ext"
}

# ::namdcph::mdFilename
#
# Get an appropriate filename for a temporary file used to launch regular MD
#
proc ::namdcph::mdFilename {ext} {
    return "[getMDBasename].$ext"
}

# ::namdcph::clearExtraBonds
#  
# Clear all bonds in the extraBondsFile.
#
proc ::namdcph::clearExtraBonds {} {
    set ExtraBondsFilename [swFilename extrabonds]
    cphPrint "clearing extraBondsFile $ExtraBondsFilename"
    set ExtraBondsFile [open $ExtraBondsFilename "w"]
    puts $ExtraBondsFile ""
    close $ExtraBondsFile
    return
}

# ::namdcph::reloadAndReinit
#
# Read in a new PSF/PDB pair and then reinitialize atoms using a basename read.
#
proc ::namdcph::reloadAndReinit {basename keepExtraBonds} {
    if {!$keepExtraBonds} {
        clearExtraBonds
    }
    structure "$basename.psf" pdb "$basename.pdb"
    reinitatoms $basename
    return
}

proc ::namdcph::cphPrint {msg} {
    print "$::namdcph::TITLE $msg"
}

proc ::namdcph::cphWarn {msg} {
    print "$::namdcph::TITLE WARNING! $msg"
}

proc ::namdcph::cphAbort {{msg ""}} {
    abort "$::namdcph::TITLE $msg"
}

proc ::namdcph::cphWarnExpt {} {
    cphWarn "THIS FEATURE IS EXPERIMENTAL!"
    cphPrint "RESULTS ARE NOT GUARANTEEED - USE AT YOUR OWN RISK"
}

# =============================================================================
# Setup Routines
# =============================================================================
# ::namdcph::initialize
#
# Initialize the system for constant pH. This requires two main things to
# happen:
#
# 1) nonequilibrium alchemical transformations must be enabled
# 2) the PSF/PDB must be rebuilt to include dummy atoms and possibly reassigned
#    protonation states
#
proc ::namdcph::initialize {} {
    variable ::namdcph::restartFilename
    variable ::namdcph::configFilenames
    variable ::namdcph::stateInfo
    variable ::namdcph::moveInfo
    variable ::namdcph::excludeList
    variable ::namdcph::SystempH
    variable ::namdcph::useHeat

    callback energyCallback
    # 1) Set up alchemical keywords and run startup. 
    #
    initializeAlch
    getThermostat
    getBarostat
    # Note that checking the stored name avoids the need to catch in case we
    # have an older version of NAMD where stochRescale is unavailable.
    if {[string match $::thermostatName "stochastic-rescaling"]} {
        if {[isset stochRescaleHeat] && [stochRescaleHeat]} {
            set useHeat 1
        }
    }

    # 2) Rebuild the PSF with dummy protons and modify protonation states as 
    # needed. Build the residue definitions, assign states to each residue, and 
    # rebuild the topology to reflect those states.
    #
    cphPrint "initializing constant pH PSF..."
    if {![llength configFilenames]} {
        cphAbort "At least one constant pH configuration file is required."
    }

    if {[string length $restartFilename]} {
        lassign [readRestart] stateList pKaiList
        lassign [list 0.0 false] temp buildH
    } else {
        lassign [list [$::thermostatTempCmd] false] temp buildH
    }
    cphSystem build [readConfig] $excludeList

    # Now that we've built the system, we can allocate the state and pKa
    # information from the restart (if present).
    if {[info exists stateList] && [llength $stateList] &&
        [info exists pKaiList] && [llength $pKaiList]} {
        set segresidnameList [cphSystem get segresidnames]

        if {[llength $stateList] != [llength $pKaiList]} {
            cphAbort "mismatch in states/pKais in $restartFilename"
        }
        if {[llength $stateList] != [llength $segresidnameList]} {
            cphAbort "Too few/many state/pKai definitions in $restartFilename"
        }

        foreach segresidname $segresidnameList\
                state $stateList pKai $pKaiList {
            if {![dict exists $stateInfo $segresidname]} {
                dict set stateInfo $segresidname state $state
                dict set stateInfo $segresidname pKai $pKai
            } else {
                if {![dict exists $stateInfo $segresidname state]} {
                    dict set stateInfo $segresidname state $state
                }
                if {![dict exists $stateInfo $segresidname pKai]} {
                    dict set stateInfo $segresidname pKai $pKai
                }
            }
        }
    }

    if {![info exists SystempH]} {
        cphAbort "A pH value is required."
    }

    # All residue state info from all sources is now in stateInfo.
    # Check that all of the keys are valid.
    foreach segresidname [dict keys $stateInfo] {
        if {[cphSystem validate $segresidname]} {
            cphAbort
        }
    }
    cphSystem initialize $SystempH $temp $buildH $stateInfo

    # 3) Build the MC move set (the "titrator").
    #
    # All neMD/MC move info from all sources is now in moveInfo.
    cphTitrator build $moveInfo

    # 4) Report to stdout.
    printSettingsSummary
    printSystemSummary
    printTitratorSummary
    return
}

proc ::namdcph::finalize {} {
    # 5) Write to disk and prepare for MD.
    if {[isset extendedSystem]} {
        set inputXSCName [extendedSystem]
    } else {
        set inputXSCName [mdFilename xsc]
        set inputXSC [open $inputXSCName "w"]
        puts $inputXSC $::dummyXSC 
        close $inputXSC
    }
    psfgenWrite [getMDBasename] $inputXSCName 
    reloadAndReinit [getMDBasename] false
    if {![isset binVelocities] && ![isset velocities]} {
        reinitvels [temperature]
    }
    return
}

# ::namdcph::initializeAlch
#
#   Initialize the alchemical settings needed for nonequilibrium alchemical
# switching of protonation states. This activates the appropriate nonbonded
# kernels and stores SimParameters settings that will be needed during 
# switches.
#
proc ::namdcph::initializeAlch {} {
    if {[isset alch]} {
        cphAbort "Constant pH is currently incompatible with alchemical"\
                 "transformations. Remove all alch settings."
    }
    cphPrint "Setting up nonequilibrium alchemical switching."
    # Here we require alchemical force computations and complete decoupling of
    # _all_ interactions. Unstaged linear coupling without softcore shifting is
    # currently assumed/enforced, but this might be ok to change.
    #
    alch on
    alchType TI
    alchDecouple off
    alchElecLambdaStart 0.0
    alchVdwLambdaEnd 1.0
    alchVdwShiftCoeff 0.0
    alchBondLambdaEnd 1.0
    alchBondDecouple on
    alchLambda 0.0
    alchLambda2 1.0
    alchOutFreq 0 ;# Suppress output - this would just be a nightmare.
    # Bonds between corresponding alchemical atoms are built dynamically when
    # a switch is started. In principle, this does not conflict with any other
    # extraBonds settings. However, indices in the extraBondsFile are likely
    # corrupted and would lead to wildly unexpected behavior - disable this for
    # now.
    #
    # TODO: Read and modify these files on they fly?
    if {[isset extraBonds]} {
        cphAbort "Constant pH is currently incompatible with extraBonds."
    }
    extraBonds on
    extraBondsFile [swFilename extrabonds]
    clearExtraBonds
    startup
    alchLambdaFreq 1
    alch off
    return
}

# ::namdcph::printSettingsSummary
#
# Print a summary of the constant pH settings.
#
proc ::namdcph::printSettingsSummary {} {
    variable ::namdcph::configFilenames
    variable ::namdcph::restartFilename
    variable ::namdcph::SystempH
    variable ::namdcph::alchFrcCons
    variable ::namdcph::useHeat

    set StarBar "***************************************"
    cphPrint $StarBar
    cphPrint "CONSTANT pH MD ACTIVE"
    if {[string length $restartFilename] > 0} {
        cphPrint "RESTART FILENAME $restartFilename"
    }
    cphPrint "SYSTEM pH $SystempH"
    cphPrint "CONSTANT pH CONFIGURATION FILE(S)"
    foreach configFilename $configFilenames {
        cphPrint "$configFilename"
    }
    cphPrint "NONEQUILIBRIUM SWITCH PARAMETERS:"
    cphPrint "ALCHEMICAL FORCE CONSTANT $alchFrcCons kcal/mol-A^2"
    cphPrint "neMD/MC CRITERION TEMPERATURE [$::thermostatTempCmd]"
    if {$useHeat} {
        cphPrint "WORK COMPUTATION: HEAT EXCHANGE W/THERMOSTAT"
        cphWarnExpt
    } else {
        cphPrint "WORK COMPUTATION: ALCHEMICAL PROTOCOL WORK"
    }
    cphPrint "TEMPORARY FILE INFORMATION:"
    cphPrint "cpH TOPOLOGY FILE BASENAME [getMDBasename]"
    cphPrint "neMD/MC TRAJECTORY BASENAME [getSWBasename]"
    cphPrint $StarBar
    return
}

# ::namdcph::printSystemSummary
#
# Print a summary of the titratable system. 
#
proc ::namdcph::printSystemSummary {} {
    set StarBar "***************************************"
    cphPrint $StarBar
    cphPrint "TITRATABLE RESIDUE DEFINITIONS:"
    cphPrint "[join [cphSystem get resdefs] " "]"
    cphPrint "TITRATABLE SYSTEM SUMMARY:"
    cphPrint "[cphSystem get numresidues] RESIDUE(S)"
    cphPrint "[llength [cphSystem get occupancy]] PROTONATION SITE(S)"
    cphPrint $StarBar
    cphPrint [format "%-19s : %5s : %s" segid:resid:resname state pKai]
    foreach segresidname [cphSystem get segresidnames]\
            state [cphSystem get state]\
            pKaiList [cphSystem get pKai] {
        cphPrint [format "%-19s : % 5s : %-s" $segresidname $state $pKaiList]
    }
    cphPrint $StarBar
    return
}

# ::namdcph::printTitratorSummary
#
# Print a summary of the MC moves.
#
proc ::namdcph::printTitratorSummary {} {
    set StarBar "*************************************************"
    set moveLabels [cphTitrator get moveLabels]
    set numMoves [llength $moveLabels]

    cphPrint $StarBar
    cphPrint "CONSTANT pH neMD/MC MOVES:"
    cphPrint "MAX. ATTEMPTS PER CYCLE: [cphTitrator get maxAttempts]"

    cphPrint "ONE RESIDUE MOVES"
    cphPrint [format "%-14s : %8s %8s" "move label" numsteps weight]
    foreach moveLabel $moveLabels {
        if {[llength [split $moveLabel "/"]] != 1} continue
        incr numEntries
        set numsteps [cphTitrator get numsteps $moveLabel]
        set weight [cphTitrator get weight $moveLabel]
        cphPrint [format "%-14s : % 8d % 8.2f" $moveLabel $numsteps $weight]
    }
    if {$numEntries == $numMoves} {
        cphPrint $StarBar
        return
    }

    cphPrint "TWO RESIDUE MOVES : proton transfer (PT) : co-titration (CT)"
    cphPrint [format "%-29s : %8s %8s %2s %2s"\
            "move label" numsteps weight PT CT]
    foreach moveLabel $moveLabels {
        if {[llength [split $moveLabel "/"]] != 2} continue
        set numsteps [cphTitrator get numsteps $moveLabel]
        set weight [cphTitrator get weight $moveLabel]
        set PT [cphTitrator get protonTransfer $moveLabel]
        set CT [cphTitrator get cotitration $moveLabel]
        cphPrint [format "%-29s : % 8d % 8.2f %2d %2d"\
                $moveLabel $numsteps $weight $PT $CT]
    }
    cphPrint $StarBar
    return
}

# ::namdcph::printProposalSummary
proc ::namdcph::printProposalSummary {segresidnameList} {
    set retList [list]
    foreach segresidname $segresidnameList {
        set state [cphSystem get state $segresidname]
        set trialState [cphSystem get trialState $segresidname]
        lappend retList "$segresidname:$state:$trialState"
    }
    return $retList
}

# ::namdcph::printTitratorReport
#
# Print a report of the titration MC statistics.
#
proc ::namdcph::printTitratorReport {} {
    set StarBar "*************************************************"
    cphPrint $StarBar
    cphPrint "CONSTANT pH MD STATISTICS:"
    cphPrint [format "%-25s : %-8s %-12s" "move label" attempts "accept. rate"]
    foreach moveLabel [cphTitrator get moveLabels] {
        set attempts [cphTitrator get attempts $moveLabel]
        set successes [cphTitrator get successes $moveLabel]
        if {$successes && $attempts} {
            set rate [expr {100.*$successes/$attempts}]
        } else {
            set rate 0.0
        }
        cphPrint [format "%-25s : %8d %12.2f" $moveLabel $attempts $rate] 
    }
    cphPrint $StarBar
    return
}

# =============================================================================
# Getter Routines
#
# These are largely unnecessary, but cut down on "variable" declarations.
# =============================================================================
proc ::namdcph::getMDBasename {} {
    return $::namdcph::MDBasename
}

proc ::namdcph::getSWBasename {} {
    return $::namdcph::SWBasename
}
