# cphtoppar.tcl
#
#   This file provides the ::cphSystem namespace, which effectively emulates an
# object containing all constant-pH specific topology and parameter (toppar)
# information.
#
package require Tcl 8.5

source [file join [file dirname [info script]] "cphpsfgen.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::cphSystem {
    #   The core information of the cphSystem is the resDict, which stores
    # the unique and mutable information for each instance of a titratable
    # residue. The keys are "segresidnames" of the form:
    #
    # <segid>:<resid>:<resname>.
    #
    # Note the <resname> is a _titratable residue_ name and may or may not be
    # the same as a recognizable residue from the force field. For example,
    # amino acid termini are titratable residues. The resDefDict stores
    # non-unique and immutable information for each residue definition. The
    # keys are residue names (e.g. ASP, HIS, etc.). Generally, when looking up
    # information for a specific residue, one simply needs look up the name in
    # resDict and then look up the info in resDefDict.
    #
    # Note that a <segid>:<resid> combination can correspond to _multiple_
    # residue names. This permits non-overlapping "sub-residues," such as amino
    # acid termini, which titrate independently of the sidechain.
    #
    variable resDict [dict create]
    variable resDefDict [dict create]

    namespace export cphSystem
}

# =============================================================================
# cphSystem Interface function
# =============================================================================
# ::cphSystem::cphSystem
#
# This is the only exported function from the cphSystem namespace and provides
# a complete interface to the cphSystem namespace. The first argument is always
# an action (e.g. "get" or "set") and the rest are variable arguments.
#
# action           description
# ---------------- -----------
# get              see cphSystemGet 
# set              see cphSystemSet
# build            see buildSystem
# initialize       see initializeSystem
# propose          propose a move of the given type
# compute          compute various quantities needed for MC
# update           accept/reject the proposed changes
# alchemifypsf     apply patches that make the system alchemical
# dealchemifypsf   remove patches so the system is no longer alchemical
# initializeState  set the initial state of a residue
#
# Some notes about Tcl switch statements for non-guru types like myself:
#
#   These can appear a bit delicate to those who are not familiar with their
# quirks (compared to C). For one, badly placed comments can break them. Also
# note that a switch returns the expression it first matches and does NOT use
# break statements.
#
proc ::cphSystem::cphSystem {action args} {
    return [switch -nocase -- $action {
        get {
            cphSystemGet {*}$args
        }
        set {
            cphSystemSet {*}$args
        }
        build {
            buildSystem {*}$args
        }
        initialize {
            initializeSystem {*}$args
        }
        propose { ;# return 0 if no valid proposal is found
            set type [lindex $args 0]
            set newArgs [lrange $args 1 end]
            switch -nocase -- $type {
                titration {
                    proposeResidueTitration {*}$newArgs
                }
                tautomerization {
                    proposeResidueTautomerization {*}$newArgs
                }
                protonTransfer {
                    proposeProtonTransfer {*}$newArgs
                }
                cotitration {
                    proposeCotitration {*}$newArgs
                }
                default {
                    abort "Invalid proposal type $type"
                }
            }
        }
        compute { ;# acceptance energy (inherent or switch correction)
            set type [lindex $args 0]
            set newArgs [lrange $args 1 end]
            switch -nocase -- $type {
                inherent {
                    computeInherentAcceptance {*}$newArgs
                }
                inherentWeights {
                    computeInherentNormedWeights {*}$newArgs
                }
                switch {
                    computeSwitchAcceptance {*}$newArgs
                }
                default {
                    abort "Invalid energy type $type"
                }
            }
        }
        update { ;# Update states based on an MC result
            updateStates {*}$args
        }
        alchemifypsf { ;#Alchemify the PSF (in memory) based on trial states
            alchemifyResidue {*}$args
        }
        dealchemifypsf { ;# Dealchemify the PSF (in memory)
            dealchemifyResidue {*}$args
        }
        initializeState { ;# Assign a state
            set method [lindex $args 0]
            switch -nocase -- $method {
                random {
                    randomizeState {*}[lrange $args 1 end]
                }
                default {
                    assignState {*}$args
                }
            }
       }
       validate { ;# Validate a segresidname
           validateSegresidname {*}$args
       }
       default {
           abort "Invalid cphSystem action $action."
       }
    }]
}

# ::cphSystem::validateSegresidname
#
# Check that a given segresidname is valid. That is:
#
# 1) is it in the correct <segid>:<resid>:<resname> format?
# 2) does it correspond to a known titratable residue?
#
# A non-zero error code is returned for each condition.
#
# NB! This should only be called _after_ buildSystem.
#
proc ::cphSystem::validateSegresidname {segresidname} {
    lassign [split $segresidname ":"] segid resid resname
    if {![info exists segid] || ![info exists resid]
        || ![info exists resname]} {
        print "segresidname selections must be of form\
                 <segid>:<resid>:<resname>!"
        return -1
    }
    if {$segresidname ni [cphSystem get segresidnames]} {
        print "Invalid segresidname $segresidname!"
        return -2
    }
    return 0
}

# =============================================================================
# Proposal Routines
# =============================================================================
# proc ::cphSystem::updateStates
#
# Update the state of one or more residues with the given
# <segid>:<resid>:<resname> specifications based on acceptance/rejection of the
# trial states. Reset the trial states to a null value.
#
proc ::cphSystem::updateStates {accept segresidnameList} {
    if {$accept} {
        foreach segresidname $segresidnameList {
            cphSystem set state $segresidname\
                    [cphSystem get trialState $segresidname]
            cphSystem set trialState $segresidname {}
        }
    } else {
        foreach segresidname $segresidnameList {
            cphSystem set trialState $segresidname {}
        }
    }
    return 0
}

# ::cphSystem::proposeResidueTitration
#
# Propose a new trial state requiring a titration - i.e. a net change in the 
# number of protons.
#
# For consistency with tautomers (see below), this always returns true in order
# to confirm that a titration was in fact found.
#
proc ::cphSystem::proposeResidueTitration {segresidname} {
    set possibleStates [cphSystem get trialStateList $segresidname]
    set numProtons [expr {lsum([cphSystem get occupancy $segresidname])}]
    while {true} {
        lassign [choice $possibleStates] state
        set occ [cphSystem get occupancy $segresidname $state]
        set numProtonsTrial [expr {lsum($occ)}]
        if {$numProtons != $numProtonsTrial} {
            cphSystem set trialState $segresidname $state
            return 1
        }
    }
    # This is an error and should never happen.
    return -1
}

# ::cphSystem::proposeProtonTransfer
#
# Propose a new trial state requiring a proton transfer - i.e. a movement of a
# proton from residue to another.
#
# This returns true if no proton transfer or else a negative error code.
#
proc ::cphSystem::proposeProtonTransfer {segresidname1 segresidname2} {
    set numProtons1 [expr {lsum([cphSystem get occupancy $segresidname1])}]
    set numProtons2 [expr {lsum([cphSystem get occupancy $segresidname2])}]
    set possibleStates1 [cphSystem get trialStateList $segresidname1]
    set possibleStates2 [cphSystem get trialStateList $segresidname2]
    set dn [expr {$numProtons1 - $numProtons2}]
    if {$dn == 0.0} {     
        # No proton transfer is possible.
        return 1
    } 
    # transfer from 1 to 2
    while {true} {
        lassign [choice $possibleStates1] state1
        set occ1 [cphSystem get occupancy $segresidname1 $state1]
        set numProtonsTrial1 [expr {lsum($occ1)}]
        lassign [choice $possibleStates2] state2
        set occ2 [cphSystem get occupancy $segresidname2 $state2]
        set numProtonsTrial2 [expr {lsum($occ2)}]
        set dnTrial [expr {$numProtonsTrial2 - $numProtonsTrial1}]
        if {$dn == $dnTrial} {
            cphSystem set trialState $segresidname1 $state1
            cphSystem set trialState $segresidname2 $state2
            return 0
        }
    }
    # This is an error and should never happen.
    return -1
}

# ::cphSystem::proposeCotitration
#
# Propose a new trial state requiring cotitration - i.e. two residues
# concurrently changing from protonated to deprotonated or vice versa.
#
# This returns true if no cotitration or else a negative error code.
#
proc ::cphSystem::proposeCotitration {segresidname1 segresidname2} {
    set maxAttempts 10

    set numProtons1 [expr {lsum([cphSystem get occupancy $segresidname1])}]
    set numProtons2 [expr {lsum([cphSystem get occupancy $segresidname2])}]
    set possibleStates1 [cphSystem get trialStateList $segresidname1]
    set possibleStates2 [cphSystem get trialStateList $segresidname2]
    while {true} {
        lassign [choice $possibleStates1] state1
        set occ1 [cphSystem get occupancy $segresidname1 $state1]
        set numProtonsTrial1 [expr {lsum($occ1)}]
        set dn1 [expr {$numProtonsTrial1 - $numProtons1}]
        lassign [choice $possibleStates2] state2
        set occ2 [cphSystem get occupancy $segresidname2 $state2]
        set numProtonsTrial2 [expr {lsum($occ2)}]
        set dn2 [expr {$numProtonsTrial2 - $numProtons2}]
        if {$dn1 == $dn2} {
            cphSystem set trialState $segresidname1 $state1
            cphSystem set trialState $segresidname2 $state2
            return 0
        }
        incr attempts
        if {$attempts >= $maxAttempts} {
            # This probably implies that no cotitration exists.
            return 1
        }
    }
    # This is an error and should never happen.
    return -1
}

# ::cphSystem::proposeResidueTautomerization
#
# Propose a new trial state requiring a tautomerization - i.e. the number of 
# protons remains unchanged. 
# 
# Unlike a state titration, a state tautomerization is not guaranteed to exist 
# for all residues. Return true if one is found, else return false.
#
# In order to ensure that a tautomer is found (if it exists), but also not
# cause an infinite loop when no tautomer exists, the number of state proposals
# is capped at maxAttempts. For a typical residue with a tautomer, there are
# three states, two of which can interconvert via tautomerization. The
# probability of selecting the tautomeric state k times in N trials is thus
# binomially distributed with p = (1 - p) = 0.5:
#
# P(k, N) = [N! / (k!(N-k)!)] 0.5^N
#
# The probability of picking that correct state at least once is thus:
#
# P(k>0, N) = 1 - P(0, N) = 1 - 0.5^N 
#
# or
#
# N = log[1 - P(k>0, N)] / log(0.5)
#
# For P(k>0, N) = 0.999, this gives N ~= 10 (the default).
#
proc ::cphSystem::proposeResidueTautomerization {segresidname} {
    set maxAttempts 10

    set possibleStates [cphSystem get trialStateList $segresidname]
    set numProtons [expr {lsum([cphSystem get occupancy $segresidname])}]
    while {true} {
        lassign [choice $possibleStates] state
        set occ [cphSystem get occupancy $segresidname $state]
        set numProtonsTrial [expr {lsum($occ)}]
        if {$numProtonsTrial == $numProtons} {
            cphSystem set trialState $segresidname $state
            return 1
        }
        incr attempts
        if {$attempts >= $maxAttempts} {
            # This probably implies that no tautomer exists. 
            return 0 
        }
    }
    # This is an error and should never happen.
    return -1
}

# ::cphSystem::computeInherentLogRatio
#
# Compute the (common) log ratio* of inherent (i.e. pH-dependent) weights of
# the given states within the given residue.
#
# Because pKa is strictly pairwise, the relative unnormed weights (Q) are most
# clearly computed as ratios (and as log ratios, for numerical stability).
# There are two cases: 1) tautomerizations and 2) titrations, the former is not
# a proper pKa since the pH dependence cancels (i.e. no protons are entering or
# leaving the bath).
#
# case 1, n = n':
#
#     P(l')
#     ----- = 10^[sgn(l, l') pKa_i(l, l')]
#     P(l)
#
# case 2, n != n':
#     
#     P(l')
#     ----- = 10^{-[sgn(n' - n) pKa_i(l, l') - (n' - n)pH]}
#     P(l)
#
# where l and l' are occupation vectors for the given states and n and n' are
# the number of protons in each state (i.e. the sums of the occupation vector
# elements). By convention, pKa_i(l, l') = pKa_i(l', l) and the antisymmetry
# of adding vs deleting protons is accounted for by the sgn function. Note 
# that, for tautomers, l and l' are converted to indices and the sgn is
# computed by an arbitrary convention.
#
# Arguments:
# ----------
# pH : float
#   The pH value to evaluate the ratio at.
# segresidname : string
#   Residue specification as "<segid>:<resid>:<resname>" 
# statei : string
#   The state whose weight is in the numerator
# statej : string
#   The state whose weight is in the denomenator
#
# Returns:
# --------
# logRatio : float
#   The natural log of the ratio Qj/Qi, i.e. ln(Qj/Qi)
#
proc ::cphSystem::computeInherentLogRatio {pH segresidname statei statej} {
    if {$statei == $statej} {
        # dn = pKai = 0 exactly, but lookup will fail.
        # Note also that, of course, multii = multij.
        return 0.0
    }

    set multii [cphSystem get multiplicity $segresidname $statei]
    set multij [cphSystem get multiplicity $segresidname $statej]
    set li [cphSystem get occupancy $segresidname $statei]
    set lj [cphSystem get occupancy $segresidname $statej]
    set dn [expr {lsum($lj) - lsum($li)}]
    set si [occupancy2Index $li]
    set sj [occupancy2Index $lj]
    set sij [index2flatindex $si $sj]
    set pKai [lindex [cphSystem get pKaiPair $segresidname] $sij]
    if {$dn == 0.0} {
        # tautomerization: "pKa" is positive in direction of lower index
        set sgn [expr {$sj > $si} ? 1.0 : -1.0]
    } else {
        # titration: pKa is positive in direction of fewer protons
        set sgn [expr {$dn > 0} ? 1.0 : -1.0]
    }
    # Note that multii and multij are integers!
    return [expr {$::LN10*($sgn*$pKai - log10((1.*$multij)/$multii) - $dn*$pH)}]
}

# ::cphSystem::computeInherentAcceptance
#
# Compute the (reduced) energy difference for a Monte Carlo move based on the 
# given segresidname, its current and trial state, and the given pH.
#
# The Metropolis criterion can be expressed as:
#
#     P(l --> l') = min{1, Q'/Q} = min{1, e^[-du(l, l')]},
#
# where
#
#     du = -ln(Q'/Q)
#
proc ::cphSystem::computeInherentAcceptance {pH segresidname} {
    set statei [cphSystem get state $segresidname]
    set statej [cphSystem get trialState $segresidname]
    set du [computeInherentLogRatio $pH $segresidname $statei $statej]
    set multij [cphSystem get multiplicity $segresidname $statej]
    # Note the sign flip for compatibility with Metropolis acceptance!
    return [expr {-($du + log($multij))}]
}

# ::cphSystem::computeInherentNormedWeights
#
# Compute the normalized inherent pKa weights (i.e. probability mass function,
# or PMF) of all states within the given residue at the given pH.
#
# The resulting PMF can be directly sampled in order to choose a new state ala
# independence sampling. However, this is often not an efficent scheme for
# Markov chain Monte Carlo and the output can also be used for other approaches
# such as Metropolized independence sampling.
#
proc ::cphSystem::computeInherentNormedWeights {pH segresidname} {
    set currState [cphSystem get state $segresidname]
    set stateList [cphSystem get stateList $segresidname] 

    # Note that all ratios are taken wrt the current state and are thus
    # equivalent to unnormed weights. This would NOT be correct if a different
    # state were in the denomenator each time.
    #
    set logQs [list]
    set logQMax 0.0
    foreach state $stateList {
        set logQ [computeInherentLogRatio $pH $segresidname $currState $state]
        lappend logQs $logQ
        if {$logQ > $logQMax} {
            set logQMax $logQ
        }
    }
    set normedWeights [normalizeLogWeights $logQs $logQMax]
    return [list $stateList $normedWeights]
}

# ::cphSystem::computeSwitchAcceptance
#
# Compute the pairwise energy difference used to shift the work value in the
# Monte Carlo acceptance criteria.
#
#   This is only the configuration _independent_ portion of that difference 
# (i.e. it only depends on the current and trial protonation states). The
# configuration _dependent_ portion is the total work W applied to effect the
# switch. As for the proposal energy, the sign of the correction depends on
# whether this is tautomerization (1) or titration (2) move:
#
# dE(s, s') = c*{-dG(s, s') + kT ln(10) [pKa(s, s') - pKa_i(s, s')]}
# 
# case 1, n = n':  c = sgn(s' - s)
#
# case 2, n != n': c = sgn(n' - n)
#
# Here T is the _reference_ temperature at which pKa is measured and dG is
# computed, not the temperature of the simulation. See 
# computeInherentAcceptance for definition of the remaining notation.
#
proc ::cphSystem::computeSwitchAcceptance {segresidnameList} {
    set du 0.0
    foreach segresidname $segresidnameList {
        set l [cphSystem get occupancy $segresidname]
        set lp [cphSystem get trialOccupancy $segresidname]
        set dn [expr {lsum($lp) - lsum($l)}]
        set s [occupancy2Index $l]
        set sp [occupancy2Index $lp]
        set ssp [index2flatindex $s $sp]
        set dG [lindex [cphSystem get dGPair $segresidname] $ssp]
        set pKa [lindex [cphSystem get pKaPair $segresidname] $ssp]
        set pKai [lindex [cphSystem get pKaiPair $segresidname] $ssp]
        set trialState [cphSystem get trialState $segresidname]
        if {$dn == 0.0} {
            # tautomerization: "pKa" is positive in direction of lower index
            set sgn [expr {$sp > $s} ? 1.0 : -1.0]
        } else {
            # titration: pKa is positive in direction of fewer protons
            set sgn [expr {$dn > 0} ? 1.0 : -1.0]
        }
        set kT [expr {$::BOLTZMANN*[cphSystem get Tref $segresidname]}]
        set du [expr {$du + $sgn*($dG - $kT*$::LN10*($pKa - $pKai))}]
    }
    return $du
}

# ---------------
# Psfgen Routines
# ---------------
# ::cphSystem::alchemifyResidue
#
# Apply a trial alchemical patch to a residue.
#
proc ::cphSystem::alchemifyResidue {segresidname frcCons temp {buildH false}} {
    lassign [cphSystem get alchAtomLists $segresidname] l0atoms l1atoms
    set patch [cphSystem get hybridPatch $segresidname]
    lassign [split $segresidname ":"] segid resid
    # NB: alchPatch uses psfgen style selections!
    alchPatch $patch "$segid:$resid" $l0atoms $l1atoms $frcCons $temp $buildH
    return 0
}

# ::cphSystem::dealchemifyResidue
#
# Remove an alchemical patch from a residue.
#
proc ::cphSystem::dealchemifyResidue {segresidname} {
    lassign [cphSystem get alchAtomLists $segresidname] l0atoms l1atoms
    lassign [split $segresidname ":"] segid resid 
    # NB: alchUnpatch uses psfgen style selections!
    alchUnpatch "$segid:$resid" $l0atoms $l1atoms
    return 0
}

# =============================================================================
# "Constructor" Routines
# =============================================================================
# ::cphSystem::initializeSystem
#
# Initialize the state of the system using:
#
# 1) input data from the user (usually from a restart file)
#
# AND/OR
#
# 2) the ensemble information (i.e. the pH and temperature)
#
proc ::cphSystem::initializeSystem {pH temperature buildH stateInfo} {
    variable ::cphSystem::resDict
    dict for {segresidname resData} $resDict {
        # Assign inherent pKa values.
        if {[dict exists $stateInfo $segresidname pKai]} {
            cphSystem set pKai $segresidname\
                    [dict get $stateInfo $segresidname pKai]
        } else { ;# Default to reference pKa.
            cphSystem set pKai $segresidname [cphSystem get pKa $segresidname]
        }

        # Assign states.
        if {[dict exists $stateInfo $segresidname state]} {
            set state [dict get $stateInfo $segresidname state]
            cphSystem initializeState $segresidname $state 
        } else { ;# default randomization based on pKai and pH
            cphSystem initializeState random $segresidname $pH
        }
    }
    # Final pass - Apply the patches
    foreach segresidname [cphSystem get segresidnames] {
        lassign [split $segresidname ":"] segid resid
        lassign [cphSystem get alchAtomLists $segresidname] atoms
        set masses [list]
        foreach atom $atoms {
            if {[catch {eval segment mass $segid $resid $atom}]} {
                lappend masses -1.0
            } else {
                lappend masses [segment mass $segid $resid $atom]
            }
        }
        # The patch statement reassigns all masses based on RTF definitions.
        # This must be undone in the case of HMR.
        patch [cphSystem get statePatch $segresidname] "$segid:$resid"
        foreach atom $atoms mass $masses {
            if {$mass > 0.0} {
                psfset mass $segid $resid $atom $mass
            } else {
                print "WARNING! Building new atom $segid:$resid:$atom with"\
                      "unknown mass."
                print "         This may cause issues, for example, with"\
                      "mass repartitioning."
            }
        }
    }
    guesscoord
    foreach segresidname [cphSystem get segresidnames] {
        cphSystem alchemifypsf $segresidname 0.0 $temperature $buildH
        cphSystem dealchemifypsf $segresidname
        cphSystem update 1 $segresidname
    }
    regenerate angles dihedrals
    # NB - These changes are only reflected in _memory_ for the cphSystem and
    # psfgen. Nothing has happened to the NAMD PSF/PDB files.
    return
}

# ::cphSystem::buildSystem
#
# Build the residue definitions and residue objects for the system based on 
# the NAMD inputs.
#
proc ::cphSystem::buildSystem {resDefs excludeList} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict [checkResDefs $resDefs]

    # Check for special sub-type residues.
    # In general, THESE WILL NOT BE MATCHED BY NAME, but rather are sub-typed
    # within one or more other residue definitions. Therefore we keep a
    # separate dict of residue names that have a sub-type (these may not even
    # be titratable otherwise, such as terminal amino acids). The keys are
    # the actual residue names and the values are a list of the possible
    # sub-types. However, even if a sub-typed residue is present, that doesn't
    # mean the sub-residue exists. There are two (optional) additional checks
    # in the sub-residue definition:
    #
    # subatoms - a list of atom names that MUST exist
    # notsubatoms - a list of atom names that MUST NOT exist
    #
    # For example, alanine is not titratable, but may contain a titratable
    # C-terminus. The terminal oxygen atoms must exist, but are only titratable
    # if they are not blocked by amination or methylation.
    #
    set subresDefDict [dict create]
    dict for {subresname resDef} $resDefDict {
        if {![dict exists $resDef subtypeof]} continue
        foreach subtypedRes [dict get $resDef subtypeof] {
            if {![dict exists $subresDefDict $subtypedRes]} {
                dict set subresDefDict $subtypedRes [list]
            }
            dict lappend subresDefDict $subtypedRes $subresname
        }
    }

    # Check for residue aliases (i.e. titratable residues whose states may span
    # multiple residue definitions, such as HIS).
    set resAliases [dict create]
    dict for {resname resDef} $resDefDict {
        if {![dict exists $resDef aliases]} continue
        dict set resAliases $resname [dict get $resDef aliases]
    }

    # Read in whatever files were specified to NAMD.
    set Args [list [structure] pdb [coordinates]]
    if {[isset binCoordinates]} {
        lappend Args namdbin [binCoordinates]
        if {[isset binVelocities]} {
            lappend Args velnamdbin [binVelocities]
        }
    }
    resetpsf
    readpsf {*}$Args

    set definedResidues [dict keys $resDefDict]
    set definedSubResidues [dict keys $subresDefDict]
    foreach segid [segment segids] {
        foreach resid [segment resids $segid] {
            set resname [segment residue $segid $resid]
            set segresid [format "%s:%s" $segid $resid]
            set segresidname [format "%s:%s" $segresid $resname]
            # Perform aliasing to fix residues with multiple names.
            dict for {realName altNames} $resAliases {
                if {$resname ni $altNames} {
                    continue
                }
                print "aliasing $segresidname to $realName"
                psfset resname $segid $resid $realName
                set resname $realName
                set segresidname [format "%s:%s" $segresid $resname]
            }

            set resIsDefined [expr {$resname in $definedResidues}]          
            set resIsSubtyped [expr {$resname in $definedSubResidues}]
            set resIsExcluded [expr {$segresidname in $excludeList}]
            # Break here if nothing to do. Be sure to remove successful
            # exclusions from the list so that we can identify bad selections
            # later.
            if {$resIsExcluded} {
                set excludeList [lsearch -inline -all -not $excludeList\
                        $segresidname]
            }
            if {(!$resIsDefined && !$resIsSubtyped) || $resIsExcluded} {
                continue
            }

            # Add the residue to the system.
            if {$resIsDefined} {
                dict set resDict $segresidname [dict create]
                dict set resDict $segresidname resname $resname
            }

            # If the residue is subtype-able, check that subresidues exist.
            # Note that this might be true even if the residue itself is not
            # titratable.
            if {!$resIsSubtyped} continue

            set atoms [segment atoms $segid $resid]
            foreach subresname [dict get $subresDefDict $resname] {
                # Atoms that MUST exist.
                set allSubatomsExist 1
                set subatoms [list]
                if {[dict exists $resDefDict $subresname subatoms]} {
                    set subatoms [dict get $resDefDict $subresname subatoms]
                }
                foreach atom $subatoms { 
                    if {$atom ni $atoms} {
                        set allSubatomsExist 0
                        break
                    }
                }

                # Atoms that must NOT exist.
                set allNotSubatomsDoNotExist 1
                set notsubatoms [list]
                if {[dict exists $resDefDict $subresname notsubatoms]} {
                    set notsubatoms\
                            [dict get $resDefDict $subresname notsubatoms]
                }
                foreach atom $notsubatoms {
                    if {$atom in $atoms} {
                        set allNotSubatomsDoNotExist 0
                        break
                    }
                }

                # Break here if anything doesn't match. 
                if {!$allSubatomsExist || !$allNotSubatomsDoNotExist} continue

                set segresidname [format "%s:%s" $segresid $subresname]
                # After all of that, the residue might be excluded!
                set resIsExcluded [expr {$segresidname in $excludeList}]
                if {$resIsExcluded} {
                    set excludeList [lsearch -inline -all -not $excludeList\
                            $segresidname]
                    continue
                }

                # Add the sub-residue to the system.
                dict set resDict $segresidname [dict create]
                dict set resDict $segresidname resname $subresname
            }
        }
    }
    # This is a bit of a hack way to check that all of the selected exclusions
    # were valid, but it works.
    if {[llength $excludeList]} {
        abort "Cannot exclude non-existant residue(s): $excludeList"
    }

    return $resDict
}

# ::cphSystem::checkResDefs
#
# Check that a dictionary of residue definitions is valid and consistent with
# the topology (RTF) files currently loaded by psfgen.
#
proc ::cphSystem::checkResDefs {resDefs} {
    dict for {resname resDef} $resDefs {
        # Check that all required fields are present.
        foreach dtype [list dG pKa states l0atoms] {
            if {![dict exists $resDef $dtype]} {
                abort "No $dtype entry for residue $resname!"
            }
        }
        set pKa [dict get $resDef pKa]
        set dG [dict get $resDef dG]
        set states [dict get $resDef states]
        set l0atoms [dict get $resDef l0atoms]
        # If not specified, use lazy method of appending a "1". If an l1atoms
        # field is present, check that it matches l0atoms.
        if {![dict exists $resDef l1atoms]} {
            set l1atoms [list]
            foreach atom $l0atoms { 
                lappend l1atoms [format "%s1" $atom]
            }
            dict set resDefs $resname l1atoms $l1atoms
        } else {
            set l1atoms [dict get $resDefs $resname l1atoms]
        }

        if {[llength $l0atoms] != [llength $l1atoms]} { 
            abort "Mismatch in atom definitions for residue $resname"
        }
        # Check that the definitions for pKa and dG are consistent.
        if {[llength $pKa] != [llength $dG]} {
            abort "Mismatch in dG/pKa definitions for residue $resname."
        }
        # Check that the state and occupancy definitions are consistent.
        set numSites [llength [dict get $states [lindex $states 0]]]
        set maxStates [expr {int(pow(2, $numSites))}]
        if {[dict size $states] > $maxStates} {
            abort "Too many states defined for residue $resname!"
        }
        dict for {state occ} $states {
            if {[llength $occ] != $numSites} {
                abort "Bad occupancy definition for $resname state $state."
            }
            # Check that the RTFs define two patches for each state:
            # 1) state patches, with the naming convention: <resname><state>
            # 2) hybrid patches, with the naming convention: <resname>H<state>
            #
            set statePatch [format "%s%s" $resname $state]
            set hybridPatch [format "%sH%s" $resname $state]
            if {$statePatch ni [topology patches]} {
                abort "No patch definition in RTFs for $statePatch!"
            }
            if {$hybridPatch ni [topology patches]} {
                abort "No patch definition in RTFs for $hybridPatch!"
            }
        }
        # Build the pairwise parameters for this definition.
        dict set resDefs $resname dGPair [resDef2Matrix $resDef $dG]
        dict set resDefs $resname pKaPair [resDef2Matrix $resDef $pKa]
    }
    return $resDefs
}

# ::cphSystem::assignState
#
# Assign a protonation state by fiat.
#
#   This is a bit of a hack. The state is assigned and then an arbitrary trial
# state is chosen by rejectionless MC. The states are then swapped and the move
# is accepted as if this had happened in reverse.
#
# Arguments:
# ----------
# segresidname : string
#   residue specification as "<segid>:<resid>:<resname>"
# state : string
#   state to assign 
#
# Returns:
# --------
# None
#
proc ::cphSystem::assignState {segresidname state} {
    cphSystem set state $segresidname $state
    cphSystem propose titration $segresidname
    cphSystem set state $segresidname [cphSystem get trialState $segresidname]
    cphSystem set trialState $segresidname $state
    return 0
}

# ::cphSystem::randomizeState
#
# Assign a protonation state by independence sampling at the given pH.
#
# Arguments:
# ----------
# segresidname : string
#   residue specification as "<segid>:<resid>:<resname>"
# pH : float
#   pH value at which to assign protonation states
#
# Returns:
# --------
# None 
#
proc ::cphSystem::randomizeState {segresidname pH} {
    # Note that computeInherentNormedWeights requires a state to be set for
    # normalization purposes. Therefore choose a random state first and then
    # rigorously sample the distribution.
    set states [cphSystem get stateList $segresidname]
    lassign [choice $states] state
    cphSystem set state $segresidname $state
    lassign [computeInherentNormedWeights $pH $segresidname] states weights
    lassign [choice $states $weights] state
    cphSystem set state $segresidname $state
    # Use the same hack as in assignState.
    cphSystem propose titration $segresidname
    cphSystem set state $segresidname [cphSystem get trialState $segresidname]
    cphSystem set trialState $segresidname $state
    return 0
}

# =============================================================================
# Getter Routines
# =============================================================================
# ::cphSystem::cphSystemGet 
#
# Getters for system and residue attributes, called as:
#
#   <attribute> [<segresidname> [<args>]]
#
# <attribute> is the name of either a system attribute (segresidname selections
# are invalid) or else a residue attribute. For SOME residue attributes, not
# specifying a segresidname will return a list for all residues. Some
# attributes also require some additional arguments (see below)
#
# system attributes  description
# -----------------  -----------
# segresidnames      list of all segresidnames
# numresidues        number of residues in the system
# resdefs            list of defined resnames
# numdefs            number of defined resnames
#
# residue attributes description
# ------------------ -----------
# resname            residue name
# state              current state
# trialState         proposed trial state
# pKai               minimal pKai list for this residue
# dG                 minimal dG list for this residue 
# pKa                minimal pKa list for this residue
# Tref               reference temperature for pKa
# occupancy          occupancy vector for the current state
# trialOccupancy     occupancy vector for the trial state
# stateList          all possible states
# equivStateList     all states equivalent to the given state
# trialStateList     all possible trial states (not the current state)
# dGPair             pair dG for current/trial states
# pKaPair            pair pKa for current/trial states
# pKaiPair           pair pKai for current/trial states 
# statePatch         patch for the current state
# hybridPatch        alchemical patch for the trial state
# alchAtomLists*     lists of alchemical atoms at 0/1
# alchBonds*^        extraBonds entries
#
# * - segresidname selection is required
# ^ - requires additional arguments
#
proc ::cphSystem::cphSystemGet {attr {segresidname {}} args} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict

    set getAll [expr {![llength $segresidname]}]
    if {!$getAll && ![dict exists $resDict $segresidname]} {
        abort "cphSystemGet: Invalid segresidname $segresidname"
    }

    return [switch -nocase -- $attr {
        segresidnames {
            dict keys $resDict
        }
        numresidues {
            dict size $resDict
        }
        resdefs {
            dict keys $resDefDict
        }
        numdefs {
            dict size $resDefDict
        }
        resname -
        state -
        trialState -
        pKai {
            if {$getAll} {
                getAllResAttr $attr
            } else {
                getResAttr $attr $segresidname
            }
        }
        dG -
        pKa -
        Tref {
            if {$getAll} {
                getAllResDefAttr $attr
            } else {
                getResDefAttr $attr $segresidname
            }
        }
        occupancy {
            if {$getAll} {
                getAllOccupancy
            } else {
                if {[llength $args]} {
                    lassign $args state
                    getOccupancy $segresidname $state
                } else {
                    getCurrOccupancy $segresidname
                }
            }
        }
        trialOccupancy {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                getTrialOccupancy $segresidname
            }
        }
        stateList {
            if {$getAll} {
                cannotGetAll $attr 
            } else {
                set resname [getResAttr resname $segresidname]
                dict keys [dict get $resDefDict $resname states]
            }
        }
        equivStateList {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                if {[llength $args]} {
                    lassign $args state
                } else {
                    # Default to current state.
                    set state [getResAttr state $segresidname]
                }
                getEquivStates $segresidname $state
            }
        }
        multiplicity {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                if {[llength $args]} {
                    lassign $args state
                } else {
                    # Default to current state.
                    set state [getResAttr state $segresidname]
                }
                llength [getEquivStates $segresidname $state]
            }
        }
        trialStateList {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                set resname [getResAttr resname $segresidname]
                set state [getResAttr state $segresidname]
                set states [dict keys [dict get $resDefDict $resname states]]
                lsearch -all -inline -not -nocase $states $state
            }
        }
        dGPair -
        pKaPair {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                getResDefAttr $attr $segresidname
            }
        }
        pKaiPair {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                getResAttr $attr $segresidname
            }
        }
        statePatch {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                set resname [getResAttr resname $segresidname]
                set state [getResAttr state $segresidname]
                format "%s%s" $resname $state
            }
        }
        hybridPatch {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                set resname [getResAttr resname $segresidname]
                set trialState [getResAttr trialState $segresidname]
                format "%sH%s" $resname $trialState
            }
        }
        alchAtomLists {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                list [getResDefAttr l0atoms $segresidname]\
                     [getResDefAttr l1atoms $segresidname]
            }
        }
        alchBonds {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                lassign $args k
                set bondEntries [list]
                foreach l0atom [getResDefAttr l0atoms $segresidname]\
                        l1atom [getResDefAttr l1atoms $segresidname] {
                    lassign [split $segresidname ":"] segid resid
                    # Note that atomid indices start at one, not zero!
                    set i [expr {[segment atomid $segid $resid $l0atom] - 1}]
                    set j [expr {[segment atomid $segid $resid $l1atom] - 1}]
                    lappend bondEntries [format "bond %d %d %f %f" $i $j $k 0]
                }
                join $bondEntries "\n"
            }
        }
        default {
            abort "cphSystemGet: Invalid attribute $attr"
        }
    }]
}

proc ::cphSystem::cannotGetAll {attr} {
    abort "cphSystemGet: Cannot get all $attr - must select a segresidname"
    return -1
}

# ------------------------------
# Getters for residue attributes
# ------------------------------
proc ::cphSystem::getResAttr {attr segresidname} {
    variable ::cphSystem::resDict
    return [dict get $resDict $segresidname $attr]
}

proc ::cphSystem::getAllResAttr {attr} {
    variable ::cphSystem::resDict
    set retList [list]
    foreach resData [dict values $resDict] {
        lappend retList [dict get $resData $attr]
    }
    return $retList
}

# -----------------------------------------
# Getters for residue definition attributes
# -----------------------------------------
proc ::cphSystem::getResDefAttr {attr segresidname} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    set resname [dict get $resDict $segresidname resname]
    return [dict get $resDefDict $resname $attr]
}

proc ::cphSystem::getAllResDefAttr {attr} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    set retList [list]
    foreach resData [dict values $resDict] { 
        set resname [dict get $resData resname]
        lappend retList [dict get $resDefDict $resname $attr]
    }
    return $retList
}

# -----------------------------------------------
# Special getters for more complicated operations
# -----------------------------------------------
proc ::cphSystem::state2Occupancy {resname state} {
    return [dict get $::cphSystem::resDefDict $resname states $state]
}

proc ::cphSystem::getOccupancy {segresidname state} {
    variable ::cphSystem::resDict
    set resname [dict get $resDict $segresidname resname]
    return [state2Occupancy $resname $state]
}

proc ::cphSystem::getCurrOccupancy {segresidname} {
    variable ::cphSystem::resDict
    set state [dict get $resDict $segresidname state]
    return [getOccupancy $segresidname $state]
}

proc ::cphSystem::getTrialOccupancy {segresidname} {
    variable ::cphSystem::resDict
    set state [dict get $resDict $segresidname trialState]
    return [getOccupancy $segresidname $state] 
}

proc ::cphSystem::getAllOccupancy {} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    set retList [list]
    foreach resData [dict values $resDict] {
        set resname [dict get $resData resname]
        set state [dict get $resData state]
        set retList [list {*}$retList {*}[state2Occupancy $resname $state]]
    }
    return $retList
}

# ::cphSystem::getEquivStates
#
# Return the list of states equivalent to the given state of the given residue. # This includes the given state such that the length of the list is exactly the
# degeneracy.
#
proc ::cphSystem::getEquivStates {segresidname testState} {
    set l [cphSystem get occupancy $segresidname $testState]
    set s [occupancy2Index $l]

    set equivStateList [list]
    foreach state [cphSystem get stateList $segresidname] {
        if {$state == $testState} {
            lappend equivStateList $state
            continue
        }

        set lp [cphSystem get occupancy $segresidname $state]
        set dn [expr {lsum($lp) - lsum($l)}]
        set sp [occupancy2Index $lp]
        set ssp [index2flatindex $s $sp]
        set pKa [lindex [cphSystem get pKaPair $segresidname] $ssp]
        if {$dn == 0.0 && $pKa == 0.0} {
            lappend equivStateList $state
        }
    }
    return $equivStateList
}

# ----------------------
# Other helper functions
# ----------------------
# Given an occupancy, get the state index.
# This is essentially a left to right conversion from binary to decimal.
# That is, 
#   index = l0*2**0 + l1*2**1 + l2*2**2 + ...
#
proc ::cphSystem::occupancy2Index {occupancy} {
    set index 0
    for {set i 0} {$i < [llength $occupancy]} {incr i} {
        incr index [expr {int([lindex $occupancy $i]*pow(2, $i))}]
    }
    return $index
}

# Map a pair index to the flat off-diagonal index.
#
# This maps to the _lower_ off-diagonal counting left-right, top-bottom.
#
# Ex. n = 4, 2,0 --> 1 and 3,2 --> 5
#
# |x x x x|
# |0 x x x|
# |1 2 x x|
# |3 4 5 x|
#
# The additional assumption here is that element i,j is the same as j,i.
# Therefore this internally swaps the two indices such that i > j _always_.
#
proc ::cphSystem::index2flatindex {i j} {
#    lassign [list [expr {max($i, $j)}] [expr {min($i, $j)}]] I J
    if {$i > $j} {
        set I $i
        set J $j
    } elseif {$i < $j} {
        set I $j
        set J $i
    } else {
        abort "invalid transition from state $i to state $i"
    }
    return [expr {($I*($I - 1) + 2*$J) / 2}]
}

# Convenience function for assignment of flat, off-diagonal, upper-triangular
# matrices using 2d indices. Note that not all index combinations are valid!
#
# Example:
#
# % set myMatrix [lrepeat 3 0.0]
# 0.0 0.0 0.0
# % mset myMatrix 0 1 10.0
# 10.0 0.0 0.0
#
proc ::cphSystem::mset {matrix i j value} {
    upvar 1 $matrix Matrix
    lset Matrix [index2flatindex $i $j] $value
    return $Matrix
}

# Convenience function - inverse of mset.
#
proc ::cphSystem::mindex {matrix i j} {
    return [lindex $matrix [index2flatindex $i $j]]
}

# =============================================================================
# Setter Routines
# =============================================================================
# ::cphSystem::cphSystemSet
#
# Setters for residue attributes, called as:
#
#   <attribute> <segresidname> <value>
#
# <attribute> is the name of a residue attribute.
#
# residue attributes description
# ------------------ -----------
# state              current state
# trialState         proposed trial state
# pKai               minimal pKai list for this residue
#
proc ::cphSystem::cphSystemSet {attr segresidname value} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict

    if {![dict exists $resDict $segresidname]} {
        abort "cphSystemSet: Invalid segresidname $segresidname"
    }

    return [switch -nocase -- $attr {
        state -
        trialState {
            set states [cphSystem get stateList $segresidname]
            if {$value ni $states} {
                if {[string match -nocase $attr trialState]
                    && ![llength $value]} {
                } else {
                    abort "Invalid state assignment $value for residue"\
                            "$segresidname"
                }
            }
            dict set resDict $segresidname $attr $value
            expr {$value}
        }
        pKai {
            set resname [cphSystem get resname $segresidname]
            set resDef [dict get $resDefDict $resname]
            set pKaiMatrix [resDef2Matrix $resDef $value]
            dict set resDict $segresidname pKai $value
            dict set resDict $segresidname pKaiPair $pKaiMatrix
            expr {$pKaiMatrix}
        }
        default {
            abort "cphSystemSet: Invalid attribute $attr"
            expr {-1}
        }
    }]
}

# ::cphSystem::resDef2Matrix
#
# Read a residue definiition for a given thermodynamic attribute (a free energy
# or pKa) and fill in the (flat)matrix elements that are derived from this.
# 
# This permits the configuration files to have "missing" information when it is
# not used. The downside is that essentially every case needs to have its own
# conventions.
#
# Arguments:
# ----------
# resDef : dict
#   The residue definition to be assigned
# data : list
#   The _minimial_ attribute data needed to fill in the matrix
#
# Returns
# -------
#  0 - The matrix was successfully updated.
# -1 - The specific combination of state definitions is not implemented
#      Ex. A state is unexpectedly missing and leaves a gap in the cycle or
#      too many definitions are given
# -2 - There is an inconsistency in the truncated definition
#      Ex. The sum of parameters around the cycle is not 0 - THIS MIGHT SIMPLY
#      IMPLY A NEGATIVE SIGN ERROR!
#
proc ::cphSystem::resDef2Matrix {resDef data} {
    #   The logic and mathematics here is a bit dreadful. Basically, we have
    # networks of states based on the number of protonation sites (numSites).
    # A given network has, at most, maxStates = 2**numSites states, but
    # generally many fewer than this, because some of the states don't exist
    # for physical reasons (e.g. doubly protonated carboxylates).
    #   The problems begin because the residue attributes (dG and pKa) do not
    # describe states, but rather transitions between pairs of states - there
    # are formally numStates**2 pairs! Fortunately, attributes between
    # identical states will always be zero and transitions between pairs of
    # states in opposite directions just have flipped signs. There are
    # therefore only numStates*(numStates - 1)/2 pairs of states that need to
    # be stored (as a flat, off-diagonal matrix, see index2flatindex, note that
    # this conversion is symmetric!). Lookup is always done with respect to the
    # full set of states (i.e. the size is computed with maxStates, not
    # numStates). The row/column indices are computed as a left-right binary
    # conversion of the occupancy vector to a scalar (see occupancy2Index).
    #   Just to make things more complicated, many systems have equivalent
    # sites such that many fewer pairwise values are needed. For example, a
    # carboxylate formally has two sites and four states, one of which is
    # neglected. The two sites are also equivalent, so this system has, at max
    # 2**2*(2**2 - 1)/2 = 6 pairs, but only three are used in practice. Two of
    # these pairs are the same, so only a single value is specified.
    #
    # Long story short: Each case is handled specially and implemented by hand
    # (YUCK!). Fortunately, it is the state pattern that is implemented, not
    # the residue name or some overly specific nonsense. If a pattern has been
    # implemented, then any residue following that pattern works also (e.g.
    # ACE, ASP, and GLU are all one pattern).
    #
    # Final note:
    #   Because only half of the diagonal elements are stored, the signs are
    # computed on-the-fly. This assumes that all proper pKa values are positive
    # - a negative pKa should never be inserted into pKa matrix. For "improper
    # pKas" (i.e. tautomerizations), the assumption is that the value is from
    # the higher index to the lower index (this may be negative).
    #   In the reduced listing all unique transitions from the highest index
    # should be specified first before going to the next highest index for
    # which unique transitions remain. In practice, that explanation is
    # probably pretty useless, so just look at the specific code for the given
    # case when making a new definition.
    #
    set states [dict get $resDef states]
    set numSites [llength [dict get $resDef states [lindex $states 0]]]
    set numStates [dict size $states]
    set maxStates [expr {int(pow(2, $numSites))}]
    set numPairs [expr {int($maxStates*($maxStates - 1) / 2)}]
    # Are any protonation counts missing?
    set protCntExists [lrepeat [expr {$numSites+1}] 0]
    dict for {state occ} $states {
        lset protCntExists [expr {int(lsum($occ))}] 1
    }
    set missingProtCnts [list]
    set i 0
    foreach exists $protCntExists {
        if {!$exists} {
            lappend missingProtCnts $i
        }
        incr i
    }

    set Matrix [lrepeat $numPairs 0.0]
    set numValues [llength $data] 
    set errorCode 0
    switch -- $numSites {
        1 {
            set errorCode [expr {$numStates == 2 ? $errorCode : -1}]
            set errorCode [expr {$numValues == 1 ? $errorCode : -1}]
            lassign $data attr10
            mset Matrix 1 0 $attr10
        } ;# END 1 site, 1 state
        2 {
            switch -- $numStates {
                3 {
                    if {0 in $missingProtCnts} {
                        # state (0,0) is deleted
                        switch -- $numValues {
                            1 { ;# Ex. H2O/H3O+
                                lassign $data attr32
                                mset Matrix 3 2 $attr32
                                mset Matrix 3 1 $attr32
                            }
                            2 { ;# Ex. HIS
                                lassign $data attr32 attr31
                                mset Matrix 3 2 $attr32
                                mset Matrix 3 1 $attr31
                            }
                            default {
                                set errorCode -1
                            }
                        }
                        if {$errorCode} break

                        mset Matrix 2 1 [expr {[mindex $Matrix 3 1]\
                                               - [mindex $Matrix 3 2]}]
                        # END 2 sites, 3 states, no state 0 
                    } elseif {2 in $missingProtCnts} {
                        # state (1,1) is deleted
                        switch -- $numValues {
                            1 { ;# Ex. carboxylates (ASP, GLU)
                                lassign $data attr20
                                mset Matrix 2 0 $attr20
                                mset Matrix 1 0 $attr20
                            }
                            2 { ;# Ex. asymmetric carboxylate?
                                lassign $data attr20 attr10
                                mset Matrix 2 0 $attr20
                                mset Matrix 1 0 $attr10
                            }
                            default {
                                set errorCode -1
                            }
                        }
                        if {$errorCode} break

                        mset Matrix 2 1 [expr {[mindex $Matrix 2 0]\
                                               - [mindex $Matrix 1 0]}]
                        # END 2 sites, 3 states, no state 3
                    } else {
                        # Why would state (0,1) or (1,0) be missing?
                        set errorCode -1
                    }
                } ;# END 2 sites, 3 states
                4 {
                    switch -- $numValues {
                        2 {
                            lassign $data attr32 attr20
                            mset Matrix 3 2 $attr32
                            mset Matrix 3 1 $attr32
                            mset Matrix 2 0 $attr20
                            mset Matrix 1 0 $attr20
                        }
                        4 {
                            lassign $data attr32 attr31 attr20 attr10
                            mset Matrix 3 2 $attr32
                            mset Matrix 3 1 $attr31
                            mset Matrix 2 0 $attr20
                            mset Matrix 1 0 $attr10
                        }
                        default {
                            set errorCode -1
                        }
                    }
                    if {$errorCode} break

                    mset Matrix 3 0 [expr {[mindex $Matrix 3 1]\
                                           + [mindex $Matrix 1 0]}]
                    mset Matrix 2 1 [expr {[mindex $Matrix 2 0]\
                                           - [mindex $Matrix 1 0]}]

                    set alt30 [expr {[mindex $Matrix 3 2]\
                                     + [mindex $Matrix 2 0]}]
                    if {[mindex $Matrix 3 0] != $alt30} {
                        set errorCode -2
                    }
                    # END 2 sites, 4 states, no states missing
                } ;# END 2 sites, 4 states
                default {
                    set errorCode -1
                }
            }
        } ;# END 2 sites
        3 {
            switch -- $numStates {
                4 {
                    if {0 in $missingProtCnts && 1 in $missingProtCnts} {
                        # state (0,0,0) and 1 proton states missing
                        switch -- $numValues {
                            1 { ;# Ex. LYS
                                lassign $data attr76
                                mset Matrix 7 6 $attr76
                                mset Matrix 7 5 $attr76
                                mset Matrix 7 3 $attr76
                            }
                            default {
                                set errorCode -1
                            }
                        }
                        if {$errorCode} break

                        mset Matrix 6 5 [expr {[mindex $Matrix 7 5]\
                                               - [mindex $Matrix 7 6]}]
                        mset Matrix 6 3 [expr {[mindex $Matrix 7 3]\
                                               - [mindex $Matrix 7 6]}]
                        mset Matrix 5 3 [expr {[mindex $Matrix 7 3]\
                                               - [mindex $Matrix 7 5]}]

                        set alt53 [expr {[mindex $Matrix 6 3]\
                                         - [mindex $Matrix 6 5]}]
                        set alt63 [expr {[mindex $Matrix 5 3]\
                                         - [mindex $Matrix 6 5]}]
                        set alt65 [expr {[mindex $Matrix 6 3]\
                                         - [mindex $Matrix 5 3]}]
                        if {[mindex $Matrix 5 3] != $alt53
                            || [mindex $Matrix 6 3] != $alt63
                            || [mindex $Matrix 6 5] != $alt65} {
                            set errorCode -2
                        }
                        # END 3 sites, 4 states, no state 0, 1, 2, 4
                    } elseif {2 in $missingProtCnts && 3 in $missingProtCnts} {
                        # state (1, 1, 1) and 2 proton states missing
                        switch -- $numValues {
                            1 { ;# Ex. phosphate monoesters
                                lassign $data attr40
                                mset Matrix 4 0 $attr40
                                mset Matrix 2 0 $attr40
                                mset Matrix 1 0 $attr40 
                            } 
                            default {
                                set errorCode -1
                            }
                        }
                        if {$errorCode} break
                        
                        mset Matrix 4 2 [expr {[mindex $Matrix 4 0]\
                                               - [mindex $Matrix 2 0]}]
                        mset Matrix 4 1 [expr {[mindex $Matrix 4 0]\
                                               - [mindex $Matrix 1 0]}]
                        mset Matrix 2 1 [expr {[mindex $Matrix 2 0]\
                                               - [mindex $Matrix 1 0]}]
                        # END 3 sites, 4 states, no state 3, 5, 6, 7
                    } else {
                        set errorCode -1
                    }
                } ;# END 3 sites, 4 states
                7 {
                    if {3 in $missingProtCnts} {
                        # state (1, 1, 1) missing
                        switch -- $numValues {
                            2 { ;# Ex. PHOsphate w/o H3PO4
                                lassign $data attr64 attr40
                                mset Matrix 6 4 $attr64
                                mset Matrix 6 2 $attr64
                                mset Matrix 6 1 $attr64
                                mset Matrix 5 4 $attr64
                                mset Matrix 5 2 $attr64
                                mset Matrix 5 1 $attr64
                                mset Matrix 3 4 $attr64
                                mset Matrix 3 2 $attr64
                                mset Matrix 3 1 $attr64

                                mset Matrix 4 0 $attr40
                                mset Matrix 2 0 $attr40
                                mset Matrix 1 0 $attr40
                            }
                            default {
                                set errorCode -1
                            }
                        }
                        if {$errorCode} break

                        mset Matrix 6 0 [expr {[mindex $Matrix 6 4]\
                                               + [mindex $Matrix 4 0]}]
                        mset Matrix 5 0 [expr {[mindex $Matrix 5 4]\
                                               + [mindex $Matrix 4 0]}]
                        mset Matrix 3 0 [expr {[mindex $Matrix 3 1]\
                                               + [mindex $Matrix 1 0]}]
                        mset Matrix 6 5 [expr {[mindex $Matrix 6 4]\
                                               - [mindex $Matrix 5 4]}]
                        mset Matrix 6 3 [expr {[mindex $Matrix 6 2]\
                                               - [mindex $Matrix 3 2]}]
                        mset Matrix 5 3 [expr {[mindex $Matrix 5 1]\
                                               - [mindex $Matrix 3 1]}]
                        mset Matrix 4 2 [expr {[mindex $Matrix 4 0]\
                                               - [mindex $Matrix 2 0]}]
                        mset Matrix 4 1 [expr {[mindex $Matrix 4 0]\
                                               - [mindex $Matrix 1 0]}]
                        mset Matrix 2 1 [expr {[mindex $Matrix 2 0]\
                                               - [mindex $Matrix 1 0]}]
                        # END 3 sites, 7 states, no state 7
                    } else {
                        set errorCode -1
                    }
                } ;# END 3 sites, 7 states
                default {
                    set errorCode -1
                }
            }
        } ;# END 3 sites
        default {
            set errorCode -1
        }
    }

    # Catch errors here.
    switch -- $errorCode {
        -1 {
            abort "Bad or unimplemented specification of $numValues values for\
                    $numSites site residue with $numStates states"
        }
        -2 {
            abort "Error in thermodynamic cycle"
        }
    }

    return $Matrix
}

