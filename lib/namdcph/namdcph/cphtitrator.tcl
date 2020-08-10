# cphtitrator.tcl
#
#   This file provides the ::cphTitrator namespace, which effectively emulates
# an object containing all constant-pH information and routines relevant to
# proposing and performing Monte Carlo titrations.
#
#   If cphSystem and cphTitrator were true objects, this would essentially be
# an object composition. That is, a cphTitrator "contains" a cphSystem and
# yields information pertaining to it depending on a context (e.g. the pH).
#
package require Tcl 8.5

source [file join [file dirname [info script]] "cphtoppar.tcl"]
source [file join [file dirname [info script]] "numtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::cphTitrator {
    #   The core information of the cphTitrator is the moveDict, which stores
    # all possible neMD/MC moves involving one or more titratable residues. The
    # keys are of the form:
    #
    # segresidname[/segresidname..[/segresidname]]
    #
    # in analogy to the keys of ::cphSystem::residueDict. Note that a "/" is
    # used as a syntactic divider between residues. Only _specific_ move
    # information is kept in moveDict. Another dict, defaultMoveParams is used
    # to store default information (similar to ::cphSystem::residueDefDict) and
    # is used to avoid saving copies of redundant information.
    #
    namespace import ::cphSystem::*

    variable moveDict [dict create]
    variable defaultMoveParams [dict create]
    variable maxAttempts

    namespace export cphSystem
    namespace export cphTitrator
}

# =============================================================================
# cphTitrator Interface function
# =============================================================================
# ::cphTitrator::cphTitrator
#
# This is the only exported function from the cphTitrator namespace and
# provides a complete interface. The first argument is always an action (e.g.
# "get" or "set") and the rest are variable arguments.
#
# action           description
# ---------------- -----------
# get              see cphTitratorGet 
# set              see cphTitratorSet
# build            see buildTitrator
# propose          select and propose a move, return switch information
# accumulateStats  accumulate mean acceptance rate for the given move
#
# Some notes about Tcl switch statements for non-guru types like myself:
#
#   These can appear a bit delicate to those who are not familiar with their
# quirks (compared to C). For one, badly placed comments can break them. Also
# note that a switch returns the expression it first matches and does NOT use
# break statements.
#
proc ::cphTitrator::cphTitrator {action args} {
    return [switch -nocase -- $action {
        get {
            cphTitratorGet {*}$args
        }
        set {
            cphTitratorSet {*}$args
        }
        build {
            buildTitrator {*}$args 
        }
        propose { ;# sample the move set and return the selected move info.
            proposeMove {*}$args
        }
        accumulateStats { ;# accumulate MC statistics
            accumulateAcceptanceRate {*}$args
        }
        default {
            abort "Invalid cphTitrator action $action"
        }
    }]
}

# ::cphTitrator::validateMoveLabel
#
# Check that a given move label is valid. That is:
#
# 1) is it in the correct segresidname[/segresidname ...] format?
# 2) are the component segresidnames valid?
#
# A non-zero error code and msg are returned for each condition. Zero is
# returned for a valid segresidname.
#
# NB! This should only be called _after_ buildSystem.
#
proc ::cphTitrator::validateMoveLabel {moveLabel} {
    set segresidnameList [split $moveLabel "/"]
    if {![llength $segresidnameList]} {
        print "Move labels should be of the form segresidname/segresidname ...]"
        return -1
    }
    foreach segresidname $segresidnameList { 
        if {[cphSystem validate $segresidname]} {
            print "Bad segresidname $segresidname in move label $moveLabel!"
            return -2
        }
    }
    return 0
}

# Pop the value from a dict with the corresponding key. If no key exists,
# return the the default value instead.
#
proc ::cphTitrator::dictPopOrDefault {myDict key {defaultValue {}}} {
    upvar 1 $myDict MyDict
    if {[dict exists $MyDict $key]} {
        set value [dict get $MyDict $key]
        dict unset MyDict $key
    }
    if {![info exists value]} {
        set value $defaultValue
    }
    return $value
}

# =============================================================================
# Proposal Routines
# =============================================================================
# ::cphTitrator::proposeMove
#
# Propose a move from all possible moves in the set. This is done by directly
# sampling the probability mass function of the weighted moves. Once this is
# done, test the inherent probability of that move given the imposed pH. If
# rejected, try again up to the given number of "maxAttempts". Always return
# the information necessary to perform a switch for the given proposal.
#
proc ::cphTitrator::proposeMove {pH} {
    set weights [cphTitrator get weight]
    lassign [choice [cphTitrator get moveLabels] $weights] moveLabel
    set proposalCmd [cphTitrator get proposalCmd $moveLabel]
    set accept [eval $proposalCmd $pH]
    set segresidnameList [cphTitrator get segresidnameList $moveLabel]
    set numsteps [cphTitrator get numsteps $moveLabel]
    return [list $accept $numsteps $segresidnameList]
}

# FOR TESTING ONLY - INDEPENDENCE SAMPLING
#
#proc ::cphTitrator::proposeResidueMove {segresidname pH} {
#    lassign [cphSystem compute inherentWeights $pH $segresidname] states weights
#    lassign [choice $states $weights] trialState
#    set currState [cphSystem get state $segresidname]
#    if {$trialState == $currState} {
#        return 0
#    } else {
#        cphSystem set trialState $segresidname $trialState
#        return 1
#    }
#}

# FOR TESTING ONLY - METROPOLIS MONTE CARLO (NON-ERGODIC IN SOME CASES!)
#
#proc ::cphTitrator::proposeResidueMove {segresidname pH} {
#    cphSystem propose titration $segresidname
#    set du [cphSystem compute inherent $pH $segresidname]
#    return [metropolisAcceptance $du]
#}

# ::cphTitrator::proposeResidueMove
#
# Propose a move involving a single residue via Metropolized independence
# sampling. Return True if such a move was proposed, accepted, and stored,
# otherwise return False.
#
proc ::cphTitrator::proposeResidueMove {segresidname pH} {
    lassign [cphSystem compute inherentWeights $pH $segresidname] states weights
    # Construct the Metropolized proposal weights for i-->j.
    set equivStates [cphSystem get equivStateList $segresidname]
    set proposalWeights [list]
    set stateiProposalNorm 1.0
    foreach state $states weight $weights {
        if {$state in $equivStates} {
            lappend proposalWeights 0.0
            set stateiProposalNorm [expr {$stateiProposalNorm - $weight}]
        } else {
            lappend proposalWeights $weight
        }
    }
    lassign [choice $states $proposalWeights] trialState
    # Construct the Metropolized proposal weights for j-->i.
    set equivStates [cphSystem get equivStateList $segresidname $trialState]
    set proposalWeights [list]
    set statejProposalNorm 1.0
    foreach state $states weight $weights {
        if {$state in $equivStates} {
            lappend proposalWeights 0.0
            set statejProposalNorm [expr {$statejProposalNorm - $weight}]
        } else {
            lappend proposalWeights $weight
        }
    }
    # NB: Invert log argument bc of negative sign in Metropolis routine.
    set du [expr {log($statejProposalNorm / $stateiProposalNorm)}]
    set accept [metropolisAcceptance $du]
    if {$accept} {
        cphSystem set trialState $segresidname $trialState
    }
    return $accept
}

# ::cphTitrator::proposeTwoResidueMove
#
# Propose a move involving two residues. This may take two forms:
#
# 1) a proton transfer, in which a proton is moved from one resiue to another
#
# 2) a co-titration, in which two residues change concurrently
#
# Return True if either kind of move is proposed, accepted, and stored.
#
proc ::cphTitrator::proposeTwoResidueMove {moveLabel pH} {
    lassign [split $moveLabel "/"] srn1 srn2 ;# srn = segresidname

    set doPT [cphTitrator get protonTransfer $moveLabel]
    set doCT [cphTitrator get cotitration $moveLabel]
    if {$doPT} {
        set doPT [expr {![cphSystem propose protonTransfer $srn1 $srn2]}]
    }
    # We don't need to bother with co-titration if proton transfer is already
    # proposed, since that means it isn't even possible. The order of these
    # could equally be switched, but it is generally faster to check that
    # proton transfer is possible, so this is done first.
    #
    if {!$doPT && $doCT} {
        set doCT [expr {![cphSystem propose cotitration $srn1 $srn2]}]
    }
    if {$doPT || $doCT} {
        set du1 [cphSystem compute inherent $pH $srn1]
        set du2 [cphSystem compute inherent $pH $srn2]
        return [metropolisAcceptance [expr {$du1 + $du2}]]
    }
    return 0
}

# =============================================================================
# "Constructor" Routines 
# =============================================================================
# ::cphTitrator::buildTitrator
#
# Construct the titrator given a set of MC move information.
#
proc ::cphTitrator::buildTitrator {moveInfo} {
    variable ::cphTitrator::moveDict
    variable ::cphTitrator::defaultMoveParams
    variable ::cphTitrator::maxAttempts

    # 1) Pop off global and default move parameters.
    #
    set maxAttempts [dictPopOrDefault moveInfo maxProposalAttempts 0]
    set maxAttempts 1
    set defaultMoveParams [dict merge $defaultMoveParams\
                                      [dictPopOrDefault moveInfo default]]
    foreach attr {weight protonTransfer cotitration}\
            hardcodedDefault {1.0 0 0} {
        if {![dict exists $defaultMoveParams $attr]} {
            cphTitrator set $attr default $hardcodedDefault
        }
    }

    # 2) Validate all remaining keys - these better be move labels!
    #
    dict for {moveLabel data} $moveInfo {
        if {[validateMoveLabel $moveLabel]} {
            abort
        }
    }

    # The default move set is to titrate each residue independently.
    foreach segresidname [cphSystem get segresidnames] {
        set moveLabel $segresidname
        cphTitrator set proposalCmd $moveLabel\
                "proposeResidueMove $segresidname"
        if {![dict exists $moveInfo $moveLabel]} continue

        dict for {attr value} [dict get $moveInfo $moveLabel] {
            cphTitrator set $attr $moveLabel $value
        }
        dict unset moveInfo $moveLabel
    }
    # Build proton transfer and co-titration moves, if present.
    dict for {moveLabel data} $moveInfo {
        set isTwoResidueMove 0
        foreach type {protonTransfer cotitration} {
            if {[dict exists $data $type] && [dict get $data $type]} {
                set isTwoResidueMove 1
                break
            }
        }
        if {!$isTwoResidueMove} continue
        if {[llength [split $moveLabel "/"]] != 2} {
            abort "Proton transfer and co-titration moves must have exactly 2\
                   segresidnames!"
        }
        cphTitrator set proposalCmd $moveLabel\
                "proposeTwoResidueMove $moveLabel"
        dict for {attr value} $data {
            cphTitrator set $attr $moveLabel $value
        }
        dict unset moveInfo $moveLabel
    }
    # Abort if extra invalid move information was provided.
    if {[dict size $moveInfo]} {
        print "ERROR! Invalid move specifications:"
        dict for {moveLabel data} $moveInfo {
            print "$moveLabel $data"
        }
        abort
    }
    # Initialize statistics for any moves that are new.
    dict for {moveLabel data} $moveDict {
        if {![dict exists $data attempts]} {
            cphTitrator set successes $moveLabel 0
            cphTitrator set attempts $moveLabel 0
        }
    }

    return
}

# =============================================================================
# Getter Routines
# =============================================================================
# ::Titrator::cphTitratorGet
#
# Getter for move attributes, called as:
#
#   <attribute> [<moveLabel>]
#
# With some specialized exceptions, <attribute> is either the name of a system
# attribute (usually a list) or else a move attribute.
#
# system attributes  description
# -----------------  -----------
# moveLabels         list of all the moveLabels
# maxAttempts        max # of attempts at move proposals
#
# move attributes  description
# ---------------  -----------
# numsteps         number of steps in a switch after successful proposal
# weight           weight for consideration in move selection
# protonTransfer   if 2 residues are involved, can they proton transfer?
# cotitration      if 2 residues are involved, can they cotitrate?
# successes        the number of successful _switches_ for this move
# attempts         the number of attempted _switches_ for this move
# proposalCmd      the command needed to set up trial states
# segresidnameList list of the residues involved (NB: this is just shorthand
#                  for [split $moveLabel "/"])
#
# specialized calls  description
# -----------------  -----------
# nodefaults         Return a dictionary of all current move information, but
#                    suppress entries where a move has the default value. This
#                    gives a minimal, but full, description of the move set.
#
# For now this is _much_ simpler than the analogous cphSystemGet. In the future
# it may be useful to generalize this though.
#
proc ::cphTitrator::cphTitratorGet {attr {moveLabel {}}} {
    variable ::cphTitrator::moveDict    
    variable ::cphTitrator::maxAttempts

    set getAll [expr {![llength $moveLabel]}]
    if {!$getAll && ![dict exists $moveDict $moveLabel]
        && ![string match -nocase $moveLabel default]} {
        abort "cphTitratorGet: Invalid moveLabel $moveLabel"
    }

    return [switch -nocase -- $attr {
        moveLabels {
            dict keys $moveDict
        }
        maxAttempts {
            expr {$maxAttempts}
        }
        numsteps -
        weight -
        protonTransfer -
        cotitration -
        successes -
        attempts -
        proposalCmd {
            if {$getAll} {
                getAllMoveAttr $attr 
            } else {
                getMoveAttr $attr $moveLabel
            }
        }
        segresidnameList {
            if {$getAll} {
                getAllSegresidnameLists
            } else {
                split $moveLabel "/"
            }
        }
        archive {
            getMinimalMoveDict
        }
        default {
            abort "cphTitratorGet: Invalid attribute $attr"
        }
    }]
}

proc ::cphTitrator::cannotGetAll {attr} {
    abort "cphTitratorGet: Cannot get all $attr - must select a moveLabel"
    return -1
}

# ---------------------------
# Getters for move attributes
# ---------------------------
# Return the default value if no specific value was set.
proc ::cphTitrator::getMoveAttr {attr moveLabel} {
    variable ::cphTitrator::moveDict
    variable ::cphTitrator::defaultMoveParams

    if {[dict exists $moveDict $moveLabel $attr]} {
        return [dict get $moveDict $moveLabel $attr]
    } elseif {[dict exists $defaultMoveParams $attr]} {
        return [dict get $defaultMoveParams $attr]
    }
    abort "cphTitratorGet: Error getting $attr for move $moveLabel"
    return -1
}

proc ::cphTitrator::getAllMoveAttr {attr} {
    set retList [list]
    foreach moveLabel [cphTitrator get moveLabels] {
        lappend retList [cphTitrator get $attr $moveLabel]
    }
    return $retList
}

proc ::cphTitrator::getAllSegresidnameLists {} {
    set retList [list]
    foreach moveLabel [cphTitrator get moveLabels] {
        lappend retList [split $moveLabel "/"]
    }
    return $retList
}

# ----------------------------
# Special getter for archiving
# ----------------------------
# Return a minimal version of moveDict. This removes inferred values such as
# command names as well as specific settings that now match the defaults
# (since these can change after a restart).
#
proc ::cphTitrator::getMinimalMoveDict {} {
    variable ::cphTitrator::defaultMoveParams

    set minMoveDict [dict create]
    dict set minMoveDict maxProposalAttempts [cphTitrator get maxAttempts]
    dict set minMoveDict default $defaultMoveParams
    foreach moveLabel [cphTitrator get moveLabels] {
        foreach attr [dict keys $defaultMoveParams] {
            set defaultValue [cphTitrator get $attr default]
            set thisValue [cphTitrator get $attr $moveLabel]
            if {$thisValue != $defaultValue} {
                dict set minMoveDict $moveLabel $attr $thisValue
            }
        }
        dict set minMoveDict $moveLabel attempts\
                [cphTitrator get attempts $moveLabel]
        dict set minMoveDict $moveLabel successes\
                [cphTitrator get successes $moveLabel]
    }
    return $minMoveDict
}

# =============================================================================
# Setter Routines
# =============================================================================
# ::cphTitrator::cphTitratorSet
#
# Setters for move attributes, called as:
#
#  <attribute> <moveLabel> <value>
#
# <attribute> is the name of a move attribute.
#
# move attributes  description
# ---------------  -----------
# numsteps         number of steps in a switch after successful proposal
# weight           weight for consideration in move selection
# protonTransfer   if 2 residues are involved, can they proton transfer?
# cotitration      if 2 residues are involved, can they cotitrate?
# successes        the number of successful _switches_ for this move
# attempts         the number of attempted _switches_ for this move
# proposalCmd      the command needed to set up trial states
#
proc ::cphTitrator::cphTitratorSet {attr moveLabel value} {
    variable ::cphTitrator::moveDict
    variable ::cphTitrator::defaultMoveParams

    if {[string match -nocase $moveLabel default]} {
        set setDefault 1
    } else {
        if {[validateMoveLabel $moveLabel]} {
            abort
        }
        set setDefault 0
    }

    return [switch -nocase -- $attr {
        numsteps -
        successes -
        attempts -
        protonTransfer -
        cotitration { ;# Require integer or boolean argument.
            set value [expr {int($value)}]
            if {$setDefault} {
                dict set defaultMoveParams $attr $value
            } else {
                dict set moveDict $moveLabel $attr $value
            }
            expr {$value}
        }
        weight { ;# Require float argument.
            set value [expr {1.*$value}]
            if {$setDefault} {
                dict set defaultMoveParams $attr $value
            } else {
                dict set moveDict $moveLabel $attr $value
            }
            expr {$value}
        }
        proposalCmd { ;# Argument is string or pure list.
            dict set moveDict $moveLabel $attr $value
            expr {$value}
        }
        default {
            abort "cphTitratorSet: Invalid attribute $attr"
            expr {-1}
        }
    }]
}

# =============================================================================
# Routines for tracking and reporting MC statistics
# =============================================================================
# ::cphTitrator::accumulateAcceptanceRate
#
# Accumulate statistics for the given move.
#
proc ::cphTitrator::accumulateAcceptanceRate {accept moveLabel} {
    # Alas, dict incr does not support nested keys.
    set successes [cphTitrator get successes $moveLabel]
    set attempts [cphTitrator get attempts $moveLabel]
    incr successes $accept
    incr attempts
    cphTitrator set successes $moveLabel $successes
    cphTitrator set attempts $moveLabel $attempts
    return
}

