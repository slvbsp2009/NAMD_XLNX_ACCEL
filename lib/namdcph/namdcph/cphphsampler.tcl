# cphphsampler.tcl
#
# Routines for expanded ensemble sampling of pH values.
#
# THESE ARE EXPERIMENTAL ROUTINES THAT ARE NOT FULLY INTEGRATED INTO THE CODE.
# In order to use these routines one would have to modify cphRun to access the
# cphpHSampler namespace. Early tests showed some promise but ultimately did
# not work as hoped or without considerable calibration.
#
source [file join [file dirname [info script]] "namdtcl.tcl"]
source [file join [file dirname [info script]] "numtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::cphpHSampler {
    variable pHLadder [list]
    variable pHWeights [list]
    variable pHIndex -1
    variable deltapHIndex
    variable gibbsSamplingCmd

    namespace export cphpHSampler
}

proc ::cphpHSampler::cphpHSampler {action args} {
    return [switch -nocase -- $action {
        get {
            cphpHSamplerGet {*}$args
        }
        set {
            cphpHSamplerSet {*}$args
        }
        build {
            buildpHSampler {*}$args
        }
        sample {
            sampleNewpH {*}$args
        }
    }]
}

# =============================================================================
# Gibbs Sampling Routines
# =============================================================================
proc ::cphpHSampler::sampleNewpH {numProtons} {
    variable ::cphpHSampler::pHIndex
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::gibbsSamplingCmd
 
    lassign [$gibbsSamplingCmd $numProtons] j accept du
    set pHi [lindex $pHLadder $pHIndex]
    set pHj [lindex $pHLadder $j]
    print [format "namdcph) pH change attempt % 5.2f --> % 5.2f : du = % 11.4f : accept %d" $pHi $pHj $du $accept] 
    if {$accept} {
        set pHIndex $j
    }
    return [lindex $pHLadder $pHIndex]
}

proc ::cphpHSampler::metropolizedConvectiveSampling {numProtons} {
    variable ::cphpHSampler::pHIndex
    variable ::cphpHSampler::deltapHIndex
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::pHWeights

    # Compute normalized weights for all states.
    set logqs [list]
    set logqMax {}
    foreach pH $pHLadder weight $pHWeights {
        lappend logqs [expr {$weight - $::LN10*$pH*$numProtons}]
        if {[lindex $logqs end] > $logqMax} {
            set logqMax [lindex $logqs end]
        }
    }
    # Correct for the fact that deltapHIndex has only ONE possible value at
    # i/j = 1, M and TWO possible values for all other i/j. Note that this is
    # the same as decrementing the non-endpoints by ln(2) due to normalization. 
    lincr logqs 0 $::LN2
    lincr logqs end $::LN2
    set normedWeights [normalizeLogWeights $logqs $logqMax]
    # Compute the marginal weight in the convective direction. If a move in the
    # opposite direction is overwhelmingly favored, then stay put. Otherwise
    # choose a new state and accept/reject with a modified Metropolis
    # criterion.
    #
    set i $pHIndex
    set numpHs [llength $normedWeights]
    set choiceWeights [lrepeat $numpHs 0.0]
    set pic 0.0
    if {$deltapHIndex > 0} {
        set ip1 [expr {$i + 1}]
        for {set k $ip1} {$k < $numpHs} {incr k} {
            lset choiceWeights $k [lindex $normedWeights $k]
            set pic [expr {$pic + [lindex $normedWeights $k]}]
        }
    } else { ;# $deltapHIndex < 0
        for {set k 0} {$k < $i} {incr k} {
            lset choiceWeights $k [lindex $normedWeights $k]
            set pic [expr {$pic + [lindex $normedWeights $k]}]
        }
    }
    # choiceWeights is now a copy of normedWeights with all weights in the
    # non-convective direction assigned to be zero. pic is now the sum of all
    # weights in the convective direction.
    #
    if {$pic < 1e-3} {
        set j $i
        set du 0.0
        set accept 0
    } else {
        lassign [choice $normedWeights $choiceWeights] pj j
        # Consider the same procedure for the oppositive convective direction
        # starting at j.
        if {$deltapHIndex > 0} {
            set jm1 [expr {$j - 1}]
            set pjc [expr {lsum([lrange $normedWeights 0 $jm1])}]
        } else { ;# $deltapHIndex < 0
            set jp1 [expr {$j + 1}]
            set pjc [expr {lsum([lrange $normedWeights $jp1 end])}]
        }
        set du [expr {log($pjc / $pic)}]
        set accept [metropolisAcceptance $du]
        set numpHsm1 [expr {$numpHs - 1}]
        if {$accept && ($j == 0 || $j == $numpHsm1)} {
            set deltapHIndex [expr {-1*$deltapHIndex}]
        } 
    }
    return [list $j $accept $du]
}

proc ::cphpHSampler::metropolizedIndependenceSampling {numProtons} {
    variable ::cphpHSampler::pHIndex
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::pHWeights

    # Compute normalized weights for all states.
    set logqs [list]
    set logqMax {}
    foreach pH $pHLadder weight $pHWeights {
        lappend logqs [expr {$weight - $::LN10*$pH*$numProtons}]
        if {[lindex $logqs end] > $logqMax} {
            set logqMax [lindex $logqs end]
        }
    }
    set normedWeights [normalizeLogWeights $logqs $logqMax]
    # Test that the current state is not overwhelmingly stable. If not, choose
    # another state state and accept/reject based on a modified Metropolis
    # criterion.
    set i $pHIndex
    set pic [expr {1. - [lindex $normedWeights $i]}]
    if {$pic < 1e-3} {
        set j $i
        set du 0.0
        set accept 0
    } else {
        lset normedWeights $i 0.0
        lassign [choice $normedWeights $normedWeights] pj j
        set du [expr {log((1. - $pj) / $pic)}]
        set accept [metropolisAcceptance $du]
    }
    return [list $j $accept $du]
}

proc ::cphpHSampler::convectiveNeighborSampling {numProtons} {
    variable ::cphpHSampler::pHIndex
    variable ::cphpHSampler::deltapHIndex
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::pHWeights

    # Propose move from index i to index j
    set i $pHIndex
    set j [expr {$i + $deltapHIndex}]
    set numpHsm1 [expr {[llength $pHLadder] - 1}]

    set dweight [expr {[lindex $pHWeights $j] - [lindex $pHWeights $i]}]
    set dpH [expr {[lindex $pHLadder $j] - [lindex $pHLadder $i]}]
    set du [expr {$::LN10*$dpH*$numProtons - $dweight}]
    # Correct for the fact that deltapHIndex has only ONE possible value at
    # i/j = 1, M and TWO possible values for all other i/j.
    if {0 < $i && $i < $numpHsm1} {
        set du [expr {$du - $::LN2}]
    }
    if {0 < $j && $j < $numpHsm1} {
        set du [expr {$du + $::LN2}]
    }
    set accept [metropolisAcceptance $du]
    if {$accept && ($j == 0 || $j == $numpHsm1)} {
        set deltapHIndex [expr {-1*$deltapHIndex}]
    }
    return [list $j $accept $du]
}

proc ::cphpHSampler::oscillatorySampling {numProtons} {
    variable ::cphpHSampler::pHIndex
    variable ::cphpHSampler::deltapHIndex
    variable ::cphpHSampler::pHLadder

    # Propose move from index i to index j
    set i $pHIndex
    set j [expr {$i + $deltapHIndex}]
    set numpHsm1 [expr {[llength $pHLadder] - 1}]
    if {0 < $i && $i < $numpHsm1} {
        # Accept 100% of up/down moves within the ladder.
        set accept 1
    } else {
        # Accept 50% of up/down moves from the endpoints.
        set accept [expr {round(rand())}]
    }
    if {$accept && ($j == 0 || $j == $numpHsm1)} {
        set deltapHIndex [expr {-1*$deltapHIndex}]
    } 
    return [list $j $accept 0.0]
}

# =============================================================================
# "Constructor" Routines 
# =============================================================================
# ::cphpHSampler::buildpHSampler
proc ::cphpHSampler::buildpHSampler {pH samplerInfo} {
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::deltapHIndex
    variable ::cphpHSampler::pHWeights
    variable ::cphpHSampler::gibbsSamplingCmd

    if {![dict exists $samplerInfo ladder]} {
        abort "Cannot construct pHSampler without pH ladder info!"
    }
    set ladder [dict get $samplerInfo ladder]
    if {![dict exists $samplerInfo weights]} {
        dict set samplerInfo weights [lrepeat [llength $ladder] 0.0]
    }
    set weights [dict get $samplerInfo weights]
    if {[llength $ladder] != [llength $weights]} {
        abort "Mismatch in size of pHSampler ladder and weights!"
    }
    set pHLadder $ladder
    set pHWeights $weights
    set numpHs [llength $pHLadder]

    if {![dict exists $samplerInfo gibbsMethodName]} {
        set gibbsMethodName "convective-neighbor"
    } else {
        set gibbsMethodName [dict get $samplerInfo gibbsMethodName]
    }
    switch -- [string tolower $gibbsMethodName] {
        convective-neighbor {
            set gibbsSamplingCmd convectiveNeighborSampling
        }
        metropolized-independence {
            set gibbsSamplingCmd metropolizedIndependenceSampling
        }
        metropolized-convective {
            set gibbsSamplingCmd metropolizedConvectiveSampling
        }
        oscillatory {
            set gibbsSamplingCmd oscillatorySampling
        }
        default {
            abort "Unknown Gibbs sampling procedure $gibbsMethodName"
        }
    }

    # Find the initial index.
    set i 0
    foreach pHValue $pHLadder {
        if {$pH == $pHValue} {
            cphpHSampler set index $i
            break
        }
        incr i
    }
    if {$i >= $numpHs} {
        abort "Cannot find pH value $pH in the ladder"
    }
    if {$i <= [expr {$numpHs / 2}]} {
        set deltapHIndex 1
    } else {
        set deltapHIndex -1
    }
    return
}

# =============================================================================
# Getter Routines
# =============================================================================
# ::cphpHSampler::cphpHSamplerGet
proc ::cphpHSampler::cphpHSamplerGet {attr} {
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::pHWeights
    variable ::cphpHSampler::pHIndex

    return [switch -nocase -- $attr {
        ladder {
            expr {$pHLadder}
        }
        weights {
            expr {$pHWeights}
        }
        default {
            abort "cphpHSamplerGet: Invalid attribute $attr"
        }
    }]
}

# =============================================================================
# Setter Routines
# =============================================================================
# ::cphpHSampler::cphpHSamplerSet
proc ::cphpHSampler::cphpHSamplerSet {attr value} {
    variable ::cphpHSampler::pHLadder
    variable ::cphpHSampler::pHIndex

    switch -nocase -- $attr {
        index {
            set pHIndex [expr {int($value)}]
            if {$pHIndex < 0 || [llength $pHLadder] <= $pHIndex} {
                abort "cphpHSamplerSet: Invalid pHIndex $pHIndex"
            }
        }
    }
}

