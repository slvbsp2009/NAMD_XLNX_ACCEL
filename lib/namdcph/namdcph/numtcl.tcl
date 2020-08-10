# numtcl.tcl
#
# Simple routines for operations on purely numerical lists. This is an 
# extremely weak attempt to regain some of my favorite functionality from numpy
# but with more Tcl-ish syntax.
#
# WARNING: USE AT YOUR OWN RISK, THERE IS GENERALLY VERY LITTLE ERROR CHECKING
# BEING PERFORMED. IF USED IMPROPERLY THESE ARE LIKELY TO FAIL WITH RATHER
# OPAQUE ERROR MESSAGES.
#

# tcl::mathfunc::lsum
#
# Sum the elements in a list.
#
proc tcl::mathfunc::lsum {list} {
    set sum 0.0
    foreach item $list {set sum [expr {$sum + $item}]}
    expr {$sum}
}

# tcl::mathfunc::lmean
#
# Average the elements in a list.
#
proc tcl::mathfunc::lmean {list} {
    set n [llength $list]
    set sum 0.0
    foreach item $list {set sum [expr {$sum + $item}]}
    expr {$sum/$n}
}

# lincr
#
# Increment a list (in place) by the given value. The modified list is also
# returned.
#
# - This supports the standard lindex syntax for nested lists (ala lindex).
# - NumPy style "broadcasting" is also supported, whereby a list plus a value
#   results in each element of the list being incremented and a list plus a
#   list results in element by element addition.
#   NB: This only works left to right (ala +=) and modifies the left list.
#
# Example:
#
# % set foo {0 2 1 7 3}
# % lincr foo 2
# 2 4 3 9 5
# % lincr foo {0 1 2 3 4}
# 2 5 5 12 9
#
# % set bar {0 1 {2 3 4} 5}
# % lincr bar 2 1
# 0 1 {3 4 5} 5
# % lincr bar 2 1 1
# 0 1 {3 5 5} 5
# 
proc lincr {list args} {
    upvar 1 $list List
    set value [lindex $args end]
    set index [lrange $args 0 end-1]
    
    set target [lindex $List {*}$index]
    set vlen [llength $value]
    set tlen [llength $target]

    if {$vlen == 1} {
        if {$tlen == 1} {
            # increment a single element by a single value
            lset List {*}$index [expr {$target + $value}]
        } elseif {$tlen > 1} {
            # increment a each element in a list by a single value
            set i 0 
            foreach item $target {
                if {[llength $item] != 1} {
                    error "lincr encountered a nested list during list + value\
                            operation"
                }
                lset List {*}$index $i [expr {$item + $value}]
                incr i
            }
        }
    } elseif {$vlen > 1} {
        if {$tlen != $vlen} {
            error "lincr cannot perform list + list operation with different\
                    lengths ($tlen != $vlen)"
        }
        set i 0
        foreach item1 $target item2 $value {
            lset List {*}$index $i [expr {$item1 + $item2}]
            incr i
        }
    }
    return $List
}

# lmultiply
#
# Multiple a list (in place) by the given value. The modified list is also
# returned.
#
# - This supports the standard lindex syntax for nested lists (ala lindex).
# - NumPy style "broadcasting" is also supported, whereby a list times a value
#   results in each element of the list being multiplied and a list times a
#   list results in element by element multiplication.
#   NB: This only works left to right (ala *=) and modifies the left list.
#
proc lmultiply {list args} {
    upvar 1 $list List
    set value [lindex $args end]
    set index [lrange $args 0 end-1]

    set target [lindex $List {*}$index]
    set vlen [llength $value]
    set tlen [llength $target]

    if {$vlen == 1} {
        if {$tlen == 1} {
            # increment a single element by a single value
            lset List {*}$index [expr {$target*$value}]
        } elseif {$tlen > 1} {
            # increment a each element in a list by a single value
            set i 0
            foreach item $target {
                if {[llength $item] != 1} {
                    error "lmultiply encountered a nested list during list *\
                            value operation"
                }
                lset List {*}$index $i [expr {$item*$value}]
                incr i
            }
        }
    } elseif {$vlen > 1} {
        if {$tlen != $vlen} {
            error "lmultiply cannot perform list * list operation with\
                    different lengths ($tlen != $vlen)"
        }
        set i 0
        foreach item1 $target item2 $value {
            lset List {*}$index $i [expr {$item1*$item2}]
            incr i
        }
    }
    return $List
}

# =============================================================================
# Probability, Statistics, and Random numbers
# =============================================================================
# Generate a Gaussian (normal) random variable with the given mean and standard
# deviation (default to zero and one, respectively). 
#
# This uses a Box-Muller-like transformation by generating two uniform variates
# within the unit circle. The cost of rejecting some random numbers generally
# outweighs the cost of the cos() and sin() needed for the standard Box-Muller
# transformation.
#
variable SecondGaussianWaiting 0
variable SecondGaussian 0.0
proc tcl::mathfunc::normal {{mu 0.0} {sigma 1.0}} {
    if {$::SecondGaussianWaiting} {
        set ::SecondGaussianWaiting 0
        expr {$mu + $sigma*$::SecondGaussian}
    } else {
        # Pick two variates within the unit circle. The lower bound is from
        # NAMD Random.h and purportedly ensures "abs result < 6".
        set r 2.0
        while {$r >=1. || $r < 1.523e-8} {
            set u1 [expr {2*rand() - 1}]
            set u2 [expr {2*rand() - 1}]
            set r [expr {$u1*$u1 + $u2*$u2}]
        }
        # The transformation produces two random variates - store one.
        set fac [expr sqrt(-2*log($r)/$r)]
        set ::SecondGaussianWaiting 1
        set ::SecondGaussian [expr {$fac*$u1}]
        expr {$mu + $sigma*$fac*$u2}
    }
}

# choice
#
# Choose a random element from a list. The index of the element is also
# returned.
#
# The default is to uniformly weight the elements, but (unnormalized)
# non-uniform weights may also be provided.
#
# Arguments:
# ----------
# list : list
#    list to choose an element from
# weights : list of floats (optional)
#   The probability weight of each choice. These need not be normalized.
#
# Returns:
# --------
# element
#    the selected list element
# index : int
#    the index of the selected element
#
proc choice {list {weights {}}} {
    # Uniform weights is really easy.
    if {![llength $weights]} {
        set j [expr {int(rand()*[llength $list])}]
        return [list [lindex $list $j] $j]
    }
    # Non-uniform weighting is more involved.
    if {[llength $list] != [llength $weights]} {
        error "Mismatch in list/weights in choice!"
    }
    set WeightSum [expr {lsum($weights)}]
    set Rand [expr {$WeightSum*rand()}]
    set j 0
    foreach Weight $weights {
        set Rand [expr {$Rand - $Weight}]
        if {$Rand <= 0.0} {
            return [list [lindex $list $j] $j]
        }
        incr j
    }
    # This may be reached if one of the weights is essentially one within
    # rounding precision. In that case, returning the most probable choice is
    # likely a good approximation and really the only stable thing to do.
    set MaxWeight 0.0
    set i 0
    foreach Weight $weights {
        if {$Weight >= $MaxWeight} {
            set MaxWeight $Weight
            set j $i
        }
        incr i
    }
    return [list [lindex $list $j] $j]
}

