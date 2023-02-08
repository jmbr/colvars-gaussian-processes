#############################################################################
# Evaluate Gaussian process regressor.
#############################################################################

namespace eval gaussian_process {
    # set epsilon ...
    # set coefficients(0) { ... }
    # set coefficients(1) { ... }
}

proc calc_gaussian_process { args } {
    upvar gaussian_process::epsilon epsilon \
	  gaussian_process::coefficients coefficients

    set n [llength $args]
    set m [array size coefficients]

    set values [lrepeat $m 0.0]

    for {set i 0} {$i < $n} {incr i} {
	set rmsd [lindex $args $i]
	set exponential [expr exp(-1 * $rmsd**2 / (2.0 * $epsilon**2))]
	for {set j 0} {$j < $m} {incr j} {
	    set coefficient [lindex $coefficients($j) $i]
	    set value [expr $coefficient * $exponential]
	    lset values $j [expr [lindex $values $j] + $value]
	}
    }

    return $values
}

proc calc_gaussian_process_gradient { args } {
    upvar gaussian_process::epsilon epsilon \
          gaussian_process::coefficients coefficients

    set n [llength $args]
    set m [array size coefficients]

    set jacobian [lrepeat $n [lrepeat $m 0.0]]

    for {set i 0} {$i < $n} {incr i} {
	set rmsd [lindex $args $i]
	set exponential [expr exp(-1 * $rmsd**2 / (2.0 * $epsilon**2))]
	for {set j 0} {$j < $m} {incr j} {
	    set coefficient [lindex $coefficients($j) $i]
	    set value [expr -1 * $coefficient * $rmsd / $epsilon**2 * $exponential]
	    lset jacobian $i $j $value
	}
    }

    return $jacobian
}
