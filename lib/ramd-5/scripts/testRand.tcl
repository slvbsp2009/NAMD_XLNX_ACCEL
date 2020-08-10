proc genRandAccel {  } {
 set pi 3.141592653589793238462643383
 set theta [expr "2*$pi*rand()"]
 set psi [expr "$pi*rand()"]
 set psi [expr "acos(1-2*rand())"]
 set sinpsi [expr "sin($psi)"]
 set rx [expr "cos($theta)*$sinpsi"]
 set ry [expr "sin($theta)*$sinpsi"]
 set rz [expr "cos($psi)"]
 set r "$rx $ry $rz"
 puts $r
}
expr "srand(100)"
for {set i 0} {$i<50} {incr i} {
  genRandAccel
}
