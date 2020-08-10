package provide orcaviewer 1.0
package require tablelist

proc enabletrace {} {
  global vmd_frame;
  trace variable vmd_frame([molinfo top]) w updateFrame
}

proc disabletrace {} {
  global vmd_frame;
  trace vdelete vmd_frame([molinfo top]) w updateFrame
}


namespace eval OrcaViewer:: {
  namespace export orcaviewer

  variable w                                          ;# handle to main window
  variable orcaFile   ""                  ;# output file
  variable fileLoaded 0
  variable orcaMol ""
  variable orbitals {}
  variable orbitalsInTable 0
  variable lowestOrbitalIndex 0
  variable reps {}
  variable orbSpan 10
}

#
# Create the window
#
proc OrcaViewer::orcaviewer {} {
  variable w
  variable orcaFile
  variable fileLoaded
  variable ::OrcaViewer::orbSpan

  # If already initialized, just turn on
  if { [winfo exists .orcaviewer] } {
    wm deiconify $w
    return
  }

  set w [toplevel ".orcaviewer"]
  wm title $w "OrcaViewer"
  wm resizable $w 0 0

  frame $w.menubar -relief raised -bd 2 ;# frame for menubar
  pack $w.menubar -padx 1 -fill x

  menubutton $w.menubar.help -text Help -underline 0 -menu $w.menubar.help.menu
  menubutton $w.menubar.file -text File -underline 0 -menu $w.menubar.file.menu

  ## help menu
  # menu $w.menubar.help.menu -tearoff no
  # $w.menubar.help.menu add command -label "Help..." -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/textview"
  # XXX - set menubutton width to avoid truncation in OS X
  # $w.menubar.help config -width 5

  menu $w.menubar.file.menu -tearoff no
  $w.menubar.file.menu add command -label "Load" -command  OrcaViewer::loadFile
  $w.menubar.file.menu add command -label "Clean" -command OrcaViewer::cleanUp
  # $w.menubar.file.menu add command -label "Save" -command  TextView::savefile
  # $w.menubar.file.menu add command -label "Save As" -command  TextView::saveasfile
  $w.menubar.file config -width 5
  pack $w.menubar.file -side left
  pack $w.menubar.help -side right


  ##
  ## main window area
  ##
  wm resizable $w 1 1
  wm geometry $w 500x200
  # frame $w.main
  scrollbar $w.scr -orient vertical -command [list $w.main  yview]
  tablelist::tablelist $w.main -columns {0 "Orb. Number" 0 "Descr." 0 "Energy"} -stretch all -background white \
    -selectmode extended \
    -exportselection true \
    -state normal -selectmode extended \
    -labelrelief groove \
    -labelbd 1 \
    -yscrollcommand [list $w.scr set] \
    -editstartcommand OrcaViewer::startOrbitalClick -editendcommand OrcaViewer::updateOrbitals

    bind $w.main <<TablelistSelect>> {
        set col [tablelist::getTablelistColumn %W]
        set row [%W curselection]
        puts $row
        OrcaViewer::updateOrbitals $row $col
    }
  pack $w.main -fill both -expand 1 -side top
  #$w.main insert end [list "first row" "another value"]
  #$w.main insert end [list "another row" "bla bla"]

  # label $w.txt.label -width 80 -relief sunken -bg White -textvariable TextView::textfile
  # text $w.txt.text -bg White -bd 2 -yscrollcommand "$::TextView::w.txt.vscr set"
  # scrollbar $w.txt.vscr -command "$::TextView::w.txt.text yview"
  # pack $w.txt.label
  # pack $w.txt.text $w.txt.vscr -side left -fill y
  pack $w.menubar $w.main -side top -expand 1
  entry $w.span -textvariable ::OrcaViewer::orbSpan -width 10
  bind $w.span <Return> { OrcaViewer::fillOrblist }
  label $w.spanLabel -text "Orbitals around HOMO: "
  pack $w.main $w.span -expand 1 -side bottom
  pack $w.spanLabel $w.span -expand 1 -side right
}

proc OrcaViewer::updateOrbitals {row col} \
{
  variable orbitals
  variable orcaMol
  variable lowestOrbitalIndex
  variable reps
  set orb [expr $row+$lowestOrbitalIndex+1]

  if {[llength $orbitals] && ![llength $reps]} {
    #puts "update orbital $row"
    puts [lindex $orbitals [expr $row]]
    mol color ColorID 0
    mol representation Orbital 0.050000 $orb 0 0 0.125 1 0 0 0 1
    mol selection all
    mol material Glossy
    mol addrep $orcaMol
    mol color ColorID 3
    mol representation Orbital -0.050000 $orb 0 0 0.125 1 0 0 0 1
    mol addrep $orcaMol
    set nr [molinfo $orcaMol get numreps]
    set reps [list [mol repname $orcaMol [expr $nr - 1]] [mol repname $orcaMol [expr $nr - 2]]]
  } elseif {[llength $orbitals] && [llength $reps]} {
    set idx1 [mol repindex $orcaMol [lindex $reps 0]]
    set idx2 [mol repindex $orcaMol [lindex $reps 1]]
    mol modstyle $idx1 $orcaMol Orbital 0.050000 $orb 0 0 0.125 1 0 0 0 1
    mol modstyle $idx2 $orcaMol Orbital -0.050000 $orb 0 0 0.125 1 0 0 0 1
  }

}

proc OrcaViewer::getOrbitalDescr {index hindex} \
{
  set diff [expr $hindex-$index]
  if {$diff == 0} {
      return "HOMO"
  } elseif {$diff > 0} {
      return "HOMO-$diff"
  } elseif {$diff < 0} {
    if {$diff == -1} {
      return "LUMO"
    } else {
      return "LUMO+[expr -$diff-1]"
    }
  }

}

proc OrcaViewer::assignPCtoUserField { method } {
  variable orcaMol
  set nf [molinfo $orcaMol get numframes]
  set sel [atomselect $orcaMol all]
  for {set i 0} {$i < $nf} {incr i} {
    $sel frame $i
    animate goto $i
    set qmcharges [molinfo $orcaMol get qmcharges]
    foreach var $qmcharges {
      if {[lindex $var 0 0] == $method} {
        $sel set user [lindex $var 0 1]
        break
      }
    }
  }
  $sel delete
}

proc updateFrame { name element op } {
  OrcaViewer::fillOrblist
}

proc OrcaViewer::fillOrblist { } \
{
  variable w
  variable orcaMol
  variable orbitals
  variable orbitalsInTable
  variable lowestOrbitalIndex
  variable ::OrcaViewer::orbSpan
  global vmd_frame

  #puts $orbitalsInTable
  $w.main delete 0 $orbitalsInTable
  set orbitals {}

  set energies [lindex [molinfo $orcaMol get orbenergies] 0 0]
  set norbs [llength $energies]
  #puts "norbs: $norbs"
  set homo [molinfo $orcaMol get homo]

  set orbitalsInTable 0
  set lowestOrbitalIndex [expr $homo - $::OrcaViewer::orbSpan]
  if {$lowestOrbitalIndex < 0} {
    set lowestOrbitalIndex 0
  }
  set limit [expr $homo + $::OrcaViewer::orbSpan]
  if {$limit > $norbs} {
    set limit [expr $norbs - 1]
  }
  for {set i $lowestOrbitalIndex} {$i <= $limit} {incr i} {
    lappend orbitals [list $i [OrcaViewer::getOrbitalDescr $i $homo] [lindex $energies $i]]
    incr orbitalsInTable
  }

  foreach orbital $orbitals {
    $w.main insert end $orbital
  }

}


proc OrcaViewer::loadFile { } {
  variable w
  variable orcaFile
  variable fileLoaded
  variable orcaMol

  set file_types {
    {"All Files" * }
  }

  set orcaFile [tk_getOpenFile -filetypes $file_types \
                -initialdir pwd \
                -defaultextension .txt]

  set rc [ catch { set fd [open $orcaFile "r"] } ]
  if { $rc == 1} {
    return
  }

  set fileLoaded 1
  puts $orcaFile
  puts $fileLoaded

  set orcaMol [mol load orca $orcaFile]
  puts $orcaMol
  #OrcaViewer::fillOrblist

  enabletrace
  OrcaViewer::assignPCtoUserField "Mulliken"
  color scale method BWR
  animate goto 0

  close $fd
}

proc OrcaViewer::cleanUp { } \
{
  variable w
  variable orcaFile
  variable fileLoaded
  variable orcaMol
  variable orbitalsInTable
  variable reps
  variable orbitals
  variable lowestOrbitalIndex

  mol delete $orcaMol
  $w.main delete 0 $orbitalsInTable
  set fileLoaded 0
  set orcaMol ""
  set orcaFile ""
  set reps {}
  set orbitals {}
  set orbitalsInTable 0
  set lowestOrbitalIndex 0
}

proc orcaviewer_tk {} {
  OrcaViewer::orcaviewer
  return $OrcaViewer::w
}

vmd_install_extension orcaviewer orcaviewer_tk "Visualization/Orca"
orcaviewer_tk
enabletrace
