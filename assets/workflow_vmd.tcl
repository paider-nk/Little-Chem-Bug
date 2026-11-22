# Input paths travel as UTF-8 hex so Windows Tcl can open non-ASCII names.
proc workflow_text {name} {
    return [encoding convertfrom utf-8 [binary format H* $::env($name)]]
}
proc workflow_display {} {
    set mid [mol new [workflow_text WORKFLOW_VMD_FILE] type $::env(WORKFLOW_VMD_TYPE) waitfor all]
    if {$mid < 0 || [molinfo $mid get numatoms] < 1} {error "No atoms loaded"}
    mol modstyle 0 $mid CPK 1.0 0.3 16 16
    if {$::env(WORKFLOW_VMD_TYPE) eq "cube"} {
        foreach sign {1 -1} color {0 1} {
            mol representation Isosurface [expr {$sign * $::env(WORKFLOW_VMD_ISO)}] 0 0 0 1 1
            mol color ColorID $color
            mol addrep $mid
        }
    }
    display resetview
    return "[molinfo $mid get numatoms] [molinfo $mid get numreps]"
}
set failed [catch {workflow_display} result]
set status [workflow_text WORKFLOW_VMD_STATUS]
set out [open ${status}.tmp w]
fconfigure $out -encoding utf-8
puts $out "$failed $result"
close $out
file rename -force ${status}.tmp $status
if {$failed} {puts stderr $result; exit}
