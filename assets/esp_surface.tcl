# Wavefunction-derived density (volume 0) colored by total ESP (volume 1).
# Paths are relative to this scene, including when reopened manually.
proc esp_scene {folder} {
    set mid [mol new [file join $folder density.cub] type cube waitfor all]
    mol addfile [file join $folder totesp.cub] type cube waitfor all molid $mid
    if {[molinfo $mid get numvolumedata] != 2} {error "Expected two volumetric datasets"}
    mol modstyle 0 $mid CPK 0.9 0.28 28 24
    mol modcolor 0 $mid Element
    mol representation Isosurface @ISO@ 0 0 0 1 1
    mol color Volume 1
    mol material Transparent
    material change opacity Transparent 0.55
    material change specular Transparent 0.14
    mol addrep $mid
    mol scaleminmax $mid 1 -@LIMIT@ @LIMIT@
    color scale method RGB
    color scale colors RGB {0.95 0.04 0.02} {0.98 0.90 0.22} {0.00 0.83 0.96}
    color scale midpoint 0.5
    color change rgb 32 0.50 0.50 0.79
    color Display Background 32
    display backgroundgradient off
    display projection Orthographic
    display depthcue off
    axes location Off
    display resetview
    molinfo $mid set rotate_matrix [list {@ROTATION@}]
    scale by 0.85
    display update
}
set folder [file dirname [file normalize [info script]]]
set failed [catch {
    esp_scene $folder
    if {[info exists ::env(ESP_RENDER)]} {
        render imagesize TachyonInternal {1600 1000}
        render TachyonInternal [file join $folder surface.tga]
        set file [open [file join $folder palette.txt] w]
        for {set i [colorinfo num]} {$i < [colorinfo max]} {incr i} {puts $file [colorinfo rgb $i]}
        close $file
    }
} detail]
if {[info exists ::env(ESP_STATUS)]} {
    set status [file join $folder $::env(ESP_STATUS)]
    set f [open ${status}.tmp w]; puts $f [list $failed $detail]; close $f
    file rename -force ${status}.tmp $status
}
if {$failed} {puts stderr $detail}
if {$failed || [info exists ::env(ESP_RENDER)]} {exit}

