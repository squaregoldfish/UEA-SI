#!/usr/bin/tclsh

set indir [lindex $argv 0]
set outdir [lindex $argv 1]

set jobList ""
set jobText ""

for {set lon 1} {$lon <= 144} {incr lon} {
    for {set lat 1} {$lat <= 72} {incr lat} {
        append jobList "${lon}_${lat} "

        append jobText "${lon}_${lat} :\n"
        append jobText "\tR --no-save --slave lon=\"${lon}\" lat=\"${lat}\" indir=\"${indir}\" outdir=\"${outdir}\" < do_interpolation.R\n\n"
    }
}

set outChan [open "Makefile" w]
puts $outChan "all : $jobList\n"
puts $outChan $jobText
close $outChan
