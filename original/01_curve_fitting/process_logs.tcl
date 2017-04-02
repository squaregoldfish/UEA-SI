#!/usr/bin/tclsh

set dir [lindex $argv 0]

set successChan [open "successCells.csv" w]

set cellCount 0
set cellTotal 10368.0

for {set lon 1} {$lon <= 144} {incr lon} {
    for {set lat 1} {$lat <= 72} {incr lat} {
        incr cellCount
        puts -nonewline "\r [format "%.2f" [expr $cellCount / $cellTotal * 100]]%    "
        flush stdout

        set logfile "${dir}/${lon}_${lat}.log"

        if {[file exists $logfile]} {

            set chan [open $logfile r]
            set data [read $chan]
            close $chan

            if {[string first "SUCCESS" $data] > -1} {
                puts $successChan "[expr ($lon * 2.5) - 1.25],[expr ($lat * 2.5) - 91.25],1"
            }
        }
    }
}    

close $successChan
puts ""

