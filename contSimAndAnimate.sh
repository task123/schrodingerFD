#!/bin/sh

#  contSimAndAnimate.sh
#  
#
#  Created by Håkon Austlid Taskén on 24.10.15.
#
#
# This shell script runs continueSchrodingerFD and plotSchrodinger and plays the movie produced
# in quicktime player.
# Require the name of the textfile spesifing the values used to run the simulation, number of interations and
# '1' or '0' representing true or false on appending the output to previous files.
# Example: $ ./runSimAndAnimate.sh free1D.txt 1000 1

filename=$1
Ni=$2
appendOldFile=$3
./continueSchrodingerFD << EOF
$filename
$Ni
$appendOldFile
EOF
python plotSchrodinger.py << EOF
$filename
EOF
open -a quicktime\ player simulationMovies/"${1:0:${#filename}-4}.mp4"