#!/bin/sh

#  runSimAndAnimate.sh
#  
#
#  Created by Håkon Austlid Taskén on 24.10.15.
#
#
# This shell script runs schrodingerFD and plotSchrodinger and plays the movie produced in quicktime player.
# Require the name of the textfile spesifing the values used to run the simulation as input.
# Example: $ ./runSimAndAnimate.sh free1D.txt

filename=$1
./schrodingerFD << EOF
$filename
EOF
python plotSchrodinger.py << EOF
$filename
EOF
open -a quicktime\ player simulationMovies/"${1:0:${#filename}-4}.mp4"