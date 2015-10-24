#!/bin/sh

#  runSimAndAnimate.sh
#  
#
#  Created by Håkon Austlid Taskén on 24.10.15.
#
#
# Takes the name of the textfile spesifing the values used to run the simulation and run schrodingerFD and
# plotSchrodinger.py with this input. Finally it playes the movie produced with quicktime player.

filename=$1
./schrodingerFD << EOF
$filename
EOF
python plotSchrodinger.py << EOF
$filename
EOF
open -a quicktime\ player simulationMovies/"${1:0:${#filename}-4}.mp4"