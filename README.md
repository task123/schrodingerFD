# schrodingerFD
Simulates the schrodinger equation using finite difference method in 1, 2 and 3 dimensions 
and animates the simulation. (The animation of 3D does not work yet.)

The program can simulate time evolution of any initial condition and any potential as long as it only depent on 
the position.

To run a simulation spesify the values you would like for the simulation in a textfile with a spesial format. 
'freeElectron1D.txt' gives an simple example of this format, and describe the format of the file.
This textfile should be in the directory 'situations', but can also be in the same directory as 'schrodingerFD'.
The textfile 'initialStatesAndPotentials.txt' describes what initial states and potential are available and what 
extra parameters these needs.

Then run the 'schrodingerFD' program and type in the name of this textfile as input.
If there exists a directory with the same name as this file then all files produced by 'schrodingerFD' is placed 
in this directory (recommended).

This simulation can then be animated by running 'plotSchrodinger.py' which produces a movie placed in the 
'simulationMovies' directory. One need ffmpeg to make this movie. If one does not have ffmpeg, it is still 
possible to plot the simulation making a minor change as descibed in 'plotSchrodinger.py'.

For mac (and probably linux) users one can simply run runSimAndAnimate.sh and type in the name of the textfile 
with the values to simulate, to both run the simulation, animate it and open the movie with quicktime player.
The runSimAndAnimate.sh must first be given permission by running '$ chmod +x runSimAndAnimate.sh'.
Example call: ./runSimAndAnimate.sh free1D.txt
Likewise contSimAndAnimation.sh run a continuation of a previous simulation, animate it open the movie with
quicktime player. contSimAndAnimation takes three input parameters; the name of the file with the values we 
simulate, number of interations and if output should be appendend to previous files as input.
Example call: ./contSimAndAnimate.sh free1D.txt 1000 1

One can make new initial states and potentials in the 'Schrodinger.cpp' file found in 'sourceFilesSchrodingerFD'
under setV() and makeInitState(). Remember to add them to the 'initialStatesAndPotentials.txt'

Authers: Håkon Taskén and Paul Thrane
