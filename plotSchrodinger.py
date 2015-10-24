# This script plot simulations run by 'Schrodinger.h' and 'Schrodinger.cpp' and produces a movie.
# One only needs to specify the filename of the textfile with the values used to run the c++ program created from 'Schrodinger.h' and 'Schrodinger.cpp',
# as well as placing these files in the same directory as this file, 'plotSchroinger.py', or in a subdirectory of this directory.
# One do need ffmpeg to make a movie of the plot
# If one wishes to plot the animation in stead of saving it as a movie because one does not have ffmpeg,
# comment out anim.save... and uncomment plt.show() close to the bottom of the file
# If one does not use a mac one might want to comment out matplotlib.use('TKAgg')

print "plotSchrodinger.py is running"
print "Type the name of the situation file:"
filename = raw_input()
print "You typed: ",
print filename
filename = filename[:-4] # remove .txt

animationTime = 10; #sec sets the length of the movie produced in seconds

# one should not be needing to do changes to the rest of the script, exept if one does not have ffmpeg and would like to show the plot instead of producing a movie
# or possibly to comment out matplotlib.use('TKAgg') if one does not use a mac

import numpy as np
import matplotlib
matplotlib.use('TKAgg') # this is done to allow blit = True in FuncAnimation on mac
# import proper graphics back-end for Mac OS X
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
import os
import time
import math
import customColormaps as ccmaps
import sys

# find location of file, 'name'
def find(name):
    for root, dirs, files in os.walk(os.path.dirname(os.path.realpath(__file__))):
        if name in files:
            return os.path.join(root, name)

# find the right values in situation file
def getValue(file):
    value = ""
    while True:
        c = file.read(1)
        if c == ':':
            c = file.read(1)
            while c != '#':
                if c == '"':
                    value = ""
                    c = file.read(1)
                    while c != '"':
                        value += c
                        c = file.read(1)
                    return value
                value += c
                c = file.read(1)
            return value
        if not c:
            break
    return value

# time the animation
startTime = time.clock()

# load variables used in simulation
try:
    situationFile = open(find(filename + ".txt"), mode = 'r')
except IOError:
     sys.exit("Could not open file.")
numOfDim = int(getValue(situationFile))
potential = getValue(situationFile)
probDistrb = getValue(situationFile)
m = float(getValue(situationFile))
hbar = float(getValue(situationFile))
Ni = int(getValue(situationFile))
numOfFrames = int(getValue(situationFile))
Nx1 = int(getValue(situationFile))
Nx2 = int(getValue(situationFile))
Nx3 = int(getValue(situationFile))
Lx1 = float(getValue(situationFile))
Lx2 = float(getValue(situationFile))
Lx3 = float(getValue(situationFile))
plotResX1 = int(getValue(situationFile))
plotResX2 = int(getValue(situationFile))
plotResX3 = int(getValue(situationFile))
situationFile.close()

simulationValues = open(find(filename + "_simulationValues.txt"), mode = 'r');
simulationValues.readline()
startEnergy = float(simulationValues.readline())
simulationValues.readline()
finalEnergy = float(simulationValues.readline())
simulationValues.readline()
finalProb = float(simulationValues.readline())
simulationValues.readline()
Vmax = float(simulationValues.readline())
simulationValues.readline()
simulationValues.readline()
simulationTime = float(simulationValues.readline())
simulationValues.readline()
plotSpacingI = int(simulationValues.readline())
simulationValues.readline()
plotSpacingX1 = int(simulationValues.readline())
simulationValues.readline()
plotSpacingX2 = int(simulationValues.readline())
simulationValues.readline()
plotSpacingX3 = int(simulationValues.readline())
simulationValues.readline()
animationTimeExists = simulationValues.readline()

# makes sure numOfFrames and plottedResolution is correct
numOfFrames = Ni / plotSpacingI
plotResX1 = Nx1 / plotSpacingX1
plotResX2 = Nx2 / plotSpacingX2
plotResX3 = Nx3 / plotSpacingX3

print "The final energy is ",
print finalEnergy / startEnergy,
print " times the start energy."
print "The final probability of finding the particle is: ",
print finalProb
print "The time used to run simulation: ",
print simulationTime

# load data we would like to plot from various files
dt = np.dtype("f8")
plotProbabilityFile = np.fromfile(find(filename + "_plot_probability"), dtype=dt)
plotPsiRFile = None
plotPsiIFile = None
if numOfDim == 1:
    plotPsiRFile = np.fromfile(find(filename + "_plot_psi_r"), dtype=dt)
    plotPsiIFile = np.fromfile(find(filename + "_plot_psi_i"), dtype=dt)
potentialFile = np.fromfile(find(filename + "_potential"), dtype=dt)

maxProb = max(plotProbabilityFile)
minProb = min(plotProbabilityFile)

# define some scaling constants
scaleConstPsi = None
if numOfDim == 1:
    scaleConstPsi = 0.8 * maxProb / (max([max(plotPsiRFile),max(plotPsiIFile)]))

scaleConstAx4 = max(plotProbabilityFile[0:plotResX1*plotResX2])/maxProb
for i in range(1,numOfFrames):
    if max(plotProbabilityFile[plotResX1*plotResX2*i:plotResX1*plotResX2*(i+1)])/maxProb < scaleConstAx4:
        scaleConstAx4 = max(plotProbabilityFile[plotResX1*plotResX2*i:plotResX1*plotResX2*(i+1)])/maxProb;
scaleConstAx2 = scaleConstAx4 + (1 - scaleConstAx4) * 2 / 3.0
scaleConstAx3 = scaleConstAx4 + (1 - scaleConstAx4) / 3.0

scaleConstEnergy = 0.5 * maxProb / startEnergy

# initialization function for 1D: plot the background of each frame
def init1D():
    plt.plot(x1, scaleConstEnergy * potentialFile[0:Nx1:plotSpacingX1], ':k', zorder=0)
    pylab.fill(x1, scaleConstEnergy * potentialFile[0:Nx1:plotSpacingX1], facecolor='y', alpha=0.2, zorder=0, label = 'Potential')
    energy = np.linspace(scaleConstEnergy * startEnergy, scaleConstEnergy * startEnergy, plotResX1)
    energyPlot, = ax.plot(x1, energy, 'g', zorder = 1, label = "Energy")
    plt.legend(loc = 'lower right', prop={'size':10})
    probPlot.set_data([], [])
    psiRPlot.set_data([], [])
    psiIPlot.set_data([], [])
    return probPlot, psiRPlot, psiIPlot, energyPlot,

# animation function for 1D.  This is called sequentially
def animate1D(i):
    probPlot.set_data(x1, plotProbabilityFile[plotResX1*i:plotResX1*(i+1)])
    psiRPlot.set_data(x1, scaleConstPsi * plotPsiRFile[plotResX1*i:plotResX1*(i+1)])
    psiIPlot.set_data(x1, scaleConstPsi * plotPsiIFile[plotResX1*i:plotResX1*(i+1)])
    return probPlot, psiRPlot, psiIPlot,

# initialization function for 2D: plot the background of each frame
def init2D():
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,plotResX1), np.linspace(0,Lx2,plotResX2))
    z = potentialFile[0:Nx1*Nx2:plotSpacingX1].reshape(Nx2,plotResX1)[0:Nx2:plotSpacingX2,:]
    colorLimit = max(abs(potentialFile))
    if colorLimit == 0:
        colorLimit = 1
    plt.sca(ax1)
    potentialPlot1 = plt.contourf(x1, x2, z, cmap=ccmaps.cmap('white_black'), zorder = 2, vmin = - colorLimit, vmax = colorLimit)
    plt.sca(ax2)
    potentialPlot3 = plt.contourf(x1, x2, z, cmap=ccmaps.cmap('white_black'), zorder = 2, vmin = - colorLimit, vmax = colorLimit)
    plt.sca(ax3)
    potentialPlot2 = plt.contourf(x1, x2, z, cmap=ccmaps.cmap('white_black'), zorder = 2, vmin = - colorLimit, vmax = colorLimit)
    plt.sca(ax4)
    potentialPlot4 = plt.contourf(x1, x2, z, cmap=ccmaps.cmap('white_black'), zorder = 2, vmin = - colorLimit, vmax = colorLimit)
    return potentialPlot1, potentialPlot2, potentialPlot3, potentialPlot4,

# animation function for 2D.  This is called sequentially
def animate2D(i):
    z = plotProbabilityFile[plotResX1*plotResX2*i:plotResX1*plotResX2*(i+1)].reshape(plotResX2,plotResX1)
    plt.sca(ax1)
    probPlot1 = plt.contourf(x1,x2,z, cmap=ccmaps.cmap('vaagen_colorscale'), zorder = 1, vmin = 0, vmax = maxProb)
    plt.sca(ax2)
    probPlot2 = plt.contourf(x1,x2,z, cmap=ccmaps.cmap('vaagen_colorscale'), zorder = 1, vmin = 0, vmax = maxProb*scaleConstAx2)
    plt.sca(ax3)
    probPlot3 = plt.contourf(x1,x2,z, cmap=ccmaps.cmap('vaagen_colorscale'), zorder = 1, vmin = 0, vmax = maxProb*scaleConstAx3)
    plt.sca(ax4)
    probPlot3 = plt.contourf(x1,x2,z, cmap=ccmaps.cmap('vaagen_colorscale'), zorder = 1, vmin = 0, vmax = maxProb*scaleConstAx4)
    return probPlot1, probPlot2, probPlot3, probPlot4

# initialization function for 3D: plot the background of each frame
def init3D():
    return

# animation function for 3D.  This is called sequentially
def animate3D(i):
    probMat = plotProbabilityFile[plotResX1*plotResX2*plotResX3*i:plotResX1*plotResX2*plotResX3*(i+1)].reshape(plotResX3,plotResX2,plotResX1)
    z11 = np.sum(probMat, axis=0)
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,plotResX1), np.linspace(0,Lx2,plotResX2))
    plt.sca(ax1)
    probPlot11 = plt.contourf(x1, x2, z11, cmap=plt.cm.gnuplot2_r, alpha = 1, zorder = 1, label = 'xy-axis')
    z12 = np.sum(probMat, axis=1)
    x1, x3 = np.meshgrid(np.linspace(0,Lx1,plotResX1), np.linspace(0,Lx3,plotResX3))
    plt.sca(ax2)
    probPlot12 = plt.contourf(x1, x3, z12, cmap=plt.cm.gnuplot2_r, alpha = 1, zorder = 1, label = 'xz-axis')
    z21 = np.sum(probMat, axis=2)
    x2, x3 = np.meshgrid(np.linspace(0,Lx2,plotResX2), np.linspace(0,Lx3,plotResX3),)
    plt.sca(ax3)
    probPlot21 = plt.contourf(x2, x3, z21, cmap=plt.cm.autumn_r, alpha = 1, zorder = 1, label = 'yz-axis')
    # http://stackoverflow.com/questions/9164950/python-matplotlib-plot3d-with-a-color-for-4d
    return

anim = None
if numOfDim == 1:
    fig = plt.figure()
    ax = plt.axes(xlim=(0, Lx1), ylim=(1.1 * scaleConstPsi * min([min(plotPsiRFile), min(plotPsiIFile)]), 1.1 * maxProb))
    probPlot, = ax.plot([], [], 'k', lw = 1, zorder = 3, label = 'Probability')
    psiRPlot, = ax.plot([], [], 'b', lw = 1, zorder = 2, label = 'Re($\Psi$)') # only used for 1D
    psiIPlot, = ax.plot([], [], 'r', lw = 1, zorder = 2, label = 'Im($\Psi$)') # only used for 1D
    x1 = np.linspace(0,Lx1,plotResX1)
    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate1D, init_func=init1D, frames = numOfFrames, interval=1000*animationTime/numOfFrames, blit=True, repeat = False)
elif numOfDim == 2:
    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    probPlot1 = None
    probPlot2 = None
    probPlot3 = None
    probPlot4 = None
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,plotResX1), np.linspace(0,Lx2,plotResX2))
    anim = animation.FuncAnimation(fig, animate2D, init_func=init2D, frames = numOfFrames, interval=1000*animationTime/numOfFrames, blit=True, repeat = False)
elif numOfDim == 3:
    anim = animation.FuncAnimation(fig, animate3D, init_func=init3D, frames = numOfFrames, interval=1000*animationTime/numOfFrames, blit=True, repeat = False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('simulationMovies/' + filename + '.mp4', fps=numOfFrames/animationTime, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])

#plt.show()

if  len(animationTimeExists) == 0:
    simulationValues = open(find(filename + "_simulationValues.txt"), mode = 'a');
    simulationValues.write(str(time.clock() - startTime))
    simulationValues.close()

print "Seconds used to animate the simulation: ",
print time.clock() - startTime

