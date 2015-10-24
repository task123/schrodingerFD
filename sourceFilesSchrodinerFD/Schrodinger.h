//
//  Schrodinger.h
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//


#ifndef __schrodingerFD__Schrodinger__
#define __schrodingerFD__Schrodinger__

#include <stdio.h>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// Uses finite difference method to simulate the schrodinger equation in 1, 2 or 3 dimensions.
// The simulation can be run with a variety of initial conditions and parameters spesified in a textfile in a spesial format given as input or a function parameter.
// If you wish to change the initial condition or other parameter do so in a textfile with this spesial format.
// If you wish to add a new potential or initial state, do so in setV and makeInitCondition respectivly. Do also add it to the textfile initialStatesAndPotentials.txt together with what parameters it require
class Schrodinger{
public:
    // run the simulation and asks for the name of a text file where the 'situation' simulated are spesified in a spesial form
    void run();
    
    // run the simulation with the 'situation' spesified by values in the text file filename
    void run(string filename);
    
    //contunue previous simulation and asks for the name of the file with the values being simulated, number of intrations and if it should append the output to the previous file
    // not appending the output to the previous file will result in a new scaling of plot in plotSchrodinger.py
    void continueSimulation();
    
    // continue previous simulation
    // appendOldFile = false will result in a new scaling of plot in plotSchrodinger.py
    void continueSimulation(string filename, int numOfIterations, bool appendOldFile);
    
    Schrodinger();
    ~Schrodinger();
private:
// MEMBER FUNCTIONS
    void storeFilename(string filename);
    void loadAndCalculateVariables();
        string getValue(ifstream &situationFile);

    // setting initial state
    void makeInitState();
        void makeInitState1D();
        void makeInitState2D();
        void makeInitState3D();
            void normalizePsi();
            double findProbability();
    
    // spesifies the potential
    void setV();
        void setV1D();
        void setV2D();
        void setV3D();
            // sets V to zero everywhere
            void setVtoZero();
    
    double findEnergy();
        double findEnergy1D();
        double findEnergy2D();
        double findEnergy3D();
    
    // check if there exists an directory named 'filename' and changes the filename to filename/filename if so. Such that all files produced to the simulation is placed here.
    void checkForDirectory(bool newSimulation);
    
    // calculating the time evolution and storing it in a tekst file
    void finiteDifference(bool newSimulation);
        void finiteDifference1D(char* fileOpenType);
        void finiteDifference2D(char* fileOpenType);
        void finiteDifference3D(char* fileOpenType);
    
    void storeFinalState(long time);
    void loadFinalState();
    
// VARIABLES FOR SIMULATION
    double* psi_r1; //real part of wave function in timesteps 3t
    double* psi_i1; //imagenary part of wave funciton in timesteps 3t
    double* psi_r2; //real part of wave function in timesteps 1 + 3t
    double* psi_i2; //imagenary part of wave funciton in timesteps 1 + 3t
    double* psi_r3; //real part of wave function in timesteps 2 + 3t
    double* psi_i3; //imagenary part of wave funciton in timesteps 2 + 3t
    double* V;
    
    string filename;
    int numOfDim;
    string potential;
    string probDistrb; // name of initial probability distribution
    double m;
    double hbar;
    int Ni; //number of iterations (timesteps = 3*Nt)
    int numOfFrames; // Nt / plotSpacingT
    int Nx1;
    int Nx2;
    int Nx3;
    double Lx1;
    double Lx2;
    double Lx3;

    double dx1;
    double dx2;
    double dx3;
    
    int plotSpacingI; //iterations between plotted
    int plotSpacingX1; //steps between x1 plotted
    int plotSpacingX2; //steps between x2 plotted
    int plotSpacingX3; //steps between x3 plotted
    
    double dt;
    
    double Vmax;
    double startEnergy;
    double finalEnergy;
    double finalProb;

    ifstream* situationFile;
};




#endif /* defined(__schrodingerFD2D__Schrodinger2D__) */
