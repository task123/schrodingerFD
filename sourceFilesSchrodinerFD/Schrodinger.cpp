//
//  Schrodinger.cpp
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include "Schrodinger.h"
#include "float.h"
#include <ctime>

using namespace std;

//PUBLIC MEMBER FUNCTIONS
void Schrodinger::run(){
    cout << "Type in the name of the situation file:" << endl;
    string typedFilename;
    cin >> typedFilename;
    cout << "You typed in: " << typedFilename << endl;
    run(typedFilename);
}

void Schrodinger::run(string filename){
    clock_t t = clock() / CLOCKS_PER_SEC;
    filename = filename.substr(0, filename.size()-4);
    this->filename = filename;
    
    ifstream sitFile(filename + ".txt");
    while (!sitFile.is_open()) {
        cout << "Could not open the file: " << filename << endl;
        cout << "Type in the name of the situation file:" << endl;
        cin >> filename;
        sitFile.open(filename);
        filename = filename.substr(0, filename.size()-4);
    }
    situationFile = &sitFile;
    loadAndCalculateVariables();
    
    V = new double [Nx1 * Nx2 * Nx3];
    psi_r1 = new double [Nx1 * Nx2 * Nx3];
    psi_i1 = new double [Nx1 * Nx2 * Nx3];
    psi_r2 = new double [Nx1 * Nx2 * Nx3];
    psi_i2 = new double [Nx1 * Nx2 * Nx3];
    psi_r3 = new double [Nx1 * Nx2 * Nx3];
    psi_i3 = new double [Nx1 * Nx2 * Nx3];
    
    makeInitState();
    setV();
    startEnergy = findEnergy();

    checkForDirectory(true);
    
    finiteDifference(true);
    
    finalEnergy = findEnergy();
    finalProb = findProbability();
    cout << "The start energy is: " << startEnergy << " and the highest potential is: " << Vmax << endl;
    cout << "The final energy is " << finalEnergy/startEnergy << " times the start energy." << endl;
    cout << "The final probability of finding the particle is: " << finalProb << endl;
    long time = static_cast<long>(clock() - t) / CLOCKS_PER_SEC;
    cout << "The simulation used " << time / 60 << " minuttes and " << time % 60 << " seconds." << endl;
    
    storeFinalState(time);
}

void Schrodinger::continueSimulation(string filename, int numOfIterations, bool appendOldFile){
    clock_t t = clock() / CLOCKS_PER_SEC;
    this->filename = filename;
    loadAndCalculateVariables();
    loadFinalState();
    Ni = numOfIterations;
    
    V = new double [Nx1 * Nx2 * Nx3];
    psi_r1 = new double [Nx1 * Nx2 * Nx3];
    psi_i1 = new double [Nx1 * Nx2 * Nx3];
    psi_r2 = new double [Nx1 * Nx2 * Nx3];
    psi_i2 = new double [Nx1 * Nx2 * Nx3];
    psi_r3 = new double [Nx1 * Nx2 * Nx3];
    psi_i3 = new double [Nx1 * Nx2 * Nx3];

    setV();
    
    checkForDirectory(!appendOldFile);

    finiteDifference(!appendOldFile);
    cout << "The start energy is: " << startEnergy << " and the highest potential is: " << Vmax << endl;
    cout << "The final energy is " << finalEnergy/startEnergy << " times the start energy." << endl;
    cout << "The final probability of finding the particle is: " << finalProb << endl;
    long time = static_cast<long>(clock() - t) / CLOCKS_PER_SEC;
    cout << "The simulation used " << time / 60 << " minuttes and " << time % 60 << " seconds." << endl;
    
    storeFinalState(time);
}

Schrodinger::Schrodinger(){
    V = nullptr;
    psi_r1 = nullptr;
    psi_i1 = nullptr;
    psi_r2 = nullptr;
    psi_i2 = nullptr;
    startEnergy = 0.0;
    finalEnergy = 0.0;
}

Schrodinger::~Schrodinger(){
    delete [] V;
    delete [] psi_r1;
    delete [] psi_i1;
    delete [] psi_r2;
    delete [] psi_i2;
}

//PRIVATE MEMBER FUNCTIONS
void Schrodinger::loadAndCalculateVariables(){
    // load variables
    numOfDim = stoi(getValue(*situationFile));
    potential = getValue(*situationFile);
    probDistrb = getValue(*situationFile);
    m = stod(getValue(*situationFile));
    hbar = stod(getValue(*situationFile));
    Ni = stoi(getValue(*situationFile));
    numOfFrames = stoi(getValue(*situationFile));
    Nx1 = stoi(getValue(*situationFile));
    Nx2 = stoi(getValue(*situationFile));
    Nx3 = stoi(getValue(*situationFile));
    Lx1 = stod(getValue(*situationFile));
    Lx2 = stod(getValue(*situationFile));
    Lx3 = stod(getValue(*situationFile));
    int plottedResolutionX1 = stoi(getValue(*situationFile));
    int plottedResolutionX2 = stoi(getValue(*situationFile));
    int plottedResolutionX3 = stoi(getValue(*situationFile));
    // do not close file on purpose to allow more variables to be loaded to spesify the potential and probDistrib further
    
    // calculate variables
    dx1 = Lx1 / Nx1;
    dx2 = Lx2 / Nx2;
    dx3 = Lx3 / Nx3;
    
    plotSpacingI = Ni / numOfFrames;
    if (plotSpacingI == 0){
        plotSpacingI = 1;
    }
    if (plottedResolutionX1 == 0 || Nx1 < plottedResolutionX1) {
        plotSpacingX1 = 1;
    } else {
        plotSpacingX1 = Nx1 / plottedResolutionX1;
    }
    if (plottedResolutionX2 == 0 || Nx2 < plottedResolutionX2){
        plotSpacingX2 = 1;
    } else {
        plotSpacingX2 = Nx2 / plottedResolutionX2;
    }
    if (plottedResolutionX3 == 0 || Nx3 < plottedResolutionX3){
        plotSpacingX3 = 1;
    } else {
        plotSpacingX3 = Nx3 / plottedResolutionX3;
    }
    
    // makes sure plotSpacing values are valid
    if (plotSpacingX1 == 0) {
        plotSpacingX1 = 1;
    } else {
        while (Nx1 % plotSpacingX1 != 0){
            plotSpacingX1--;
        }
    }
    if (plotSpacingX2 == 0) {
        plotSpacingX2 = 1;
    } else {
        while (Nx2 % plotSpacingX2 != 0){
            plotSpacingX2--;
        }
    }
    if (plotSpacingX2 == 0) {
        plotSpacingX2 = 1;
    } else {
        while (Nx1 % plotSpacingX2 != 0){
            plotSpacingX2--;
        }
    }
    
    dt = hbar/(2 * hbar * hbar/(m*dx1*dx1) + Vmax) * 1;
    if (numOfDim ==2){
        dt = hbar/(2 * hbar * hbar/(m*dx1*dx1) + 2 * hbar * hbar/(m*dx2*dx2) + Vmax) * 1;
    }
}

// gets values from the 'situation' textfile with special format
string Schrodinger::getValue(ifstream &situationFile){
    char c;
    string value;
    while (situationFile.get(c)){
        if (c == ':') {
            situationFile.get(c);
            while (c != '#') {
                if (c == '"'){
                    value = "";
                    situationFile.get(c);
                    while (c != '"'){
                        value += c;
                        situationFile.get(c);
                    }
                    return value;
                }
                value += c;
                situationFile.get(c);
            }
            return value;
        }
    }
    return value;
}

void Schrodinger::makeInitState(){
    if (numOfDim == 1){
        makeInitState1D();
    } else if (numOfDim == 2){
        makeInitState2D();
    } else {
        makeInitState3D();
    }
    normalizePsi();
}

void Schrodinger::makeInitState1D(){
    if (probDistrb == "sinusoidalGaussian") {
        double Lx1OverSDx1 = stod(getValue(*situationFile));
        double SDx1 = Lx1/Lx1OverSDx1;
        double p = stod(getValue(*situationFile));
        double k = p / hbar;
        double startX1 = Nx1 / 4;
        for (int x1 = 0; x1 < Nx1; x1++){
            psi_r1[x1] = exp(-pow(dx1 * (x1 - startX1) / SDx1, 2) / 2) * cos(k * dx1 * x1);
            psi_r2[x1] = 0;
            psi_r3[x1] = 0;
            psi_i1[x1] = exp(-pow(dx1 * (x1 - startX1) / SDx1, 2) / 2) * sin(k * dx1 * x1);
            psi_i2[x1] = 0;
            psi_i3[x1] = 0;
        }
    }
}

void Schrodinger::makeInitState2D(){
    if (probDistrb == "sinusiodalGaussian") {
        double Lx1OverSDx1 = stod(getValue(*situationFile));
        double Lx2OverSDx2 = stod(getValue(*situationFile));
        double p = stod(getValue(*situationFile));
        double SDx1 = Lx1 / Lx1OverSDx1;
        double SDx2 = Lx2 / Lx2OverSDx2;
        double k = p / hbar;
        int startX1 = Nx1 / 4;
        int startX2 = Nx2 / 2;
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                psi_r1[Nx1*x2 + x1] = exp(-pow(dx1 * (x1 - startX1) / SDx1, 2) / 2 - pow(dx2 * (x2 - startX2) / SDx2, 2) / 2) * cos(k * (dx1 * x1));
                psi_r2[Nx1*x2 + x1] = 0;
                psi_r3[Nx1*x2 + x1] = 0;
                psi_i1[Nx1*x2 + x1] = exp(-pow(dx1 * (x1 - startX1) / SDx1, 2) / 2 - pow(dx2 * (x2 - startX2) / SDx2, 2) / 2) * sin(k * (dx1 * x1));
                psi_i2[Nx1*x2 + x1] = 0;
                psi_i3[Nx1*x2 + x1] = 0;
            }
        }
    }
}

void Schrodinger::makeInitState3D(){
    
}

void Schrodinger::normalizePsi(){
    double probability = findProbability();
    for (int x3 = 0; x3 < Nx3; x3++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] /= sqrt(probability);
                psi_i1[Nx1*Nx2*x3+Nx1*x2+x1] /= sqrt(probability);
            }
        }
    }
}

double Schrodinger::findProbability(){
    double probability = 0.0;
    for (int x3 = 0; x3 < Nx3; x3++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                probability += (psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] * psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] + psi_i1[Nx1*Nx2*x3+Nx1*x2+x1] * psi_i1[Nx1*Nx2*x3+Nx1*x2+x1]) * dx1 * dx2 * dx3;
            }
        }
    }
    return probability;
}

void Schrodinger::setV(){
    if (numOfDim == 1){
        setV1D();
    } else if (numOfDim == 2){
        setV2D();
    } else {
        setV3D();
    }
}

void Schrodinger::setV1D(){
    if (potential == "free"){
        setVtoZero();
        Vmax = 0;
    } else if (potential == "constBarrier"){
        double V0 = stod(getValue(*situationFile));
        double VThickness = stod(getValue(*situationFile));
        setVtoZero();
        for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
            V[x1] = V0;
        }
        Vmax = V0;
    } else if (potential == "triangle"){
        double V0 = stod(getValue(*situationFile));
        double VThickness = stod(getValue(*situationFile));
        double VdistanceToMax = stod(getValue(*situationFile)); // must be between 0 and VThickness
        
        setVtoZero();
        for (int x1 = (Nx1/2); x1 < Nx1/2 + VdistanceToMax/dx1; x1++){
            V[x1] = V0/(VdistanceToMax/dx1)*(x1-Nx1/2);
        }
        for (int x1 = (Nx1/2) + VdistanceToMax/dx1; x1 < Nx1/2 + VThickness/dx1; x1++){
            V[x1] = V0 - V0/((VThickness-VdistanceToMax)/dx1)* (x1-Nx1/2 - VdistanceToMax/dx1);
        }
        Vmax = V0;
    }
}

void Schrodinger::setV2D(){
    if (potential == "free"){
        setVtoZero();
        Vmax = 0;
    } else if (potential == "constBarrier"){
        double V0 = stod(getValue(*situationFile));
        double VThickness = stod(getValue(*situationFile));
        setVtoZero();
        for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
            for (int x2 = 0; x2 < Nx2; x2++){
                V[Nx1*x2+x1] = V0;
            }
        }
        Vmax = V0;
    } else if (potential == "multiSlit"){
        double V0 = stod(getValue(*situationFile));
        double VThickness = stod(getValue(*situationFile));
        int slitNumber = stoi(getValue(*situationFile));
        double slitWidth = stod(getValue(*situationFile));
        double slitDistance = stod(getValue(*situationFile));
        
        setVtoZero();
        // Making constant potential barrier
        for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
            for (int x2 = 0; x2 < Nx2; x2++){
                V[Nx1*x2+x1] = V0;
            }
        }
        // Making slits in constant potential barrier
        int slitsPlaced = 0;
        double nextSlitX2 = 0; // to keep track of where to place next slit
        // Placing first slit, if odd number of slits.
        if (slitNumber % 2 == 0){
            nextSlitX2 = slitDistance/(dx2*2);
        }else if (slitsPlaced % 2 == 1){
            for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                for (int x2 = Nx2/2 - slitWidth/(dx2*2); x2 < Nx2/2 + slitWidth/(dx2*2); x2++){
                    V[Nx1*x2+x1] = 0;
                }
            }
            slitsPlaced++;
            nextSlitX2 = slitWidth/dx2/2 + slitDistance/dx2;
        }
        // Placing remaining slits, two at a time
        while (slitsPlaced < slitNumber && (nextSlitX2 + slitWidth < Lx2/2)){
            for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                for (int x2 = nextSlitX2/dx2; x2 < nextSlitX2/dx2 + slitWidth/dx2; x2++){
                    V[Nx1*(Nx2/2 + x2) +x1] = 0;
                    V[Nx1*(Nx2/2 - x2) +x1] = 0;
                }
            }
            slitsPlaced += 2;
            nextSlitX2 += slitWidth/dx2 + slitDistance/dx2;
        }
        Vmax = V0;

    } else if (potential == "circle"){
        double V0 = stod(getValue(*situationFile));
        int radius = static_cast<int>(round(stod(getValue(*situationFile)) / dx1));
        int centerX1 = static_cast<int>(round(stod(getValue(*situationFile)) / dx1));
        int centerX2 = static_cast<int>(round(stod(getValue(*situationFile)) / dx2));
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                if (pow(x1-centerX1,2)+pow(x2-centerX2,2)<pow(radius,2)) {
                    V[Nx1*x2 + x1] = V0;
                } else {
                    V[Nx1*x2 + x1] = 0;
                }
            }
        }
        Vmax = V0;
    } else if (potential == "ball"){
        double V0 = stod(getValue(*situationFile));
        int radius = static_cast<int>(round(stod(getValue(*situationFile)) / dx1));
        int centerX1 = static_cast<int>(round(stod(getValue(*situationFile)) / dx1));
        int centerX2 = static_cast<int>(round(stod(getValue(*situationFile)) / dx2));
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                if (pow(x1-centerX1,2)+pow(x2-centerX2,2)<pow(radius,2)) {
                    V[Nx1*x2 + x1] = V0*(pow(radius,2)-(pow(x1-centerX1,2)+pow(x2-centerX2,2)));
                } else {
                    V[Nx1*x2 + x1] = 0;
                }
            }
        }
        Vmax = V0;
    }
}

void Schrodinger::setV3D(){
    
}

void Schrodinger::setVtoZero(){
    for (int x3 = 0; x3 < Nx3; x3++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                V[Nx1*Nx2*x3+Nx1*x2+x1] = 0;
            }
        }
    }
}

double Schrodinger::findEnergy(){
    if (numOfDim == 1) {
        return findEnergy1D();
    } else if (numOfDim == 2){
        return findEnergy2D();
    } else {
        return findEnergy3D();
    }
}

double Schrodinger::findEnergy1D(){
    double cx1 = -hbar * hbar / 2 / m / dx1 / dx1;
    double energy = 0.0;
    int i;
    for (int x1 = 1; x1 < Nx1 - 1; x1++){
        i = x1;
        energy += (cx1*(psi_r1[i]*(psi_r1[i+1]+psi_r1[i-1]) + psi_i1[i]*(psi_i1[i+1]+psi_i1[i-1])) + (V[i]-2*cx1)*(psi_r1[i]*psi_r1[i] + psi_i1[i]*psi_i1[i]))*dx1;
    }
    return energy;
}

double Schrodinger::findEnergy2D(){
    double cx1 = -hbar * hbar / 2 / m / dx1 / dx1;
    double cx2 = -hbar * hbar / 2 / m / dx2 / dx2;
    double energy = 0.0;
    int i;
    for (int x2 = 1; x2 < Nx2 - 1; x2++){
        for (int x1 = 1; x1 < Nx1 - 1; x1++){
            i = Nx1*x2 + x1;
            energy += (cx1*(psi_r1[i]*(psi_r1[i+1]+psi_r1[i-1]) + psi_i1[i]*(psi_i1[i+1]+psi_i1[i-1])) + cx2*(psi_r1[i]*(psi_r1[i+Nx1]+psi_r1[i-Nx1]) + psi_i1[i]*(psi_i1[i+Nx1]+psi_i1[i-Nx1])) + (V[i]-2*cx1-2*cx2)*(psi_r1[i]*psi_r1[i] + psi_i1[i]*psi_i1[i]))*dx1*dx2;
        }
    }
    return energy;
}

double Schrodinger::findEnergy3D(){
    double cx1 = -hbar * hbar / 2 / m / dx1 / dx1;
    double cx2 = -hbar * hbar / 2 / m / dx2 / dx2;
    double cx3 = -hbar * hbar / 2 / m / dx3 / dx3;
    double energy = 0.0;
    int i;
    for (int x3 = 1; x3 < Nx3 - 1; x3++){
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*Nx2*x3 + Nx1*x2 + x1;
                energy += (cx1*(psi_r1[i]*(psi_r1[i+1]+psi_r1[i-1]) + psi_i1[i]*(psi_i1[i+1]+psi_i1[i-1])) + cx2*(psi_r1[i]*(psi_r1[i+Nx1]+psi_r1[i-Nx1]) + psi_i1[i]*(psi_i1[i+Nx1]+psi_i1[i-Nx1]))  + cx3*(psi_r1[i]*(psi_r1[i+Nx1*Nx2]+psi_r1[i-Nx1*Nx2]) + psi_i1[i]*(psi_i1[i+Nx1*Nx2]+psi_i1[i-Nx1*Nx2])) + (V[i]-2*cx1-2*cx2-2*cx3)*(psi_r1[i]*psi_r1[i] + psi_i1[i]*psi_i1[i]))*dx1*dx2*dx3;
            }
        }
    }
    return energy;
}

void Schrodinger::checkForDirectory(bool newSimulation){
    char fileOpenType[] = "wb";
    if (!newSimulation) {
        fileOpenType[0] = 'a';
        fileOpenType[1] = 'b';
    }
    FILE* plotProbabilityFile = fopen((filename + "/" + filename + "_plot_probability").c_str(), fileOpenType);
    if (plotProbabilityFile != 0){
        filename = filename + "/" + filename;
    }
    fclose(plotProbabilityFile);
}

void Schrodinger::finiteDifference(bool newSimulation){
    char fileOpenType[] = "wb";
    if (!newSimulation) {
        fileOpenType[0] = 'a';
        fileOpenType[1] = 'b';
    }
    if (numOfDim == 1){
        finiteDifference1D(fileOpenType);
    } else if (numOfDim == 2){
        finiteDifference2D(fileOpenType);
    } else {
        finiteDifference3D(fileOpenType);
    }
}

void Schrodinger::finiteDifference1D(char* fileOpenType){
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), fileOpenType);
    FILE* plotPsiRFile = fopen((filename + "_plot_psi_r").c_str(), fileOpenType);
    FILE* plotPsiIFile = fopen((filename + "_plot_psi_i").c_str(), fileOpenType);
    double c1 = hbar * dt / m / dx1 / dx1;
    double c2 = 2 * dt / hbar;
    if (c1 == 0){
        cout << "c1 is 0 and the wavefunction will not change. You will probably need to make dt bigger." << endl;
    }
    for (int x = 1; x < Nx1 - 1; x++){
        psi_r2[x] = psi_r1[x] + (c1 + c2 / 2 * V[x]) * psi_i1[x] - c1 / 2 * (psi_i1[x + 1] + psi_i1[x - 1]);
        psi_i2[x] = psi_i1[x] - (c1 + c2 / 2 * V[x]) * psi_r1[x] + c1 / 2 * (psi_r1[x + 1] + psi_r1[x - 1]);
    }
    for (int t = 0; t < Ni; t++){
        double possibility;
        if (t % plotSpacingI == 0){;
            for (int x = 0; x < Nx1; x += plotSpacingX1){
                fwrite(&psi_r1[x], sizeof(double), 1, plotPsiRFile);
                fwrite(&psi_i1[x], sizeof(double), 1, plotPsiIFile);
                possibility = psi_r1[x] * psi_r1[x] + psi_i1[x] * psi_i1[x];
                fwrite(&possibility, sizeof(double), 1, plotProbabilityFile);
            }
        }
        for (int x = 1; x < Nx1 - 1; x++){
            psi_r3[x] = psi_r1[x] + (2 * c1 + c2 * V[x]) * psi_i2[x] - c1 * (psi_i2[x + 1] + psi_i2[x - 1]);
            psi_i3[x] = psi_i1[x] - (2 * c1 + c2 * V[x]) * psi_r2[x] + c1 * (psi_r2[x + 1] + psi_r2[x - 1]);
        }
        for (int x = 1; x < Nx1 - 1; x++){
            psi_r1[x] = psi_r2[x] + (2 * c1 + c2 * V[x]) * psi_i3[x] - c1 * (psi_i3[x + 1] + psi_i3[x - 1]);
            psi_i1[x] = psi_i2[x] - (2 * c1 + c2 * V[x]) * psi_r3[x] + c1 * (psi_r3[x + 1] + psi_r3[x - 1]);
        }
        for (int x = 1; x < Nx1 - 1; x++){
            psi_r2[x] = psi_r3[x] + (2 * c1 + c2 * V[x]) * psi_i1[x] - c1 * (psi_i1[x + 1] + psi_i1[x - 1]);
            psi_i2[x] = psi_i3[x] - (2 * c1 + c2 * V[x]) * psi_r1[x] + c1 * (psi_r1[x + 1] + psi_r1[x - 1]);
        }
    }
    fclose(plotProbabilityFile);
    fclose(plotPsiRFile);
    fclose(plotPsiIFile);
}

void Schrodinger::finiteDifference2D(char* fileOpenType){
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), fileOpenType);
    double c1x1 = hbar * dt / m / dx1 / dx1;
    double c1x2 = hbar * dt / m / dx2 / dx2;
    double c2 = 2 * dt / hbar;
    if (c1x1 == 0 || c1x2 == 0){
        cout << "c1x1 or c1x2 is 0 and the wavefunction will not change. You will probably need to make dt bigger." << endl;
    }
    int i;
    for (int x2 = 1; x2 < Nx2 - 1; x2++){
        for (int x1 = 1; x1 < Nx1 - 1; x1++){
            i = Nx1*x2+x1;
            psi_r2[i] = psi_r1[i] + (c1x1 + c1x2 + c2/2*V[i])*psi_i1[i] - c1x1/2*(psi_i1[i+1] + psi_i1[i-1]) - c1x2/2*(psi_i1[i+Nx1] + psi_i1[i-Nx1]);
            psi_i2[i] = psi_i1[i] - (c1x1 + c1x2 + c2/2*V[i])*psi_r1[i] + c1x1/2*(psi_r1[i+1] + psi_r1[i-1]) + c1x2/2*(psi_r1[i+Nx1] + psi_r1[i-Nx1]);
        }
    }
    for (int t = 0; t < Ni; t++){
        double possibility;
        if (t % plotSpacingI == 0){
            for (int x2 = 0; x2 < Nx2; x2+=plotSpacingX2){
                for (int x1 = 0; x1 < Nx1; x1+=plotSpacingX1){
                    possibility = psi_r1[x2*Nx1 + x1] * psi_r1[x2*Nx1 + x1] + psi_i1[x2*Nx1 + x1] * psi_i1[x2*Nx1 + x1];
                    fwrite(&possibility, sizeof(double), 1, plotProbabilityFile);
                }
            }
        }
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*x2+x1;
                psi_r3[i] = psi_r1[i] + (2*c1x1 + 2*c1x2 + c2*V[i])*psi_i2[i] - c1x1*(psi_i2[i+1] + psi_i2[i-1]) - c1x2*(psi_i2[i+Nx1] + psi_i2[i-Nx1]);
                psi_i3[i] = psi_i1[i] - (2*c1x1 + 2*c1x2 + c2*V[i])*psi_r2[i] + c1x1*(psi_r2[i+1] + psi_r2[i-1]) + c1x2*(psi_r2[i+Nx1] + psi_r2[i-Nx1]);
            }
        }
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*x2+x1;
                psi_r1[i] = psi_r2[i] + (2*c1x1 + 2*c1x2 + c2*V[i])*psi_i3[i] - c1x1*(psi_i3[i+1] + psi_i3[i-1]) - c1x2*(psi_i3[i+Nx1] + psi_i3[i-Nx1]);
                psi_i1[i] = psi_i2[i] - (2*c1x1 + 2*c1x2 + c2*V[i])*psi_r3[i] + c1x1*(psi_r3[i+1] + psi_r3[i-1]) + c1x2*(psi_r3[i+Nx1] + psi_r3[i-Nx1]);
            }
        }
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*x2+x1;
                psi_r2[i] = psi_r3[i] + (2*c1x1 + 2*c1x2 + c2*V[i])*psi_i1[i] - c1x1*(psi_i1[i+1] + psi_i1[i-1]) - c1x2*(psi_i1[i+Nx1] + psi_i1[i-Nx1]);
                psi_i2[i] = psi_i3[i] - (2*c1x1 + 2*c1x2 + c2*V[i])*psi_r1[i] + c1x1*(psi_r1[i+1] + psi_r1[i-1]) + c1x2*(psi_r1[i+Nx1] + psi_r1[i-Nx1]);
            }
        }
    }
    fclose(plotProbabilityFile);
}

void Schrodinger::finiteDifference3D(char* fileOpenType){
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), fileOpenType);
    double c1x1 = hbar * dt / m / dx1 / dx1;
    double c1x2 = hbar * dt / m / dx2 / dx2;
    double c1x3 = hbar * dt / m / dx3 / dx3;
    double c2 = 2 * dt / hbar;
    int i;
    for (int x3 = 1; x3 < Nx3 - 1; x3++){
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*Nx2*x3+Nx1*x2+x1;
                psi_r2[i] = psi_r1[i] + (c1x1 + c1x2 + c1x3 + c2/2*V[i])*psi_i1[i] - c1x1/2*(psi_i1[i+1] + psi_i1[i-1]) - c1x2/2*(psi_i1[i+Nx1] + psi_i1[i-Nx1]) - c1x3/2*(psi_i1[i+Nx1*Nx2] + psi_i1[i-Nx1*Nx2]);
                psi_i2[i] = psi_i1[i] - (c1x1 + c1x2 + c1x3 + c2/2*V[i])*psi_r1[i] + c1x1/2*(psi_r1[i+1] + psi_r1[i-1]) + c1x2/2*(psi_r1[i+Nx1] + psi_r1[i-Nx1]) + c1x3/2*(psi_r1[i+Nx1*Nx2] + psi_r1[i-Nx1*Nx2]);
            }
        }
    }
    for (int t = 0; t < Ni; t++){
        double possibility;
        if (t % plotSpacingI == 0){
            for (int x3 = 0; x3 < Nx3; x3+= plotSpacingX3){
                for (int x2 = 0; x2 < Nx2; x2+=plotSpacingX2){
                    for (int x1 = 0; x1 < Nx1; x1+=plotSpacingX1){
                        possibility = psi_r1[Nx1*Nx2*x3 + Nx1*x2 + x1] * psi_r1[Nx1*Nx2*x3 + Nx1*x2 + x1] + psi_i1[Nx1*Nx2*x3 + Nx1*x2 + x1] * psi_i1[Nx1*Nx2*x3 + Nx1*x2 + x1];
                        fwrite(&possibility, sizeof(double), 1, plotProbabilityFile);
                    }
                }
            }
        }
        for (int x3 = 1; x3 < Nx3 - 1; x3++){
            for (int x2 = 1; x2 < Nx2 - 1; x2++){
                for (int x1 = 1; x1 < Nx1 - 1; x1++){
                    i = Nx1*Nx2*x3+Nx1*x2+x1;
                    psi_r3[i] = psi_r1[i] + (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_i2[i] - c1x1*(psi_i2[i+1] + psi_i2[i-1]) - c1x2*(psi_i2[i+Nx1] + psi_i2[i-Nx1]) - c1x3*(psi_i2[i+Nx1*Nx2] + psi_i2[i-Nx1*Nx2]);
                    psi_i3[i] = psi_i1[i] - (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_r2[i] + c1x1*(psi_r2[i+1] + psi_r2[i-1]) + c1x2*(psi_r2[i+Nx1] + psi_r2[i-Nx1]) + c1x3*(psi_r2[i+Nx1*Nx2] + psi_r2[i-Nx1*Nx2]);
                }
            }
        }
        for (int x3 = 1; x3 < Nx3 - 1; x3++){
            for (int x2 = 1; x2 < Nx2 - 1; x2++){
                for (int x1 = 1; x1 < Nx1 - 1; x1++){
                    i = Nx1*Nx2*x3+Nx1*x2+x1;
                    psi_r1[i] = psi_r2[i] + (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_i3[i] - c1x1*(psi_i3[i+1] + psi_i3[i-1]) - c1x2*(psi_i3[i+Nx1] + psi_i3[i-Nx1]) - c1x3*(psi_i3[i+Nx1*Nx2] + psi_i3[i-Nx1*Nx2]);
                    psi_i1[i] = psi_i2[i] - (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_r3[i] + c1x1*(psi_r3[i+1] + psi_r3[i-1]) + c1x2*(psi_r3[i+Nx1] + psi_r3[i-Nx1]) + c1x3*(psi_r3[i+Nx1*Nx2] + psi_r3[i-Nx1*Nx2]);
                }
            }
        }
        for (int x3 = 1; x3 < Nx3 - 1; x3++){
            for (int x2 = 1; x2 < Nx2 - 1; x2++){
                for (int x1 = 1; x1 < Nx1 - 1; x1++){
                    i = Nx1*Nx2*x3+Nx1*x2+x1;
                    psi_r2[i] = psi_r3[i] + (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_i1[i] - c1x1*(psi_i1[i+1] + psi_i1[i-1]) - c1x2*(psi_i1[i+Nx1] + psi_i1[i-Nx1]) - c1x3*(psi_i1[i+Nx1*Nx2] + psi_i1[i-Nx1*Nx2]);
                    psi_i2[i] = psi_i3[i] - (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_r1[i] + c1x1*(psi_r1[i+1] + psi_r1[i-1]) + c1x2*(psi_r1[i+Nx1] + psi_r1[i-Nx1]) + c1x3*(psi_r1[i+Nx1*Nx2] + psi_r1[i-Nx1*Nx2]);
                }
            }
        }
    }
    fclose(plotProbabilityFile);
}

void Schrodinger::storeFinalState(long time){
    FILE* finalStateFile = fopen((filename + "_finalState").c_str(), "wb");
    fwrite(&psi_r1[0], sizeof(double), Nx1*Nx2*Nx3, finalStateFile);
    fwrite(&psi_i1[0], sizeof(double), Nx1*Nx2*Nx3, finalStateFile);
    fclose(finalStateFile);
    FILE* potentialFile = fopen((filename + "_potential").c_str(), "wb");
    fwrite(&V[0], sizeof(double), Nx1*Nx2*Nx3, potentialFile);
    fclose(potentialFile);
    finalProb = findProbability();
    finalEnergy = findEnergy();
    ofstream usefulSimValues;
    usefulSimValues.open(filename + "_simulationValues.txt");
    usefulSimValues << startEnergy << endl << finalEnergy << endl << finalProb << endl << Vmax << endl << time << endl << plotSpacingI << endl << plotSpacingX1 << endl << plotSpacingX2 << endl << plotSpacingX3 << endl;
    usefulSimValues.close();
}

void Schrodinger::loadFinalState(){
    FILE* finalStateFile = fopen((filename + "_finalState").c_str(), "rb");
    fread(psi_r1, sizeof(double), Nx1*Nx2*Nx3, finalStateFile);
    fread(psi_i1, sizeof(double), Nx1*Nx2*Nx3, finalStateFile);
    fclose(finalStateFile);
}