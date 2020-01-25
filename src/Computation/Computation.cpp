//
// Created by Julia Pelzer on 26.10.2019.
//

#include <memory>
#include "Computation.h"
#include "Settings.h"
#include <cmath>
#include <iostream>
#include <StaggeredGrid/StaggeredGrid.h>

using namespace std;

void Computation::initialize(int argc, char **argv) {
    std::cout << "Running with" << argv[0] << std::endl;
    Settings settings;
    settings_ = settings;
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();
    array<int, 2> nCellsBoundary = {settings_.nCells[0] + 2, settings_.nCells[1] + 2}; // Mit Ghost cells
    int nVelo = settings_.nVelo; //TODO in settings

    geometry_ = settings_.geometry;

    //initialize initial values
    fInit = settings_.fInit; //TODO in settings

    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] / (settings_.nCells[0]);
    meshWidth_[1] = settings_.physicalSize[1] / (settings_.nCells[1]);

    //initialize staggeredGrid
    StaggeredGrid grid(nCellsBoundary, nVelo, meshWidth_);
    staggeredGrid_ = make_shared<StaggeredGrid>(grid);

    //initialize outputWriters
    OutputWriterText outText(staggeredGrid_);
    outputWriterText_ = make_unique<OutputWriterText>(outText);

    OutputWriterParaview outPara(staggeredGrid_);
    outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);
}

void Computation::runSimulation() {
    double t = 0;
    applyInitialConditions();
    double lastOutputTime = 0;
    for (int timeStepNumber = 0;
         std::abs(t - settings_.endTime) > 1e-10 && settings_.endTime - t > 0; timeStepNumber++) {
        
        computeTimeStepWidth();

        if (t + dt_ > settings_.endTime) {
            dt_ = settings_.endTime - t;
        }
        outputWriterText_->writeFile(t);

        applyLatticeVelocities(); //c_i
        applyWeights();

        computeMacroscopicQuantities(); //DensityPressureAndVelocities
        computeFtempFeq(); //Collision step
        applyBoundaryValuesF();
        computeF(); //Streaming step

        t += dt_;
        if (t - lastOutputTime > settings_.outputFileEveryDt - 1e-4) {
            cout << "current time: " << t << " dt: " << dt_ << endl;
            outputWriterParaview_->writeFile(t);
            lastOutputTime = t;
        }
    }


    if (std::fabs(t - lastOutputTime) > 1e-4) {
        outputWriterParaview_->writeFile(t);
        lastOutputTime = t;
    }

}

void Computation::applyInitialConditions() {
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
            
                staggeredGrid_.get()->f(i, j, k) = fInit;
        }
    }
    }
    tau_ = settings_.viscosity/(settings_.rhoInit * pow(settings_.cs,2)) + 0.5;
}

void Computation::applyLatticeVelocities() { //c_i(k,l)
    staggeredGrid_.get()->c(0,0) = 0;
    staggeredGrid_.get()->c(0,1) = 0;

    staggeredGrid_.get()->c(1,0) = 0;
    staggeredGrid_.get()->c(1,1) = 1;

    staggeredGrid_.get()->c(2,0) = 1;
    staggeredGrid_.get()->c(2,1) = 1;

    staggeredGrid_.get()->c(3,0) = 1;
    staggeredGrid_.get()->c(3,1) = 0;

    staggeredGrid_.get()->c(4,0) = 1;
    staggeredGrid_.get()->c(4,1) = -1;

    staggeredGrid_.get()->c(5,0) = 0;
    staggeredGrid_.get()->c(5,1) = -1;

    staggeredGrid_.get()->c(6,0) = -1;
    staggeredGrid_.get()->c(6,1) = -1;

    staggeredGrid_.get()->c(7,0) = -1;
    staggeredGrid_.get()->c(7,1) = 0;

    staggeredGrid_.get()->c(8,0) = -1;
    staggeredGrid_.get()->c(8,1) = 1;

    for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
        staggeredGrid_.get()->c(k,0) = staggeredGrid_.get()->c(k,0) * staggeredGrid_.get()->dx() / dt_; //TODO *meshWidth/dt correct?
        staggeredGrid_.get()->c(k,1) = staggeredGrid_.get()->c(k,1) * staggeredGrid_.get()->dy() / dt_;
    }
}

void Computation::applyWeights(){ // Für 3D anders!
    if (settings_.nVelo == 9){
    staggeredGrid_.get()->w(0) = 4/9;
    staggeredGrid_.get()->w(1) = 1/9;
    staggeredGrid_.get()->w(2) = 1/36;
    staggeredGrid_.get()->w(3) = 1/9;
    staggeredGrid_.get()->w(4) = 1/36;
    staggeredGrid_.get()->w(5) = 1/9;
    staggeredGrid_.get()->w(6) = 1/36;
    staggeredGrid_.get()->w(7) = 1/9;
    staggeredGrid_.get()->w(8) = 1/36;
    } else {
        std::cout << "Please choose D2Q9." << std::endl;
    }
}

void Computation::applyBoundaryValuesF() { //TODO settings: rausschmeißen
// type: "NSW" (no slip wall)
    int i_low;
    int i_high;
    int j_low;
    int j_high;

    //bottom bound
    j_low = staggeredGrid_.get()->jBegin() - 1; 
    for (int i = staggeredGrid_.get()->iBegin() - 1; i <= staggeredGrid_.get()->iEnd() + 1; i++) {
        if(i > staggeredGrid_.get()->iBegin()){
           staggeredGrid_.get()->f(i, j_low, 8) = staggeredGrid_.get()->f(i - 1, j_low + 1, 4);
        }
        if(i >= staggeredGrid_.get()->iBegin() && i <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i, j_low, 1) = staggeredGrid_.get()->f(i, j_low + 1, 5);
        }
        if(i < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i ,j_low, 2) = staggeredGrid_.get()->f(i + 1, j_low + 1, 6);
        }
    }

    //top bound
    j_high = staggeredGrid_.get()->jEnd() + 1; //Schleife kürzer
    for (int i = staggeredGrid_.get()->iBegin() - 1; i <= staggeredGrid_.get()->iEnd() - 1; i++) {
        if(i > staggeredGrid_.get()->iBegin()){
           staggeredGrid_.get()->f(i, j_high, 6) = staggeredGrid_.get()->f(i - 1, j_high - 1, 2);
        }
        if(i >= staggeredGrid_.get()->iBegin() && i <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i, j_high, 5) = staggeredGrid_.get()->f(i, j_high - 1, 1);
        }
        if(i < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i ,j_high, 4) = staggeredGrid_.get()->f(i + 1, j_high - 1, 8);
        }
    }

    //left bound
    i_low = staggeredGrid_.get()->iBegin() - 1;
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        if(j > staggeredGrid_.get()->jBegin()){
           staggeredGrid_.get()->f(i_low, j, 4) = staggeredGrid_.get()->f(i_low + 1, j - 1, 8);
        }
        if(j >= staggeredGrid_.get()->iBegin() && j <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i_low, j, 3) = staggeredGrid_.get()->f(i_low + 1, j, 7);
        }
        if(j < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i_low, j, 2) = staggeredGrid_.get()->f(i_low + 1, j + 1, 6);
        }
    }

    //right bound
    i_high = staggeredGrid_.get()->iEnd() + 1;
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        if(j > staggeredGrid_.get()->jBegin()){
           staggeredGrid_.get()->f(i_high, j, 6) = staggeredGrid_.get()->f(i_high - 1, j - 1, 2);
        }
        if(j >= staggeredGrid_.get()->iBegin() && j <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i_high, j, 7) = staggeredGrid_.get()->f(i_high - 1, j, 3);
        }
        if(j < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->f(i_high, j, 8) = staggeredGrid_.get()->f(i_high - 1, j + 1, 4);
        }
    }
}

void Computation::computeTimeStepWidth() {
    /*
    double uMaximum = 0;
    for (int j = staggeredGrid_.get()->jBegin() - 1; j <= staggeredGrid_.get()->jEnd() + 1; j++) {
        for (int i = staggeredGrid_.get()->iBegin() - 1; i <= staggeredGrid_.get()->iEnd() + 1; i++) {
            if (uMaximum < abs(staggeredGrid_.get()->u(i, j))) {
                uMaximum = abs(staggeredGrid_.get()->u(i, j));
            }
        }
    }
    double vMaximum = 0;
    for (int j = staggeredGrid_.get()->jBegin() - 1; j <= staggeredGrid_.get()->jEnd() + 1; j++) {
        for (int i = staggeredGrid_.get()->iBegin() - 1; i <= staggeredGrid_.get()->iEnd() + 1; i++) {
            if (vMaximum < abs(staggeredGrid_.get()->v(i, j))) {
                vMaximum = abs(staggeredGrid_.get()->v(i, j));
            }
        }
    }
    double condition_diffusion = pow(staggeredGrid_.get()->dx(), 2) * pow(staggeredGrid_.get()->dy(), 2) /
                                 (pow(staggeredGrid_.get()->dx(), 2) + pow(staggeredGrid_.get()->dy(), 2)) *
                                 settings_.re /
                                 2;
    double condition_convection1 = staggeredGrid_.get()->dx() / uMaximum;
    double condition_convection2 = staggeredGrid_.get()->dy() / vMaximum;

    // Min time requirements for temperature
    double condition_temp = settings_.re * settings_.prandtl * 0.5 *
                            pow((1 / pow(staggeredGrid_.get()->dx(), 2) + 1 / pow(staggeredGrid_.get()->dy(), 2)),
                                -1);

    dt_ = min(condition_convection1, condition_convection2);
    dt_ = min(condition_diffusion, dt_);
    dt_ = min(condition_temp, dt_);
    dt_ = min(settings_.maximumDt, dt_) * settings_.timeStepRelaxation;
    */
}

void Computation::computeMacroscopicQuantities(){ //DensityPressureAndVelocities
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            double fSum = 0;
            double fSumWeightedX = 0;
            double fSumWeightedY = 0;
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++){
                fSum += staggeredGrid_.get()->f(i,j,k);
                fSumWeightedX += staggeredGrid_.get()->f(i,j,k) * staggeredGrid_.get()->c(k,0);
                fSumWeightedY += staggeredGrid_.get()->f(i,j,k) * staggeredGrid_.get()->c(k,1);
            }
            staggeredGrid_.get()->rho(i, j) = fSum;
            staggeredGrid_.get()->u(i, j) = fSumWeightedX;
            staggeredGrid_.get()->v(i, j) = fSumWeightedY;
            staggeredGrid_.get()->p(i, j) = pow(settings_.cs,2) * staggeredGrid_.get()->rho(i,j);
            }
        }
}

void Computation::computeFtempFeq(){
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
                //w * rho* (1 + (c0 * u  + c1 * v) / cs^2 + (c0 * u  + c1 * v)^2 / (2 * cs^4) + (u^2  + v^2) / (2 * cs^2))
                staggeredGrid_.get()->feq(i,j,k) = staggeredGrid_.get()->w(k) * staggeredGrid_.get->rho(i,j) * ( 1 +
                (staggeredGrid_.get()->c(k,0) * staggeredGrid_.get()->u(i,j)  + staggeredGrid_.get()->c(k,1) * staggeredGrid_.get()->v(i,j)) / pow(settings_.cs,2) 
                + pow(staggeredGrid_.get()->c(k,0) * staggeredGrid_.get()->u(i,j)  + staggeredGrid_.get()->c(k,1) * staggeredGrid_.get()->v(i,j),2) / (2 * pow(settings_.cs,4))
                + (staggeredGrid_.get()->u(i,j) * staggeredGrid_.get()->u(i,j)  + staggeredGrid_.get()->v(i,j) * staggeredGrid_.get()->v(i,j)) / (2 * pow(settings_.cs,2)) 
                );

                staggeredGrid_.get()->ftmp(i,j,k) = staggeredGrid_.get()->f(i,j,k) + 1 / tau_ * (staggeredGrid_.get()->f(i,j,k) - staggeredGrid_.get()->feq(i,j,k));
            }
        }
    }
}


void Computation::computeF(){
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
                staggeredGrid_.get()->f(i,j,0) = staggeredGrid_.get()->ftmp(i,j,0);
                staggeredGrid_.get()->f(i,j,1) = staggeredGrid_.get()->ftmp(i,j-1,1);
                staggeredGrid_.get()->f(i,j,2) = staggeredGrid_.get()->ftmp(i-1,j-1,2);
                staggeredGrid_.get()->f(i,j,3) = staggeredGrid_.get()->ftmp(i-1,j,3);
                staggeredGrid_.get()->f(i,j,4) = staggeredGrid_.get()->ftmp(i-1,j+1,4);
                staggeredGrid_.get()->f(i,j,5) = staggeredGrid_.get()->ftmp(i,j+1,5);
                staggeredGrid_.get()->f(i,j,6) = staggeredGrid_.get()->ftmp(i+1,j+1,6);
                staggeredGrid_.get()->f(i,j,7) = staggeredGrid_.get()->ftmp(i+1,j,7);
                staggeredGrid_.get()->f(i,j,8) = staggeredGrid_.get()->ftmp(i+1,j-1,8);
            }
        }
    }
}