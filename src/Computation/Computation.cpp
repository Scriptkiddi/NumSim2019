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
    int nVelo = settings_.nVelo; 


    //geometry_ = settings_.geometry;

    //initialize initial values
    fInit = settings_.fInit; 
    
    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] / (settings_.nCells[0]);
    meshWidth_[1] = settings_.physicalSize[1] / (settings_.nCells[1]);

    //initialize staggeredGrid
    StaggeredGrid grid(nCellsBoundary, nVelo, meshWidth_);
    staggeredGrid_ = make_shared<StaggeredGrid>(grid);

    //initialize outputWriters
    //OutputWriterText outText(staggeredGrid_);
    //outputWriterText_ = make_unique<OutputWriterText>(outText);

    OutputWriterParaview outPara(staggeredGrid_);
    outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);
}

void Computation::runSimulation() {
    double t = 0;
    applyInitialConditions();
    applyWeights(); //w
    applyBoundaryValuesF();

    double lastOutputTime = 0;
    for (int timeStepNumber = 0;
         std::abs(t - settings_.endTime) > 1e-10 && settings_.endTime - t > 0; timeStepNumber++) {
        
        computeTimeStepWidth();

        if (t + dt_ > settings_.endTime) {
            dt_ = settings_.endTime - t;
        }
        //outputWriterText_->writeFile(t);

        applyLatticeVelocities(); //c_i
        cout << "before compMacro -- current time: " << t << " dt: " << dt_ << ", rho = " << staggeredGrid_.get()->rho(3,3) << ", u = " << staggeredGrid_.get()->u(3,3) <<
        ", v = " << staggeredGrid_.get()->v(3,3) << ", p = " << staggeredGrid_.get()->p(3,3) << endl;
        
        computeMacroscopicQuantities(timeStepNumber); //DensityPressureAndVelocities
        cout << "after compMacro -- current time: " << t << " dt: " << dt_ << ", rho = " << staggeredGrid_.get()->rho(3,3) << ", u = " << staggeredGrid_.get()->u(3,3) <<
        ", v = " << staggeredGrid_.get()->v(3,3) << ", p = " << staggeredGrid_.get()->p(3,3) << endl;
        computeFtempFeq(timeStepNumber); //Collision step
        applyBoundaryValuesF();
        computeF(timeStepNumber); //Streaming step

        t += dt_;
        if (t - lastOutputTime > settings_.outputFileEveryDt - 1e-4) {
            //cout << "current time: " << t << " dt: " << dt_ << endl;
            outputWriterParaview_->writeFile(t);
            lastOutputTime = t;
        }
    }

    if (std::fabs(t - lastOutputTime) > 1e-4) {
        outputWriterParaview_->writeFile(t);
        lastOutputTime = t;
    }

}

void Computation::applyInitialConditions() { // für f, tau
    fInit = 0.1338;
    settings_.rhoInit = fInit * 9;
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
            
                staggeredGrid_.get()->f(i, j, k) = fInit;
            }
        }
    }
    for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
        staggeredGrid_.get()->f(3, 3, k) = .2;
    }
    tau_ = settings_.viscosity/(settings_.rhoInit * pow(settings_.cs,2)) + 0.5; 
    //std::cout << "tau: " << tau_ << endl;
}

void Computation::applyLatticeVelocities(){
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
        //std::cout << "c nach korrektur durch konstanten term: " << staggeredGrid_.get()->c(k,0) << endl;
    }
}

void Computation::applyWeights(){ // Für 3D anders!
    //if (settings_.nVelo == 9){
    staggeredGrid_.get()->w(0) = 4.0/9.0;
    staggeredGrid_.get()->w(1) = 1.0/9.0;
    staggeredGrid_.get()->w(2) = 1.0/36.0;
    staggeredGrid_.get()->w(3) = 1.0/9.0;
    staggeredGrid_.get()->w(4) = 1.0/36.0;
    staggeredGrid_.get()->w(5) = 1.0/9.0;
    staggeredGrid_.get()->w(6) = 1.0/36.0;
    staggeredGrid_.get()->w(7) = 1.0/9.0;
    staggeredGrid_.get()->w(8) = 1.0/36.0;
    //} else {
    //    std::cout << "Please choose D2Q9." << std::endl;
    //}
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
           staggeredGrid_.get()->ftmp(i, j_low, 8) = staggeredGrid_.get()->ftmp(i - 1, j_low + 1, 4);
        }
        if(i >= staggeredGrid_.get()->iBegin() && i <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i, j_low, 1) = staggeredGrid_.get()->ftmp(i, j_low + 1, 5);
        }
        if(i < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i ,j_low, 2) = staggeredGrid_.get()->ftmp(i + 1, j_low + 1, 6);
        }
    }

    //top bound
    j_high = staggeredGrid_.get()->jEnd() + 1; //Schleife kürzer
    for (int i = staggeredGrid_.get()->iBegin() - 1; i <= staggeredGrid_.get()->iEnd() - 1; i++) {
        if(i > staggeredGrid_.get()->iBegin()){
           staggeredGrid_.get()->ftmp(i, j_high, 6) = staggeredGrid_.get()->ftmp(i - 1, j_high - 1, 2);
        }
        if(i >= staggeredGrid_.get()->iBegin() && i <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i, j_high, 5) = staggeredGrid_.get()->ftmp(i, j_high - 1, 1);
        }
        if(i < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i ,j_high, 4) = staggeredGrid_.get()->ftmp(i + 1, j_high - 1, 8);
        }
    }

    //left bound
    i_low = staggeredGrid_.get()->iBegin() - 1;
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        if(j > staggeredGrid_.get()->jBegin()){
           staggeredGrid_.get()->ftmp(i_low, j, 4) = staggeredGrid_.get()->ftmp(i_low + 1, j - 1, 8);
        }
        if(j >= staggeredGrid_.get()->iBegin() && j <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i_low, j, 3) = staggeredGrid_.get()->ftmp(i_low + 1, j, 7);
        }
        if(j < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i_low, j, 2) = staggeredGrid_.get()->ftmp(i_low + 1, j + 1, 6);
        }
    }

    //right bound
    i_high = staggeredGrid_.get()->iEnd() + 1;
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        if(j > staggeredGrid_.get()->jBegin()){
           staggeredGrid_.get()->ftmp(i_high, j, 6) = staggeredGrid_.get()->ftmp(i_high - 1, j - 1, 2);
        }
        if(j >= staggeredGrid_.get()->iBegin() && j <= staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i_high, j, 7) = staggeredGrid_.get()->ftmp(i_high - 1, j, 3);
        }
        if(j < staggeredGrid_.get()->iEnd()){
            staggeredGrid_.get()->ftmp(i_high, j, 8) = staggeredGrid_.get()->ftmp(i_high - 1, j + 1, 4);
        }
    }
}

void Computation::computeTimeStepWidth() {
    
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
    /*
    double condition_diffusion = pow(staggeredGrid_.get()->dx(), 2) * pow(staggeredGrid_.get()->dy(), 2) /
                                 (pow(staggeredGrid_.get()->dx(), 2) + pow(staggeredGrid_.get()->dy(), 2)) *
                                 settings_.re /
                                 2;
    */   
    double condition_convection1 = staggeredGrid_.get()->dx() / uMaximum;
    double condition_convection2 = staggeredGrid_.get()->dy() / vMaximum;
    /*
    // Min time requirements for temperature
    double condition_temp = settings_.re * settings_.prandtl * 0.5 *
                            pow((1 / pow(staggeredGrid_.get()->dx(), 2) + 1 / pow(staggeredGrid_.get()->dy(), 2)),
                                -1);
    */
    dt_ = min(condition_convection1, condition_convection2);
    //dt_ = min(condition_diffusion, dt_);
    //dt_ = min(condition_temp, dt_);
    dt_ = min(settings_.maximumDt, dt_) * settings_.timeStepRelaxation;
    
}

void Computation::computeMacroscopicQuantities(int t){ //DensityPressureAndVelocities
    double fSum = 0.0;
    double fSumWeightedX = 0.0;
    double fSumWeightedY = 0.0;
            
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++){
                fSum += staggeredGrid_.get()->f(i,j,k);
                /*
                if (t<=1 && i == 3 && j == 3){
                    std::cout << "f: " << staggeredGrid_.get()->f(i,j,k) << ", fSum so far: " << fSum << endl;
                }
                */
                fSumWeightedX += staggeredGrid_.get()->f(i,j,k) * staggeredGrid_.get()->c(k,0);
                fSumWeightedY += staggeredGrid_.get()->f(i,j,k) * staggeredGrid_.get()->c(k,1);
            }
            staggeredGrid_.get()->rho(i, j) = fSum;
            staggeredGrid_.get()->u(i, j) = fSumWeightedX;
            fSum = 0.0;
            fSumWeightedX = 0.0;
            fSumWeightedY = 0.0;

            staggeredGrid_.get()->p(i, j) = pow(settings_.cs,2) * fSum;
            
            
            if (t <= 1 && i < 5 && j < 5) {
                std::cout << "rho: " << fSum << ", u: " << fSumWeightedX << ", v:  " << fSumWeightedY << ", p: " << staggeredGrid_.get()->p(i,j) << endl;
            }
            
            }
        }
}

void Computation::computeFtempFeq(int t){
    for (int j = staggeredGrid_.get()->jBegin(); j <= staggeredGrid_.get()->jEnd(); j++) {
        for (int i = staggeredGrid_.get()->iBegin(); i <= staggeredGrid_.get()->iEnd(); i++) {
            for (int k = staggeredGrid_.get()->kBegin(); k <= staggeredGrid_.get()->kEnd(); k++) {
                //w * rho* (1 + (c0 * u  + c1 * v) / cs^2 + (c0 * u  + c1 * v)^2 / (2 * cs^4) + (u^2  + v^2) / (2 * cs^2))
                staggeredGrid_.get()->feq(i,j,k) = staggeredGrid_.get()->w(k) * staggeredGrid_.get()->rho(i,j) * ( 1 +
                (staggeredGrid_.get()->c(k,0) * staggeredGrid_.get()->u(i,j)  + staggeredGrid_.get()->c(k,1) * staggeredGrid_.get()->v(i,j)) / pow(settings_.cs,2) 
                + pow(staggeredGrid_.get()->c(k,0) * staggeredGrid_.get()->u(i,j)  + staggeredGrid_.get()->c(k,1) * staggeredGrid_.get()->v(i,j),2) / (2 * pow(settings_.cs,4))
                - (staggeredGrid_.get()->u(i,j) * staggeredGrid_.get()->u(i,j)  + staggeredGrid_.get()->v(i,j) * staggeredGrid_.get()->v(i,j)) / (2 * pow(settings_.cs,2)) 
                );

                staggeredGrid_.get()->ftmp(i,j,k) = staggeredGrid_.get()->f(i,j,k) + 1 / tau_ * (- staggeredGrid_.get()->f(i,j,k) + staggeredGrid_.get()->feq(i,j,k));
                
                if (t <= 1 && i == 3 && j == 3) {
                    
                    //std::cout << "w at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->w(k) << endl;
                    std::cout << "rho at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->rho(i,j) << endl;
                    //std::cout << "c0 at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->c(k,0) << endl;
                    //std::cout << "c1 at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->c(k,1) << endl;
                    std::cout << "u at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->u(i,j) << endl;
                    std::cout << "v at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->v(i,j) << endl;
                    //std::cout << "cs at (" << i << " - " << j << " - " << k << "): " << settings_.cs << endl;
                
                    std::cout << "feq at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->feq(i,j,k) << endl;
                    std::cout << "ftmp at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->ftmp(i,j,k) << endl;
                }
            }
        }
    }
}


void Computation::computeF(int t){
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

                if (t <= 1 && i == 3 && j == 3) {
                    //std::cout << "f at (" << i << " - " << j << " - " << k << "): " << staggeredGrid_.get()->f(i,j,k) << endl;
                }
            }
        }
    }
}