//
// Created by Julia Pelzer on 26.10.2019.
//

#include <memory>
#include "Computation_solid.h"
#include "Settings.h"
#include <cmath>
#include <iostream>
#include <StaggeredGrid/CentralDifferences.h>
#include <StaggeredGrid/DonorCell.h>
#include <PressureSolver/SOR.h>
#include <TemperatureSolver/GaussSeidel.h>
#include <PressureSolver/GaussSeidel.h>

using namespace std;

void Computation_solid::initialize(int argc, char **argv) {
    std::cout << "Running with" << argv[0] << std::endl;
    Settings settings;
    settings_ = settings;
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();
    array<int, 2> nCellsBoundary = {settings_.nCells[0] + 2, settings_.nCells[1] + 2}; // Mit Ghost cells

    geometry_ = settings_.geometry;

    //initialize initial values
    tInit = settings_.tInit;

    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] / (settings_.nCells[0]);
    meshWidth_[1] = settings_.physicalSize[1] / (settings_.nCells[1]);

    //initialize discretization
    if (!settings_.useDonorCell) {
        CentralDifferences grid(nCellsBoundary, meshWidth_, settings_.gamma);
        discretization_ = make_shared<CentralDifferences>(grid);
    } else {
        DonorCell grid(nCellsBoundary, meshWidth_, settings_.alpha, settings_.gamma);
        discretization_ = make_shared<DonorCell>(grid);
    }

    //initialize IMPLICIT temperatureSolver
    GaussSeidel tSolver(discretization_, geometry_, settings_.epsilon, settings_.maximumNumberOfIterations,dt_, settings_.heatDiffusivity);
    temperatureSolver_ = make_unique<GaussSeidel>(tSolver); //TODO not working 

    //initialize outputWriters
    OutputWriterText outText(discretization_);
    outputWriterText_ = make_unique<OutputWriterText>(outText);

    OutputWriterParaview outPara(discretization_);
    outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);
}

void Computation_solid::runSimulation() {
    double t = 0;
    Array2D tTmp(discretization_.get()->t().size()); //previous temperature field (last timestep)
    applyInitialConditions();
    double lastOutputTime = 0;
    for (int timeStepNumber = 0;
        std::abs(t - settings_.endTime) > 1e-10 && settings_.endTime - t > 0; timeStepNumber++) {

        //getTimeStephWidth from preCICE

        if (t+dt_ > settings_.endTime){
            dt_ = settings_.endTime-t;
        }
        computeTemperature(tTmp);
        t += dt_;
        if (t - lastOutputTime > settings_.outputFileEveryDt - 1e-4) {
            cout << "current time: " << t << " dt: " << dt_ << " pressure solver iterations: " << endl;
            outputWriterParaview_->writeFile(t);
            outputWriterText_->writeFile(t);
            lastOutputTime = t;
        }
    }

    if (std::fabs(t - lastOutputTime) > 1e-4) {
        outputWriterParaview_->writeFile(t);
        lastOutputTime = t;
    }

}


//TODO 
/*
add 
(i) resetting of the time step, 
(ii) counting of data points at the couplig interface,
(iii) writing fluxes at the interface to a separate data structure to be sent by preCICE,
(iv) reading of temperature data at the grd points of the coupling interface that have 
 been received from the flow solver via preCICE similar to these changes in the flow solver.
*/



void Computation_solid::computeTemperature(Array2D tTmp) {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if (geometry_.get()->isFluid(i, j)) {
                tTmp(i,j) = discretization_.get()->t(i, j);
            }
        }
    }
    
    temperatureSolver_->solve();
}

/* former computeTemperature
void Computation_solid::computeTemperature() {
    Array2D tTmp(discretization_.get()->t().size());
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if (geometry_.get()->isFluid(i, j)) { //TODO nur dann oder überall?
                tTmp(i,j) = discretization_.get()->t(i, j) +
                                                 dt_ * (1 / settings_.re * 1 / settings_.prandtl * (
                                                         discretization_.get()->computeD2TDx2(i, j)
                                                         +
                                                         discretization_.get()->computeD2TDy2(i, j)
                                                 )
                                                        - discretization_.get()->computeDuTDx(i, j)
                                                        - discretization_.get()->computeDvTDy(i, j)
                                                 );
            }
        }
    }
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if (geometry_.get()->isFluid(i, j)) {
                discretization_.get()->t(i, j) = tTmp(i,j);
            }
        }
    }
}
*/

void Computation_solid::applyInitialConditions() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if (!geometry_.get()->isFluid(i, j)) {
                if (!geometry_.get()->isFluid(i + 1, j) && !geometry_.get()->isFluid(i - 1, j) &&
                    !geometry_.get()->isFluid(i, j + 1) && !geometry_.get()->isFluid(i, j - 1)) {
                    discretization_.get()->u(i, j) = std::nan("");
                    discretization_.get()->v(i, j) = std::nan("");
                    discretization_.get()->p(i, j) = std::nan("");
                    discretization_.get()->t(i, j) = std::nan("");
                } else {
                    discretization_.get()->u(i, j) = 0;
                    discretization_.get()->v(i, j) = 0;
                    discretization_.get()->p(i, j) = 0; //TODO what to do?
                    discretization_.get()->t(i, j) = tInit;
                }
            } else {
                discretization_.get()->u(i, j) = 0;
                discretization_.get()->v(i, j) = 0;
                discretization_.get()->p(i, j) = 0; //TODO what to do?
                discretization_.get()->t(i, j) = tInit;
            }
        }
    }
}


void Computation_solid::copyOldValues() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++){
        for (int i  = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++){
            discretization_.get()->tOld(i,j) = discretization_.get()->t(i,j);
        }
    }
}