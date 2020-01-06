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
#include <precice/SolverInterface.hpp>

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
    applyInitialConditions();
    double lastOutputTime = 0;

    //INIT Precice
    precice::SolverInterface solverInterface(this->settings_.participantName, 0, 1);
    //load precice config
    solverInterface.configure(settings_.preciceConfigFile);

    //configure ids for read and write data attrs
    int meshID = solverInterface.getMeshID(settings_.meshName);
    int writeDataID = solverInterface.getDataID(settings_.writeDataName, meshID);
    int readDataID = solverInterface.getDataID(settings_.readDataName, meshID);

    // Get vertex size
    int vertexSize = 0;
    int dim = solverInterface.getDimensions();
    for (int j = discretization_.get()->tJBegin()-1; j <= discretization_.get()->tJEnd()+1; j++) {
        for (int i = discretization_.get()->tIBegin()-1; i <= discretization_.get()->tIEnd()+1; i++) {
            if(geometry_.get()->get_temperature(i,j).first == "TPD" || geometry_.get()->get_temperature(i,j).first == "TPN"){
                vertexSize++;
            }
        }
    }

    // init strucutres for data transfer
    double coords[vertexSize*dim];
    double* temperature = new double[vertexSize];
    double* heatFlow = new double[vertexSize];
    int* vertexIDs = new int[vertexSize];

    // init coords TODO are i and j coords?
    int k=0;
    for (int j = discretization_.get()->tJBegin()-1; j <= discretization_.get()->tJEnd()+1; j++) {
        for (int i = discretization_.get()->tIBegin()-1; i <= discretization_.get()->tIEnd()+1; i++) {
            if(geometry_.get()->get_temperature(i,j).first == "TPD" || geometry_.get()->get_temperature(i,j).first == "TPN"){
                coords[k] = i;
                coords[+1] = j;
                k += 2;
            }
        }
    }

    // Set Verticies
    solverInterface.setMeshVertices(meshID,vertexSize,coords,vertexIDs);

    // Create Checkpoints
    const std::string& coric = precice::constants::actionReadIterationCheckpoint();
    const std::string& cowic = precice::constants::actionWriteIterationCheckpoint();
    static const std::string& cowid = precice::constants::actionWriteInitialData();

    // finally call init
    solverInterface.initialize();
    if ( solverInterface.isActionRequired(cowid)){
        computeTemperature();
        int k=0;
        for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
            for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
                if(geometry_.get()->get_temperature(i,j).first == "TDP" || geometry_.get()->get_temperature(i,j).first == "TNP"){
                    temperature[k] = geometry_.get()->get_temperature(i,j).second[0];
                    k++;
                }
            }
        }
        solverInterface.writeBlockScalarData(writeDataID, vertexSize, vertexIDs, temperature);
        solverInterface.fulfilledAction(cowid);

    }

    solverInterface.initializeData();
    for (int timeStepNumber = 0;
        std::abs(t - settings_.endTime) > 1e-10 && settings_.endTime - t > 0; timeStepNumber++) {

        // Save old state and acknowledge checkpoint
        if(solverInterface.isActionRequired(cowic)){
            saveOldState(); // save checkpoint
            solverInterface.fulfilledAction(cowic);
        }

        if (t+dt_ > settings_.endTime){
            dt_ = settings_.endTime-t;
        }

        // Since we are the solid solver we wait till the fluid has completed the temp calculation and has transfered
        // his data over to us
        solverInterface.readBlockScalarData(writeDataID, vertexSize, vertexIDs, temperature);
        cout << "asdf" << endl;

        // write recieved data for boundary values to our geometry mesh
        // TODO not sure if index by k for the data transfer works, it should but better think about it again
        k=0;
        for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
            for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
                if(geometry_.get()->get_temperature(i,j).first == "TDP" || geometry_.get()->get_temperature(i,j).first == "TNP"){
                    geometry_.get()->get_temperature(i,j).second[0] = temperature[k];
                    k++;
                }
            }
        }
        // got all the required data calculate time step width
        computeTimeStepWidth();
        computeTemperature();
        solverInterface.writeBlockScalarData(readDataID, vertexSize, vertexIDs, heatFlow);
        k=0;
        for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
            for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
                if(geometry_.get()->get_temperature(i,j).first == "TDP" || geometry_.get()->get_temperature(i,j).first == "TNP"){
                    heatFlow[k] = geometry_.get()->get_temperature(i,j).second[0];
                    k++;
                }
            }
        }
        dt_ = solverInterface.advance(dt_);
        if(solverInterface.isActionRequired(coric)) { // timestep not converged
            reloadOldState(); // set variables back to checkpoint
            solverInterface.fulfilledAction(coric);
        } else { // timestep converged
            // e.g. update variables, increment time
            t += dt_;
        }

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


void Computation_solid::computeTemperature() {

    temperatureSolver_->solve();
}

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

void Computation_solid::reloadOldState() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            discretization_.get()->t(i,j) = discretization_.get()->tOld(i,j);
        }
    }
}

void Computation_solid::saveOldState() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++){
        for (int i  = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++){
            discretization_.get()->tOld(i,j) = discretization_.get()->t(i,j);
        }
    }
}


void Computation_solid::computeTimeStepWidth() {
    // Min time requirements for temperature
    double condition_temp = settings_.re * settings_.prandtl * 0.5 *
                            pow((1 / pow(discretization_.get()->dx(), 2) + 1 / pow(discretization_.get()->dy(), 2)),
                                -1);

    dt_ = min(settings_.maximumDt, condition_temp) * settings_.tau;
}
