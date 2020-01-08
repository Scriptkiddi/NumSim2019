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
    GaussSeidel tSolver(discretization_, geometry_, settings_.epsilon, settings_.maximumNumberOfIterations, dt_,
                        settings_.heatDiffusivity);
    temperatureSolver_ = make_unique<GaussSeidel>(tSolver); //TODO not working 

    //initialize outputWriters

    OutputWriterText outText(discretization_, settings);
    outputWriterText_ = make_unique<OutputWriterText>(outText);

    OutputWriterParaview outPara(discretization_, settings);
    outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);
}

void Computation_solid::runSimulation() {
    double t = 0;
    applyInitialConditions();
    double lastOutputTime = 0;

#ifdef WITHPRECICE
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
    for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
        for (int i = discretization_.get()->tIBegin() - 1; i <= discretization_.get()->tIEnd() + 1; i++) {
            if (geometry_.get()->get_temperature(i, j).first == "TPD" ||
                geometry_.get()->get_temperature(i, j).first == "TPN") {
                vertexSize++;
            }
        }
    }

    // init strucutres for data transfer
    double coords[vertexSize * dim];
    double temperature[vertexSize];
    double heatFlow[vertexSize];
    int *vertexIDs = new int[vertexSize];

    // init coords
    int k = 0;
    for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
        for (int i = discretization_.get()->tIBegin() - 1; i <= discretization_.get()->tIEnd() + 1; i++) {
            if (geometry_.get()->get_temperature(i, j).first == "TPD" ||
                geometry_.get()->get_temperature(i, j).first == "TPN") {
                if (j > 0 && geometry_.get()->isFluid(i, j - 1)) {
                    coords[k] = settings_.origin[0] + discretization_.get()->dx() * (i - 1) +
                                discretization_.get()->dx() / 2;
                    coords[k + 1] = settings_.origin[1] + discretization_.get()->dy() * (j - 1);
                } else if (geometry_.get()->isFluid(i, j + 1)) {
                    coords[k] = settings_.origin[0] + discretization_.get()->dx() * (i - 1) +
                                discretization_.get()->dx() / 2;
                    coords[k + 1] =
                            settings_.origin[1] + discretization_.get()->dy() * (j - 1) + discretization_.get()->dy();
                } else if (i > 0 && geometry_.get()->isFluid(i - 1, j)) {
                    coords[k] = settings_.origin[0] + discretization_.get()->dx() * (i - 1);
                    coords[k + 1] =
                            settings_.origin[1] + discretization_.get()->dy() * (j - 1) +
                            discretization_.get()->dy() / 2;
                } else if (geometry_.get()->isFluid(i + 1, j)) {
                    coords[k] = settings_.origin[0] + discretization_.get()->dx() * (i) + discretization_.get()->dx();
                    coords[k + 1] =
                            settings_.origin[1] + discretization_.get()->dy() * (j - 1) +
                            discretization_.get()->dy() / 2;
                }
                std::cout << "Compute coupling point " << k / 2 << " at (x,y) = (" << coords[k] << ", " << coords[k + 1]
                          << ")" << std::endl;
                k += 2;
            }
        }
    }

    // Set Verticies
    solverInterface.setMeshVertices(meshID, vertexSize, coords, vertexIDs);

    // Create Checkpoints
    static const std::string &cowid = precice::constants::actionWriteInitialData();
    const std::string &coric = precice::constants::actionReadIterationCheckpoint();
    const std::string &cowic = precice::constants::actionWriteIterationCheckpoint();

    // finally call init
    double precice_dt = solverInterface.initialize();

    if (solverInterface.isActionRequired(cowid)) {
        // TODO is this needed
        //computeTemperature();
        int k = 0;
        for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
            for (int i = discretization_.get()->tIBegin() - 1; i <= discretization_.get()->tIEnd() + 1; i++) {
                if (geometry_.get()->get_temperature(i, j).first == "TPN") {
                    if (i >= discretization_.get()->tIBegin() && geometry_.get()->isFluid(i - 1, j)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i - 1, j));
                    } else if (i <= discretization_.get()->tIEnd() && geometry_.get()->isFluid(i + 1, j)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i + 1, j));
                    } else if (j >= discretization_.get()->tJBegin() && geometry_.get()->isFluid(i, j - 1)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i, j - 1));
                    } else if (j <= discretization_.get()->tJEnd() && geometry_.get()->isFluid(i, j + 1)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i, j + 1));
                    }
                    k++;
                } else if (geometry_.get()->get_temperature(i, j).first == "TPD") {
                    if (i >= discretization_.get()->tIBegin() && geometry_.get()->isFluid(i - 1, j)) {
                        heatFlow[k] = 1 / (discretization_.get()->dx() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i - 1, j) - discretization_.get()->t(i, j));
                    } else if (i <= discretization_.get()->tIEnd() && geometry_.get()->isFluid(i + 1, j)) {
                        heatFlow[k] = 1 / (discretization_.get()->dx() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i + 1, j) - discretization_.get()->t(i, j));
                    } else if (j >= discretization_.get()->tJBegin() && geometry_.get()->isFluid(i, j - 1)) {
                        heatFlow[k] = 1 / (discretization_.get()->dy() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i, j - 1) - discretization_.get()->t(i, j));
                    } else if (j <= discretization_.get()->tJEnd() && geometry_.get()->isFluid(i, j + 1)) {
                        heatFlow[k] = 1 / (discretization_.get()->dy() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i, j + 1) - discretization_.get()->t(i, j));
                    }
                    k++;
                }
            }
        }

        solverInterface.writeBlockScalarData(writeDataID, vertexSize, vertexIDs, heatFlow);
        solverInterface.fulfilledAction(cowid);
    }

    solverInterface.initializeData();
        //Do preCICE calls here
#endif

    int timeStepNumber = 0;
#ifdef WITHPRECICE
    while (solverInterface.isCouplingOngoing()) {
#else
    dt_ = 0.001;
    while (true) {
#endif
#ifdef WITHPRECICE
        // Save old state and acknowledge checkpoint
        if (solverInterface.isActionRequired(cowic)) {
            //saveOldState(); // save checkpoint
            solverInterface.fulfilledAction(cowic);
        }
#endif

        // Calculate fluid time step
#ifdef WITHPRECICE
        dt_ =  precice_dt;
        // Read Temperature
        solverInterface.readBlockScalarData(readDataID, vertexSize, vertexIDs, temperature);
        k = 0;
        for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
            for (int i = discretization_.get()->tIBegin() - 1; i <= discretization_.get()->tIEnd() + 1; i++) {
                if (geometry_.get()->get_temperature(i, j).first == "TPD") {
                    std::pair<std::string, std::vector<double>> value = geometry_.get()->get_temperature(i, j);
                    value.second[0] = temperature[k];
                    geometry_.get()->set_temperature(i, j, value);
                    k++;
                } else if (geometry_.get()->get_temperature(i, j).first == "TPN") {
                    std::pair<std::string, std::vector<double>> value = geometry_.get()->get_temperature(i, j);
                    value.second[0] = -temperature[k];
                    geometry_.get()->set_temperature(i, j, value);
                    k++;
                }
            }
        }
#endif
        computeTemperature();
        computeRightHandSide();

        // Coupling

        // write to preCice Buffers
#ifdef WITHPRECICE
        int k = 0;
        for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
            for (int i = discretization_.get()->tIBegin() - 1; i <= discretization_.get()->tIEnd() + 1; i++) {
                if (geometry_.get()->get_temperature(i, j).first == "TPN") {
                    if (i >= discretization_.get()->tIBegin() && geometry_.get()->isFluid(i - 1, j)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i - 1, j));
                    } else if (i <= discretization_.get()->tIEnd() && geometry_.get()->isFluid(i + 1, j)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i + 1, j));
                    } else if (j >= discretization_.get()->tJBegin() && geometry_.get()->isFluid(i, j - 1)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i, j - 1));
                    } else if (j <= discretization_.get()->tJEnd() && geometry_.get()->isFluid(i, j + 1)) {
                        heatFlow[k] = 0.5 * (discretization_.get()->t(i, j) + discretization_.get()->t(i, j + 1));
                    }
                    k++;
                } else if (geometry_.get()->get_temperature(i, j).first == "TPD") {
                    if (i >= discretization_.get()->tIBegin() && geometry_.get()->isFluid(i - 1, j)) {
                        heatFlow[k] = 1 / (discretization_.get()->dx() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i - 1, j) - discretization_.get()->t(i, j));
                    } else if (i <= discretization_.get()->tIEnd() && geometry_.get()->isFluid(i + 1, j)) {
                        heatFlow[k] = 1 / (discretization_.get()->dx() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i + 1, j) - discretization_.get()->t(i, j));
                    } else if (j >= discretization_.get()->tJBegin() && geometry_.get()->isFluid(i, j - 1)) {
                        heatFlow[k] = 1 / (discretization_.get()->dy() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i, j - 1) - discretization_.get()->t(i, j));
                    } else if (j <= discretization_.get()->tJEnd() && geometry_.get()->isFluid(i, j + 1)) {
                        heatFlow[k] = 1 / (discretization_.get()->dy() * settings_.heatDiffusivity) *
                                         (discretization_.get()->t(i, j + 1) - discretization_.get()->t(i, j));
                    }
                    k++;
                }
            }
        }
        solverInterface.writeBlockScalarData(writeDataID, vertexSize, vertexIDs, temperature);

        // Advance fluid solver
        precice_dt = solverInterface.advance(dt_);
#endif

        // reset if required
#ifdef WITHPRECICE
        if (solverInterface.isActionRequired(coric)) { // timestep not converged
            //reloadOldState(); // set variables back to checkpoint
            solverInterface.fulfilledAction(coric);
        } else { // timestep converged
#else
        if (true) {
#endif
            std::cout << "Solid: Advancing int time!" << std::endl;
            // e.g. update variables, increment time
            t += dt_;
            if (!(std::abs(t - settings_.endTime) > 1e-10 && settings_.endTime - t > 0)) {
                break;
            }
            timeStepNumber++;
            outputWriterParaview_->writeFile(t, "solid");
            if (t - lastOutputTime > settings_.outputFileEveryDt - 1e-4) {
                cout << "current time: " << t << " dt: " << dt_ << " pressure solver iterations: " << endl;
                outputWriterParaview_->writeFile(t, "solid");
                lastOutputTime = t;
            }
        }
    }


    if (std::fabs(t - lastOutputTime) > 1e-4) {
        outputWriterParaview_->writeFile(t, "solid");
    }
#ifdef WITHPRECICE
    solverInterface.finalize();
#endif


}


//TODO 
/*
add 
(ii) counting of data points at the couplig interface,
(iii) writing fluxes at the interface to a separate data structure to be sent by preCICE,
(iv) reading of temperature data at the grd points of the coupling interface that have 
 been received from the flow solver via preCICE similar to these changes in the flow solver.
*/


void Computation_solid::computeTemperature() {

    temperatureSolver_->solve(dt_);
}

void Computation_solid::computeRightHandSide() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if (geometry_.get()->isFluid(i, j)) {
                discretization_->rhs(i, j) = discretization_.get()->t(i, j) / dt_;
            } else {
                discretization_->rhs(i, j) = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
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
                    discretization_.get()->p(i, j) = 0;
                    discretization_.get()->rhs(i, j) = 0;
                    discretization_.get()->t(i, j) = settings_.tInit;
                    discretization_.get()->tOld(i, j) = settings_.tInit;
                }
            } else {
                discretization_.get()->u(i, j) = 0;
                discretization_.get()->v(i, j) = 0;
                discretization_.get()->p(i, j) = 0;
                discretization_.get()->rhs(i, j) = 0;
                discretization_.get()->t(i, j) = settings_.tInit;
                discretization_.get()->tOld(i, j) = settings_.tInit;
            }
        }
    }
}

void Computation_solid::reloadOldState() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            discretization_.get()->t(i, j) = discretization_.get()->tOld(i, j);
        }
    }
}

void Computation_solid::saveOldState() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            discretization_.get()->tOld(i, j) = discretization_.get()->t(i, j);
        }
    }
}


