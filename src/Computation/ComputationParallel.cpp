//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include "ComputationParallel.h"

ComputationParallel::ComputationParallel(std::string settingsFilename) : Computation(settingsFilename),
                                                                  partitioning_(
                                                                          Partitioning(settings_.nCells)) {

}

void ComputationParallel::initialize(int argc, char **argv) {
    partitioning_.getRank()
    array<int, 2> nCellsBoundary = {settings_.nCells[0] + 2, settings_.nCells[1] + 2}; // Mit Ghost cells

    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] /
                    (nCellsBoundary[0]-2);
    meshWidth_[1] = settings_.physicalSize[1] / (nCellsBoundary[1]-2);

    //initialize discretization
    if (!settings_.useDonorCell) {
        CentralDifferences grid(nCellsBoundary, meshWidth_);
        discretization_ = make_shared<CentralDifferences>(grid);
    } else {
        DonorCell grid(nCellsBoundary, meshWidth_, settings_.alpha);
        discretization_ = make_shared<DonorCell>(grid);
    }

    //initialize explicit pressureSolver
    if (settings_.pressureSolver == "SOR") {
        SOR pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
        pressureSolver_ = make_unique<SOR>(pSolver);
    } else {
        //GaussSeidel pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
        //pressureSolver_ = make_unique<PressureSolver>(pSolver);
        std::cout << "Please select SOR-solver" << std::endl;
    }
    //initialize outputWriters
    OutputWriterText outText(discretization_);
    outputWriterText_ = make_unique<OutputWriterText>(outText);

    OutputWriterParaview outPara(discretization_);
    outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);


}


void ComputationParallel::computeTimeStepWidth() {
    double uMaximum = 0;
    for (int j = discretization_.get()->uJBegin()-1; j <= discretization_.get()->uJEnd()+1; j++) {
        for (int i = discretization_.get()->uIBegin()-1; i <= discretization_.get()->uIEnd()+1; i++) {
            if (uMaximum < abs(discretization_.get()->u(i, j))) {
                uMaximum = abs(discretization_.get()->u(i, j));
            }
        }
    }
    double vMaximum = 0;
    for (int j = discretization_.get()->vJBegin()-1; j <= discretization_.get()->vJEnd()+1; j++) {
        for (int i = discretization_.get()->vIBegin()-1; i <= discretization_.get()->vIEnd()+1; i++) {
            if (vMaximum < abs(discretization_.get()->v(i, j))) {
                vMaximum = abs(discretization_.get()->v(i, j));
            }
        }
    }
    double condition_diffusion = pow(discretization_.get()->dx(), 2) * pow(discretization_.get()->dy(), 2) /
                                 (pow(discretization_.get()->dx(), 2) + pow(discretization_.get()->dy(), 2)) *
                                 settings_.re /
                                 2;
    double condition_convection1 = discretization_.get()->dx() / uMaximum;
    double condition_convection2 = discretization_.get()->dy() / vMaximum;

    dt_ = min(condition_convection1, condition_convection2);
    dt_ = min(condition_diffusion, dt_);
    dt_ = min(settings_.maximumDt, dt_) * settings_.tau;

    MPI_Allreduce(&dt_, &dtAll_, partitioning_.getSize(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_ = dtAll_;

}
