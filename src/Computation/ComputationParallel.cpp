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
