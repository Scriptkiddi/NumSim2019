//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include <StaggeredGrid/CentralDifferences.h>
#include <StaggeredGrid/DonorCell.h>
#include <Communication.h>
#include "ComputationParallel.h"

ComputationParallel::ComputationParallel(std::string settingsFilename) : Computation(settingsFilename),
                                                                         partitioning_(
                                                                                 Partitioning(settings_.nCells)) {

}

void ComputationParallel::runSimulation() {
    double t = 0;
    applyBoundaryValues();

    while (t < settings_.endTime) {
        computeTimeStepWidth();
        applyBoundaryValues();
        PreliminaryVelocities();
        communication_.get()->communicate(discretization_.get()->f(), "f");
        communication_.get()->communicate(discretization_.get()->g(), "g");
        computeRightHandSide();
        //computePressure();
        computeVelocities();
        communication_.get()->communicate(discretization_.get()->f(), "u");
        communication_.get()->communicate(discretization_.get()->g(), "v");
        t += dt_;
        //outputWriterParaview_.get()->writeFile(t);
        //outputWriterText_.get()->writeFile(t);
        //cout << partitioning_.getRank()  << "|current time: " << t << " dt: " << dt_ << " pressure solver iterations: " << endl;
    }


}

void ComputationParallel::initialize(int argc, char **argv) {
    array<int, 2> nCellsBoundary = {partitioning_.getNCells()[0],
                                    partitioning_.getNCells()[1]};

    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] /
                    (nCellsBoundary[0]);
    meshWidth_[1] = settings_.physicalSize[1] / (nCellsBoundary[1]);


    //init grid

    int nX = partitioning_.getNCells()[0];
    int nY = partitioning_.getNCells()[1];
    FieldVariable u = FieldVariable({nX + 2, nY + 2}, {0 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    FieldVariable v = FieldVariable({nX + 2, nY + 2}, {-0.5 * meshWidth_[0], -0 * meshWidth_[1]}, meshWidth_);
    if (partitioning_.getRankOfLeftNeighbour() == -1) {
        FieldVariable u = FieldVariable({nX + 1, nY + 2}, {0 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    }
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        FieldVariable v = FieldVariable({nX + 1, nY + 2}, {-0.5 * meshWidth_[0], 0 * meshWidth_[1]}, meshWidth_);
    }
    FieldVariable p = FieldVariable({nX + 2, nY + 2}, {-0.5 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    FieldVariable rhs = FieldVariable({nX, nY}, {-0.5 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    FieldVariable f = FieldVariable({nX + 2, nY + 2}, {-0.5 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    FieldVariable g = FieldVariable({nX + 2, nY + 2}, {-0.5 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    if (partitioning_.getRankOfLeftNeighbour() == -1) {
        FieldVariable f = FieldVariable({nX + 1, nY + 2}, {-0.5 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    }
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        FieldVariable g = FieldVariable({nX + 1, nY + 2}, {-0.5 * meshWidth_[0], -0.5 * meshWidth_[1]}, meshWidth_);
    }


    //initialize discretization
    if (!settings_.useDonorCell) {
        CentralDifferences grid(nCellsBoundary, meshWidth_, make_shared<Partitioning>(partitioning_), u, v, p, f, g, rhs);
        discretization_ = make_shared<CentralDifferences>(grid);
    } else {
        DonorCell grid(nCellsBoundary, meshWidth_, settings_.alpha, make_shared<Partitioning>(partitioning_), u, v, p, f, g, rhs);
        discretization_ = make_shared<DonorCell>(grid);
    }
    communication_ = make_shared<Communication>(make_shared<Partitioning>(partitioning_), discretization_);

    //initialize explicit pressureSolver
    //if (settings_.pressureSolver == "SOR") {
    //    SOR pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    //    pressureSolver_ = make_unique<SOR>(pSolver);
    //} else {
    //    //GaussSeidel pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    //    //pressureSolver_ = make_unique<PressureSolver>(pSolver);
    //    std::cout << "Please select SOR-solver" << std::endl;
    //}
    //initialize outputWriters
    //OutputWriterText outText(discretization_);
    //outputWriterText_ = make_unique<OutputWriterText>(outText);

    //OutputWriterParaview outPara(discretization_);
    //outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);


}


void ComputationParallel::computeTimeStepWidth() {
    double uMaximum = 0;
    for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
        for (int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd() + 1; i++) {
            if (uMaximum < abs(discretization_.get()->u(i, j))) {
                uMaximum = abs(discretization_.get()->u(i, j));
            }
        }
    }
    double vMaximum = 0;
    for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
        for (int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++) {
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

    MPI_Allreduce(&dt_, &dtAll_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_ = dtAll_;

}

void ComputationParallel::applyBoundaryValues() {
    // U
    //Linker Rand
    if (partitioning_.getRankOfLeftNeighbour() == -1) {
        int i_low = discretization_.get()->uIBegin() - 1;
        for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
            discretization_.get()->u(i_low, j) = settings_.dirichletBcLeft[0];
            discretization_.get()->f(i_low, j) = discretization_.get()->u(i_low, j);
        }
    }

    // Rechter Rand
    if (partitioning_.getRankOfRightNeighbour() == -1) {
        int i_high = discretization_.get()->uIEnd() + 1;
        for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
            discretization_.get()->u(i_high, j) = settings_.dirichletBcRight[0];
            discretization_.get()->f(i_high, j) = discretization_.get()->u(i_high, j);
        }
    }

    //unterer Rand
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        int j_low = discretization_.get()->uJBegin() - 1;
        for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
            discretization_.get()->u(i, j_low) =
                    2 * settings_.dirichletBcBottom[0] - discretization_.get()->u(i, j_low + 1);
            discretization_.get()->f(i, j_low) = discretization_.get()->u(i, j_low);
        }
    }
    // oberer Rand
    if (partitioning_.getRankOfTopNeighbour() == -1) {
        int j_high = discretization_.get()->uJEnd() + 1;
        for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
            discretization_.get()->u(i, j_high) =
                    2 * settings_.dirichletBcTop[0] - discretization_.get()->u(i, j_high - 1);
            discretization_.get()->f(i, j_high) = discretization_.get()->u(i, j_high);
        }
    }

    // v
    //unterer Rand
    // oberer Rand
    if (partitioning_.getRankOfTopNeighbour() == -1) {
        int j_high = discretization_.get()->vJEnd() + 1;
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            discretization_.get()->v(i, j_high) = settings_.dirichletBcTop[1];
            discretization_.get()->g(i, j_high) = discretization_.get()->v(i, j_high);
        }
    }
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        int j_low = discretization_.get()->vJBegin() - 1;
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            discretization_.get()->v(i, j_low) = settings_.dirichletBcBottom[1];
            discretization_.get()->g(i, j_low) = discretization_.get()->v(i, j_low);
        }
    }

    //rechter und linker Rand
    if (partitioning_.getRankOfRightNeighbour() == -1) {
        int i_high = discretization_.get()->vIEnd() + 1;
        for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
            discretization_.get()->v(i_high, j) =
                    2 * settings_.dirichletBcRight[1] - discretization_.get()->v(i_high - 1, j);
            discretization_.get()->g(i_high, j) = discretization_.get()->v(i_high, j);
        }
    }
    if (partitioning_.getRankOfRightNeighbour() == -1) {
        int i_low = discretization_.get()->vIBegin() - 1;
        for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
            discretization_.get()->v(i_low, j) =
                    2 * settings_.dirichletBcLeft[1] - discretization_.get()->v(i_low + 1, j);
            discretization_.get()->g(i_low, j) = discretization_.get()->v(i_low, j);

        }
    }


}
