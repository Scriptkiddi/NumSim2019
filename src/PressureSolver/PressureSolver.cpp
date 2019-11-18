//
// Created by Julia Pelzer on 26.10.2019.
//

#include "PressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization,std::shared_ptr<Communication> communication, double epsilon,
                               int maximumNumberOfIterations, Partitioning partitioning) :
        discretization_(discretization),
        communication_(communication),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations),partitioning(partitioning) {

}

void PressureSolver::setBoundaryValues() {
    //rechter und linker Rand ohne ecken
    if (partitioning.getRankOfRightNeighbour() == -1) {
        int i_high = discretization_.get()->pIEnd() + 1;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            discretization_.get()->p(i_high, j) = discretization_.get()->p(i_high - 1, j);
        }
    }
    if (partitioning.getRankOfLeftNeighbour() == -1) {
        int i_low = discretization_.get()->pIBegin() - 1;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            discretization_.get()->p(i_low, j) = discretization_.get()->p(i_low + 1, j);
        }
    }
    //unterer Rand mit ecken
    if (partitioning.getRankOfTopNeighbour() == -1) {
        int j_high = discretization_.get()->pJEnd() + 1;
        for (int i = discretization_.get()->pIBegin() - 1; i <= discretization_.get()->pIEnd() + 1; i++) {
            discretization_.get()->p(i, j_high) = discretization_.get()->p(i, j_high - 1);
        }
    }
    if (partitioning.getRankOfBottomNeighbour() == -1) {
        int j_low = discretization_.get()->pJBegin() - 1;
        for (int i = discretization_.get()->pIBegin() - 1; i <= discretization_.get()->pIEnd() + 1; i++) {
            discretization_.get()->p(i, j_low) = discretization_.get()->p(i, j_low + 1);
        }
    }
}