//
// Created by Julia Pelzer on 26.10.2019.
//
#include <StaggeredGrid/Discretization.h>
#include <cmath>
#include <iostream>
#include "GaussSeidel.h"

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization, std::shared_ptr<Communication> communication,
                         Partitioning partitioning, double epsilon, int maximumNumberOfIterations) :

        PressureSolver(discretization, communication, epsilon, maximumNumberOfIterations), partitioning(partitioning) {
}


void GaussSeidel::setBoundaryValues() {

    //rechter und linker Rand ohne ecken
    if (partitioning.getRankOfRightNeighbour() == -1) {
        int i_high = discretization_.get()->pIEnd() + 1;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            discretization_.get()->p(i_high, j) = -discretization_.get()->p(i_high - 1, j);
        }
    }
    if (partitioning.getRankOfLeftNeighbour() == -1) {
        int i_low = discretization_.get()->pIBegin() - 1;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            discretization_.get()->p(i_low, j) = -discretization_.get()->p(i_low + 1, j);
        }
    }
    //unterer Rand mit ecken
    if (partitioning.getRankOfTopNeighbour() == -1) {
        int j_high = discretization_.get()->pJEnd() + 1;
        for (int i = discretization_.get()->pIBegin() - 1; i <= discretization_.get()->pIEnd() + 1; i++) {
            discretization_.get()->p(i, j_high) = -discretization_.get()->p(i, j_high - 1);
        }
    }
    if (partitioning.getRankOfBottomNeighbour() == -1) {
        int j_low = discretization_.get()->pJBegin() - 1;
        for (int i = discretization_.get()->pIBegin() - 1; i <= discretization_.get()->pIEnd()+1; i++) {
            discretization_.get()->p(i, j_low) = -discretization_.get()->p(i, j_low + 1);
        }
    }
}

void GaussSeidel::solve() {
    int iter = 0;
    double eps = 1;
    double epsAll = 1;
    double factor = 0.5
                    * std::pow(discretization_.get()->dx(), 2) * std::pow(discretization_.get()->dy(), 2)
                    / (std::pow(discretization_.get()->dx(), 2) + std::pow(discretization_.get()->dy(), 2));
    while (iter < maximumNumberOfIterations_ && epsAll > pow(epsilon_, 2)) {
        setBoundaryValues();

        // first iteration over "red" values
        int decision = (partitioning.nodeOffset()[1] + partitioning.nodeOffset()[1]) % 2;
        // (j % 2 + decision) % 2
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin() + j % 2; i <= discretization_.get()->pIEnd(); i += 2) {
                discretization_.get()->p(i, j) = factor
                                                 * ((discretization_.get()->p(i - 1, j) +
                                                     discretization_.get()->p(i + 1, j))
                                                    / std::pow(discretization_.get()->dx(), 2)
                                                    + (discretization_.get()->p(i, j - 1) +
                                                       discretization_.get()->p(i, j + 1))
                                                      / std::pow(discretization_.get()->dy(), 2)
                                                    - discretization_.get()->rhs(i, j));
            }
        }
        //// communication of new "red" values
        ////std::cout << "Communicating p1" << std::endl;
        communication_.get()->communicate(discretization_.get()->p(), "p");
        setBoundaryValues();

        // second iteration over "black" values
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin() + 1 - j % 2; i <= discretization_.get()->pIEnd(); i += 2) {
                discretization_.get()->p(i, j) = factor
                                                 * ((discretization_.get()->p(i - 1, j) +
                                                     discretization_.get()->p(i + 1, j))
                                                    / std::pow(discretization_.get()->dx(), 2)
                                                    + (discretization_.get()->p(i, j - 1) +
                                                       discretization_.get()->p(i, j + 1))
                                                      / std::pow(discretization_.get()->dy(), 2)
                                                    - discretization_.get()->rhs(i, j));
            }
        }
        // communication over new "black" values
        //std::cout << "Communicating p2" << std::endl;
        communication_.get()->communicate(discretization_.get()->p(), "p");
        setBoundaryValues();

        // local residual
        eps = 0;
        epsAll = 0;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
                eps = eps + pow(
                        discretization_->rhs(i, j) -
                        (discretization_.get()->p(i - 1, j) - 2 * discretization_.get()->p(i, j) +
                         discretization_.get()->p(i + 1, j)) / pow(discretization_.get()->dx(), 2) -
                        (discretization_.get()->p(i, j - 1) - 2 * discretization_.get()->p(i, j) +
                         discretization_.get()->p(i, j + 1)) / pow(discretization_.get()->dy(), 2), 2);
            }
        }
        eps = eps / (discretization_.get()->nCells()[0] * discretization_.get()->nCells()[1]);

        // collection of global residual
        MPI_Allreduce(&eps, &epsAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        eps = epsAll;
        iter++;
    }
    std::cout << "pressure solver iterations: " << iter << " eps :" << eps << " epsilon² " << pow(epsilon_, 2)
              << std::endl;
    setBoundaryValues();
}



