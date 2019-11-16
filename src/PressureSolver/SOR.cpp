//
// Created by Julia Pelzer on 26.10.2019.
//

#include "SOR.h"
#include <cmath>
#include <iostream>

SOR::SOR(std::shared_ptr<Discretization> discretization, std::shared_ptr<Communication> communication,
         double epsilon, int maximumNumberOfIterations, double omega) :
        PressureSolver(discretization, communication, epsilon, maximumNumberOfIterations),
        omega(omega) {
}

void SOR::solve() {
    int iter = 0;
    double eps = 1;
    double epsAll = 1;
    while (iter < maximumNumberOfIterations_ && epsAll > pow(epsilon_, 2)) {
        setBoundaryValues();
        // first iteration over "red" values
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin() + j % 2; i <= discretization_.get()->pIEnd(); i += 2) {
                discretization_.get()->p(i, j) = (1 - omega) *
                                                 discretization_.get()->p(i, j) +
                                                 omega *
                                                 pow(discretization_.get()->dx(), 2) *
                                                 pow(discretization_.get()->dy(), 2) /
                                                 (2 * (pow(discretization_.get()->dx(), 2) +
                                                       pow(discretization_.get()->dy(), 2))) *
                                                 ((discretization_.get()->p(i - 1, j) +
                                                   discretization_.get()->p(i + 1, j)) /
                                                  pow(discretization_.get()->dx(), 2) +
                                                  (discretization_.get()->p(i, j - 1) +
                                                   discretization_.get()->p(i, j + 1)) /
                                                  pow(discretization_.get()->dy(), 2) -
                                                  discretization_.get()->rhs(i, j));
            }
        }
        // communication of new "red" values
        //std::cout << "Communicating p1" << std::endl;
        communication_.get()->communicate(discretization_.get()->p(), "p");

        // second iteration over "black" values
        for (int j = discretization_.get()->pJBegin() + 1; j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin() + 1 - j % 2; i <= discretization_.get()->pIEnd(); i += 2) {
                discretization_.get()->p(i, j) = (1 - omega) *
                                                 discretization_.get()->p(i, j) +
                                                 omega *
                                                 pow(discretization_.get()->dx(), 2) *
                                                 pow(discretization_.get()->dy(), 2) /
                                                 (2 * (pow(discretization_.get()->dx(), 2) +
                                                       pow(discretization_.get()->dy(), 2))) *
                                                 ((discretization_.get()->p(i - 1, j) +
                                                   discretization_.get()->p(i + 1, j)) /
                                                  pow(discretization_.get()->dx(), 2) +
                                                  (discretization_.get()->p(i, j - 1) +
                                                   discretization_.get()->p(i, j + 1)) /
                                                  pow(discretization_.get()->dy(), 2) -
                                                  discretization_.get()->rhs(i, j));
            }
        }
        // communication over new "black" values
        //std::cout << "Communicating p2" << std::endl;
        communication_.get()->communicate(discretization_.get()->p(), "p");

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
        //epsAll = eps;
        MPI_Allreduce(&eps, &epsAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        iter++;
    }
    std::cout << "pressure solver iterations: " << iter << " eps :" << eps << " epsilon² " << pow(epsilon_, 2)
              << std::endl;
    setBoundaryValues();
}
