//
// Created by Julia Pelzer on 26.10.2019.
//

#include "SOR.h"
#include <cmath>
#include <iostream>

SOR::SOR(std::shared_ptr<Discretization> discretization, std::shared_ptr<Communication> communication,
         Partitioning partitioning, double epsilon, int maximumNumberOfIterations) :
        PressureSolver(discretization, communication, epsilon, maximumNumberOfIterations, partitioning),
        omega(omega) {

}

void SOR::solve() {
    int iter = 0;
    double eps = 1;
    double epsAll = 1;
    setBoundaryValues();
    while (iter < maximumNumberOfIterations_ && eps > pow(epsilon_, 2)) {
        // first iteration over "red" values
        int decision = (partitioning.nodeOffset()[0] + partitioning.nodeOffset()[1]) % 2;

        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin() + (j % 2 + decision) % 2;
                 i <= discretization_.get()->pIEnd(); i += 2) {
                discretization_.get()->p(i, j) = (1 - omega) *
                discretization_.get()->p(i, j) +
                omega *
                pow(discretization_.get()->dx(), 2) *
                pow(discretization_.get()->dy(), 2)
                /
                (2 * (pow(discretization_.get()->dx(), 2) +
                      pow(discretization_.get()->dy(), 2)))
                *
                (
                        (discretization_.get()->p(i - 1, j) +
                         discretization_.get()->p(i + 1, j)) /
                        pow(discretization_.get()->dx(), 2)
                        +
                        (discretization_.get()->p(i, j - 1) +
                         discretization_.get()->p(i, j + 1)) /
                        pow(discretization_.get()->dy(), 2) -
                        discretization_.get()->rhs(i, j));

            }
        }
        // communication of new "red" values
        communication_.get()->communicate(discretization_.get()->p(), "p");
        setBoundaryValues();

        // second iteration over "black" values
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin() + 1 - (j % 2 + decision) % 2;
                 i <= discretization_.get()->pIEnd(); i += 2) {
                discretization_.get()->p(i, j) = (1 - omega) *
                discretization_.get()->p(i, j) +
                omega *
                pow(discretization_.get()->dx(), 2) *
                pow(discretization_.get()->dy(), 2)
                /
                (2 * (pow(discretization_.get()->dx(), 2) +
                      pow(discretization_.get()->dy(), 2)))
                *
                (
                        (discretization_.get()->p(i - 1, j) +
                         discretization_.get()->p(i + 1, j)) /
                        pow(discretization_.get()->dx(), 2)
                        +
                        (discretization_.get()->p(i, j - 1) +
                         discretization_.get()->p(i, j + 1)) /
                        pow(discretization_.get()->dy(), 2) -
                        discretization_.get()->rhs(i, j));
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

        // collection of global residual
        MPI_Allreduce(&eps, &epsAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        eps = epsAll / (partitioning.nCellsGlobal()[0] * partitioning.nCellsGlobal()[1]);
        iter++;
        setBoundaryValues();
    }
    std::cout << "pressure solver iterations: " << iter << " eps :" << eps << " epslion² " << pow(epsilon_,2) <<std::endl;
    setBoundaryValues();
}
