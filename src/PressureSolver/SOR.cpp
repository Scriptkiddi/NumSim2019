//
// Created by Julia Pelzer on 26.10.2019.
//

#include "SOR.h"
#include <cmath>
#include <iostream>

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations),
        omega(omega) {

}

void SOR::solve() {
    int iter = 0;
    double eps = 1;
    while (iter <= maximumNumberOfIterations_ && eps > epsilon_) {
        //setBoundaryValues();
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
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
        //residuum
        eps = 0;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
                eps = eps + pow(
                        discretization_->rhs(i, j)
                        - (discretization_.get()->p(i - 1, j) - 2 * discretization_.get()->p(i, j) +
                           discretization_.get()->p(i + 1, j))
                          / pow(discretization_.get()->meshWidth()[0], 2)
                        - (discretization_.get()->p(i, j - 1) - 2 * discretization_.get()->p(i, j) +
                           discretization_.get()->p(i, j + 1))
                          / pow(discretization_.get()->meshWidth()[1], 2), 2);
            }
        }
        eps = sqrt(eps);
        iter++;
    }
    std::cout << "pressure solver iterations: " << iter << std::endl;
    setBoundaryValues();
}
