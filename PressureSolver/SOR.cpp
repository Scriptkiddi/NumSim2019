//
// Created by Julia Pelzer on 26.10.2019.
//

#include "SOR.h"
#include <cmath>

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega) :
        PressureSolver(discretization, epsilon, maximumNumberOfIterations),
        omega(omega) {

}

void SOR::solve() {
    int iter = 1;
    double eps = 1;
    while (iter <= maximumNumberOfIterations_ && eps > epsilon_) {
        iter++;
        setBoundaryValues();
        for (int j = discretization_.get()->pJBegin() + 1; j <= discretization_.get()->pJEnd() - 1; j++) {
            for (int i = discretization_.get()->pIBegin() + 1; i <= discretization_.get()->pIEnd() - 1; i++) {
                discretization_.get()->p(i, j) = (1 - omega) * discretization_.get()->p(i, j) +
                                                 omega *
                                                 pow(discretization_.get()->dx(), 2) *
                                                 pow(discretization_.get()->dy(), 2) /
                                                 2 * (pow(discretization_.get()->dx(), 2) +
                                                      pow(discretization_.get()->dy(), 2)) *
                                                 ((discretization_.get()->p(i - 1, j) +
                                                   discretization_.get()->p(i + 1, j)) /
                                                  pow(discretization_.get()->dx(), 1) +
                                                  (discretization_.get()->p(i, j - 1) +
                                                   discretization_.get()->p(i, j + 1)) /
                                                  pow(discretization_.get()->dy(), 2) -
                                                  discretization_.get()->rhs(i, j));
            }
        }
        //residuum
        eps = 0;
        for (int j = discretization_.get()->pJBegin() + 1; j <= discretization_.get()->pJEnd() - 1; j++) {
            for (int i = discretization_.get()->pIBegin() + 1; i <= discretization_.get()->pIEnd() - 1; i++) {
                eps = eps + pow((discretization_.get()->p(i - 1, j) -
                                 2 * discretization_.get()->p(i, j) +
                                 discretization_.get()->p(i + 1, j)) /
                                pow(discretization_.get()->dx(), 2) +
                                (discretization_.get()->p(i, j - 1) -
                                 2 * discretization_.get()->p(i, j) +
                                 discretization_.get()->p(i, j + 1)) /
                                pow(discretization_.get()->dy(), 2) -
                                discretization_.get()->rhs(i, j),
                                2);
            }
        }
        eps = sqrt(eps);
    }
    setBoundaryValues();
}
