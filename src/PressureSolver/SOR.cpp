//
// Created by Julia Pelzer on 26.10.2019.
//

#include "SOR.h"
#include <cmath>
#include <iostream>

SOR::SOR(std::shared_ptr<Discretization> discretization, std::shared_ptr<Geometry> geometry,double epsilon, int maximumNumberOfIterations, double omega) :
        PressureSolver(discretization, geometry, epsilon, maximumNumberOfIterations),
        omega(omega) {

}

void SOR::solve() {
    int iter = 0;
    double eps = 1;
    while (iter <= maximumNumberOfIterations_ && eps > pow(epsilon_,2)) {
        setBoundaryValues();
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
                if(geometry_.get()->isFluid(i,j)){
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
        }
        //residuum
        eps = 0;
        for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
            for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
                if(geometry_.get()->isFluid(i,j)){
                    eps = eps + pow(
                            discretization_->rhs(i, j)
                            - (discretization_.get()->p(i - 1, j) - 2 * discretization_.get()->p(i, j) +
                            discretization_.get()->p(i + 1, j))
                            / pow(discretization_.get()->dx(), 2)
                            - (discretization_.get()->p(i, j - 1) - 2 * discretization_.get()->p(i, j) +
                            discretization_.get()->p(i, j + 1))
                            / pow(discretization_.get()->dy(), 2), 2);
                }
            }
        }
        eps = eps/geometry_.get()->nCellsFluid();
        //eps = eps / ((discretization_.get()->nCells()[0] - 2) * (discretization_.get()->nCells()[1] - 2));
        iter++;
    }
    std::cout << "pressure solver iterations: " << iter << " eps :" << eps << " epslion² " << pow(epsilon_,2) <<std::endl;
    setBoundaryValues();
}
