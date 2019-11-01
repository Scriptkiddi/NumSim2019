//
// Created by Julia Pelzer on 26.10.2019.
//

#include "PressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon,
                               int maximumNumberOfIterations) :
        discretization_(discretization),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations) {

}

void PressureSolver::setBoundaryValues() {
    //todo: Zellen in den Ecken werden noch ignoriert, evtl. nach Bedarf ergänzen
    //unterer Rand
    int j = discretization_.get()->pJBegin();
    for (int i = discretization_.get()->pIBegin() + 1; i <= discretization_.get()->pIEnd() - 1; i++) {
        discretization_.get()->p(i, j) = discretization_.get()->p(i, j + 1);
    }
    //rechter und linker Rand
    int i_low = discretization_.get()->pIBegin();
    int i_high = discretization_.get()->pIEnd();
    for (j = discretization_.get()->pJBegin() + 1; j <= discretization_.get()->pJEnd() - 1; j++) {
        discretization_.get()->p(i_low, j) = discretization_.get()->p(i_low + 1, j);
        discretization_.get()->p(i_high, j) = discretization_.get()->p(i_high - 1, j);
    }
    // oberer Rand
    j = discretization_.get()->pJEnd();
    for (int i = discretization_.get()->pIBegin() + 1; i <= discretization_.get()->pIEnd() - 1; i++) {
        discretization_.get()->p(i, j) = discretization_.get()->p(i, j - 1);
    }
    discretization_.get()->p(discretization_.get()->pIBegin(), discretization_.get()->pJBegin()) =
            discretization_.get()->p(discretization_.get()->pIBegin()+1, discretization_.get()->pJBegin()+1);

    discretization_.get()->p(discretization_.get()->pIBegin(), discretization_.get()->pJEnd()) =
            discretization_.get()->p(discretization_.get()->pIBegin()+1, discretization_.get()->pJEnd()-1);

    discretization_.get()->p(discretization_.get()->pIEnd(), discretization_.get()->pJBegin()) =
            discretization_.get()->p(discretization_.get()->pIEnd()-1, discretization_.get()->pJBegin()+1);

    discretization_.get()->p(discretization_.get()->pIEnd(), discretization_.get()->pJEnd()) =
            discretization_.get()->p(discretization_.get()->pIEnd()-1, discretization_.get()->pJEnd()-1);
}