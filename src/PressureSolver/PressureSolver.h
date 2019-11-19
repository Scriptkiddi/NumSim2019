//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_PRESSURESOLVER_H
#define CODE_NUMSIM_PRESSURESOLVER_H

#include <memory>
#include "../StaggeredGrid/Discretization.h"
#include "Communication.h"
class PressureSolver {
public:

    PressureSolver(std::shared_ptr<Discretization> discretization, std::shared_ptr<Communication> communication, Partitioning partitioning,
                   double epsilon, int maximumNumberOfIterations);

    virtual void solve() = 0;

protected:
    void setBoundaryValues();

    std::shared_ptr <Discretization> discretization_;

    std::shared_ptr<Communication> communication_;
    Partitioning partitioning;

    double epsilon_;

    int maximumNumberOfIterations_;
};


#endif //CODE_NUMSIM_PRESSURESOLVER_H
