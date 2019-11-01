//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_PRESSURESOLVER_H
#define CODE_NUMSIM_PRESSURESOLVER_H

#include <memory>
#include "../StaggeredGrid/Discretization.h"

class PressureSolver {
public:
    PressureSolver(std::shared_ptr <Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    virtual void solve() = 0;

protected:
    void setBoundaryValues();

    std::shared_ptr <Discretization> discretization_;

    double epsilon_;

    int maximumNumberOfIterations_;
};


#endif //CODE_NUMSIM_PRESSURESOLVER_H
