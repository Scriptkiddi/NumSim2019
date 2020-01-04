//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_TEMPERATURESOLVER_H
#define CODE_NUMSIM_TEMPERATURESOLVER_H

#include <memory>
#include "../StaggeredGrid/Discretization.h"
#include "../Geometry.h"

class TemperatureSolver {
public:
    TemperatureSolver(std::shared_ptr <Discretization> discretization, std::shared_ptr<Geometry> geometry,double epsilon, int maximumNumberOfIterations, double dt, double heatDiffusivity);

    virtual void solve(Array2D tTmp) = 0;

protected:
    void applyBoundaryValuesTemperature();
    //void setBoundaryValues();

    std::shared_ptr <Discretization> discretization_;

    std::shared_ptr <Geometry> geometry_;

    double epsilon_;

    int maximumNumberOfIterations_;

    double dt_;

    double alpha_;
};


#endif //CODE_NUMSIM_TEMPERATURESOLVER_H
