//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_GAUSSSEIDEL_H
#define CODE_NUMSIM_GAUSSSEIDEL_H

#include "TemperatureSolver.h"

class GaussSeidel : public TemperatureSolver {
public:
    GaussSeidel(std::shared_ptr<Discretization> discretization, std::shared_ptr<Geometry> geometry, double epsilon, int maximumNumberOfIterations, double dt, double heatDiffusivity);
    void solve(Array2D tTmp) override;
};


#endif //CODE_NUMSIM_GAUSSSEIDEL_H