//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_SOR_H
#define CODE_NUMSIM_SOR_H

#include "PressureSolver.h"

class SOR : public PressureSolver
{
public:
    SOR(std::shared_ptr<Discretization> discretization, std::shared_ptr<Communication> communication, double epsilon, int maximumNumberOfIterations, double omega);

    void solve() override;

private:
    double omega;
};

#endif //CODE_NUMSIM_SOR_H
