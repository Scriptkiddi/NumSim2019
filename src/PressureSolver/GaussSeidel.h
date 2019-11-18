//
// Created by Julia Pelzer on 26.10.2019.
//
#include <memory>
#include "PressureSolver.h"

#ifndef CODE_NUMSIM_GAUSSSEIDEL_H
#define CODE_NUMSIM_GAUSSSEIDEL_H


class GaussSeidel : public PressureSolver
{
public:
    GaussSeidel(std::shared_ptr<Discretization> sharedPtr, std::shared_ptr<Communication> sharedPtr1,
                Partitioning partitioning, double epsilon, int maximumNumberOfIterations);

    void solve() override;

private:
    double omega;

    void setBoundaryValues();
};


#endif //CODE_NUMSIM_GAUSSSEIDEL_H
