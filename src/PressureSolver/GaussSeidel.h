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
    GaussSeidel(std::shared_ptr<Discretization> discretization, std::shared_ptr<Communication> communication, double epsilon,
                int maximumNumberOfIterations);

    void solve() override;

private:
    double omega;
};


#endif //CODE_NUMSIM_GAUSSSEIDEL_H
