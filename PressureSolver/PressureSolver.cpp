//
// Created by Julia Pelzer on 26.10.2019.
//

#include "PressureSolver.h"

PressureSolver(std::shared_ptr <Discretization>
discretization,
double epsilon,
int maximumNumberOfIterations)

virtual void solve() = 0

protected:

void setBoundaryValues()

std::shared_ptr <Discretization> discretization_

double epsilon_

int maximumNumberOfIterations_
