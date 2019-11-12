//
// Created by fritz on 11/12/19.
//

#ifndef NUMSIM2019_COMPUTATIONPARALLEL_H
#define NUMSIM2019_COMPUTATIONPARALLEL_H

#include "Computation.h"
#include "../Partitioning/Partitioning.h"
#include <string>

class ComputationParallel : public Computation
{
public:
    ComputationParallel(string settingsFilename);

    virtual void initialize(int argc, char *argv[]);

    //virtual void runSimulation();

private:
    double dtAll_;
    void computeTimeStepWidth();

    //virtual void applyBoundaryValues();

    //virtual void PreliminaryVelocities();

    //virtual void computeRightHandSide();

    //virtual void computePressure();

    //virtual void computeVelocities();

    //Attributes

    Partitioning partitioning_;

};

#endif //NUMSIM2019_COMPUTATIONPARALLEL_H
