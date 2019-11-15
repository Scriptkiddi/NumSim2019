//
// Created by fritz on 11/12/19.
//

#ifndef NUMSIM2019_COMPUTATIONPARALLEL_H
#define NUMSIM2019_COMPUTATIONPARALLEL_H

#include "Computation.h"
#include "Communication.h"
#include "../Partitioning/Partitioning.h"
#include <string>

class ComputationParallel : public Computation
{
public:
    ComputationParallel(string settingsFilename);

    void initialize(int argc, char *argv[]);

    void runSimulation();

private:
    double dtAll_;
    void computeTimeStepWidth();

    void applyBoundaryValues();

    //virtual void PreliminaryVelocities();

    //virtual void computeRightHandSide();

    //virtual void computePressure();

    //virtual void computeVelocities();

    //Attributes

    Partitioning partitioning_;

    std::shared_ptr<Communication> communication_;
};

#endif //NUMSIM2019_COMPUTATIONPARALLEL_H
