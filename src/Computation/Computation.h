//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_COMPUTATION_H
#define CODE_NUMSIM_COMPUTATION_H

#include "StaggeredGrid/StaggeredGrid.h"
#include "Settings.h"
#include <memory>
#include "StaggeredGrid/Discretization.h"
#include "PressureSolver/PressureSolver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"

class Computation
{
public:
    virtual void initialize(int argc, char *argv[]);

    virtual void runSimulation();

private:
    virtual void computeTimeStepWidth();

    virtual void applyBoundaryValues();

    virtual void PreliminaryVelocities();

    virtual void computeRightHandSide();

    virtual void computePressure();

    virtual void computeVelocities();

    // Attributes
    Settings settings_;

    std::shared_ptr<Discretization> discretization_;

    std::unique_ptr<PressureSolver> pressureSolver_;

    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;

    std::unique_ptr<OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_;

    double dt_;
};

#endif //CODE_NUMSIM_COMPUTATION_H
