//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_COMPUTATION_H
#define CODE_NUMSIM_COMPUTATION_H

#include "Settings.h"
#include "memory"
#include "StaggeredGrid/Discretization.h"
#include "PressureSolver/PressureSolver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"


class Computation {
public:
    void intitialize(int argc, char *argv[]);

    void runSimulation();

private:
    void computeTimeStepWidth();

    void applyBoundaryValues();

    void PreliminaryVelocities();

    void computeRightHandSide();

    void computePressure();

    void computeVelocities();

    // Attributes
    Settings settings_;

    std::shared_ptr <Discretization> discretization_;

    std::unique_ptr <PressureSolver> pressureSolver_;

    std::unique_ptr <OutputWriterParaview> outputWriterParaview_;

    std::unique_ptr <OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_;

    double dt_;
};


#endif //CODE_NUMSIM_COMPUTATION_H
