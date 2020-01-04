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
#include "Geometry.h"

class Computation {
public:
    void initialize(int argc, char *argv[]);

    void runSimulation();

private:
    void computeTimeStepWidth();

    void applyBoundaryValuesVelocities();

    void applyBoundaryValuesTemperature();

    void PreliminaryVelocities();

    void computeRightHandSide();

    void computePressure();

    void computeVelocities();

    void computeTemperature();

    void applyInitialConditions();

    void copyOldValues();

    // Attributes
    Settings settings_;

    std::shared_ptr<Geometry> geometry_;

    std::shared_ptr <Discretization> discretization_;

    std::unique_ptr <PressureSolver> pressureSolver_;

    std::unique_ptr <OutputWriterParaview> outputWriterParaview_;

    std::unique_ptr <OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_;

    double dt_;

    double t_c;
    double t_h;
    
    double uInit;
    double vInit;
    double pInit; 
    double tInit;
};


#endif //CODE_NUMSIM_COMPUTATION_H
