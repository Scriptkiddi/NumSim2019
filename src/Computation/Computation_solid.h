//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_COMPUTATION_SOLID_H
#define CODE_NUMSIM_COMPUTATION_SOLID_H

#include "StaggeredGrid/StaggeredGrid.h"
#include "Settings.h"
#include <memory>
#include "StaggeredGrid/Discretization.h"
#include "TemperatureSolver/TemperatureSolver.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "Geometry.h"

class Computation_solid {
public:
    void initialize(int argc, char *argv[]);

    void runSimulation();

    void computerightHandSide();

private:

    void computeTemperature();

    void applyInitialConditions();
    void computeRightHandSide();

    void applyBoundaryValuesTemperature();

    void computeTimeStepWidth();

    void saveOldState();

    void reloadOldState();

    // Attributes
    Settings settings_;

    std::shared_ptr<Geometry> geometry_;

    std::shared_ptr <Discretization> discretization_;
    
    std::unique_ptr <TemperatureSolver> temperatureSolver_;

    std::unique_ptr <OutputWriterParaview> outputWriterParaview_;

    std::unique_ptr <OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_;

    double dt_;
    
    double tInit;

    void calcHeatFlow(double pDouble[-1]);
};


#endif //CODE_NUMSIM_COMPUTATION_SOLID_H
