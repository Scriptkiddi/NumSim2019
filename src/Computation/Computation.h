//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_COMPUTATION_H
#define CODE_NUMSIM_COMPUTATION_H

#include "StaggeredGrid/StaggeredGrid.h"
#include "Settings.h"
#include <memory>
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"
#include "Geometry.h"

class Computation {
public:
    void initialize(int argc, char *argv[]);

    void runSimulation();

private:
    void computeTimeStepWidth();

    void applyBoundaryValuesF();

    void computeMacroscopicQuantities(int t);

    void applyInitialConditions();

    void applyLatticeVelocities();

    void applyWeights();

    void computeFtempFeq(int t);

    void computeF(int t);

    double computeBoundaryValue(int i, int j, int k, std::string boundary);

    // Attributes
    Settings settings_;

    std::shared_ptr<Geometry> geometry_;

    std::shared_ptr <StaggeredGrid> staggeredGrid_;

    std::unique_ptr <OutputWriterParaview> outputWriterParaview_;

    std::unique_ptr <OutputWriterText> outputWriterText_;

    std::array<double, 2> meshWidth_;

    double dt_;

    double t_c;
    double t_h;
    
    double fInit;

    double tau_;
};


#endif //CODE_NUMSIM_COMPUTATION_H
