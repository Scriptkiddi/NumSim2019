//
// Created by fritz on 11/12/19.
//

#ifndef NUMSIM2019_COMPUTATIONPARALLEL_H
#define NUMSIM2019_COMPUTATIONPARALLEL_H

#include "Computation.h"
#include "Communication.h"
#include "../Partitioning/Partitioning.h"
#include <string>
#include <output_writer/output_writer_text_parallel.h>
#include "output_writer/output_writer_paraview_parallel.h"

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




    std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaview_;

    std::unique_ptr<OutputWriterTextParallel> outputWriterText_;

    //Attributes

    Partitioning partitioning_;

    std::shared_ptr<Communication> communication_;
};

#endif //NUMSIM2019_COMPUTATIONPARALLEL_H
