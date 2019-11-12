//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include "ComputationParallel.h"

ComputationParallel::ComputationParallel(std::string settingsFilename) : Computation(settingsFilename),
                                                                  partitioning_(
                                                                          Partitioning(settings_.nCells)) {

}

void ComputationParallel::initialize(int argc, char **argv) {


}

