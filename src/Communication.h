//
// Created by fritz on 11/12/19.
//

#ifndef NUMSIM2019_COMMUNICATION_H
#define NUMSIM2019_COMMUNICATION_H

#include "Partitioning/Partitioning.h"
#include <memory>
#include <mpi.h>

class Communication
{
public:
    Communication(std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Discretization> ptr);

    void communicate_fg();

private:
    std::shared_ptr<Partitioning> partitioning_;
    std::shared_ptr<Discretization> discretization_;


    void exchangeValues(std::vector<double> &values, int receiverRank, std::vector<MPI_Request> &requests);
};

#endif //NUMSIM2019_COMMUNICATION_H
