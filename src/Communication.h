//
// Created by fritz on 11/12/19.
//

#ifndef NUMSIM2019_COMMUNICATION_H
#define NUMSIM2019_COMMUNICATION_H

#include "Partitioning/Partitioning.h"
#include <memory>

class Communication
{
public:
    Communication(Partitioning partitioning);

    void communicate_f();

private:
    std::unique_ptr<Partitioning> partitioning_;
};

#endif //NUMSIM2019_COMMUNICATION_H
