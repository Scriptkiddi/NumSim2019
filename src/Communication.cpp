#include "Communication.h"

Communication::Communication(Partitioning partitioning)
{
    partitioning_ = std::make_unique<Partitioning>(partitioning);
}

void Communication::communicate_f()
{
    
}