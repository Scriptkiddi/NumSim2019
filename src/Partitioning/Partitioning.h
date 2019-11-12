//
// Created by fritz on 11/12/19.
//

#include <array>
#ifndef NUMSIM2019_PARTITIONING_H
#define NUMSIM2019_PARTITIONING_H

class Partitioning
{
public:
    explicit Partitioning(std::array<int, 2> nCells);


    int getRank();
    int getRankOfLeftNeighbour();
    int getRankOfRightNeighbour();
    int getRankOfBottomNeighbour();
    int getRankOfTopNeighbour();
    std::array<int, 2> nCellsLocal;
};

#endif //NUMSIM2019_PARTITIONING_H
