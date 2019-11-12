//
// Created by fritz on 11/12/19.
//

#ifndef NUMSIM2019_PARTITIONING_H
#define NUMSIM2019_PARTITIONING_H

class Partitioning
{
public:
    int getRank();
    int getRankOfLeftNeighbour();
    int getRankOfRightNeighbour();
    int getRankOfBottomNeighbour();
    int getRankOfTopNeighbour();
    int nCellsLocal;
};

#endif //NUMSIM2019_PARTITIONING_H
