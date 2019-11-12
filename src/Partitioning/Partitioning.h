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

private:
    int nCellsLocal;
    const int rank;
    const int rankLeft;
    const int rankRight;
    const int rankTop;
    const int rankBottom;
};

#endif //NUMSIM2019_PARTITIONING_H
