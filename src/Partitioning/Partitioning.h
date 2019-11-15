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
    int getSize();
    std::array<int, 2> getNCells();


private:
    std::array<int, 2> nCellsLocal;
    int rank;
    int rankLeft;
    int rankRight;
    int rankTop;
    int rankBottom;
    int size;


};

#endif //NUMSIM2019_PARTITIONING_H
