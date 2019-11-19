//
// Created by fritz on 11/12/19.
//

#include <array>
#ifndef NUMSIM2019_PARTITIONING_H
#define NUMSIM2019_PARTITIONING_H

class Partitioning
{
public:
    explicit Partitioning(std::array<int, 2> nCells, std::array<double,2>);

    int getRank();
    int getRankOfLeftNeighbour();
    int getRankOfRightNeighbour();
    int getRankOfBottomNeighbour();
    int getRankOfTopNeighbour();
    int getSize();
    const std::array<int,2> nCellsGlobal();
    const bool ownPartitionContainsRightBoundary();
    const bool ownPartitionContainsTopBoundary();
    std::array<int, 2> getNCells();
    std::array<double, 2> getMeshWidth();
    int ownRankNo();

    std::array<int, 2> nodeOffset();

private:
    std::array<int, 2> nCellsLocal;
    int rank;
    int rankLeft;
    int rankRight;
    int rankTop;
    int rankBottom;
    int size;
    const std::array<int,2> nCellsGlobal_;
    std::array<int,2> nodeOffset_;
    std::array<double,2> physicalSize_;


};

#endif //NUMSIM2019_PARTITIONING_H
