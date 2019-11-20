//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include "Partitioning.h"
#include <array>
#include <vector>
#include <iostream>
#include <cmath>
#include <climits>


// Partitions are numbered from left to right starting with 0. After a line is finished you start for the next line again from the left side

Partitioning::Partitioning(std::array<int, 2> nCells, std::array<double, 2> physicalSize): 
                            nCellsGlobal_(nCells),
                            physicalSizeGlobal_(physicalSize){

    //Determine partition
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout << "test" << std::endl;

    double dx = physicalSize[0] / nCells[0];
    double dy = physicalSize[1] / nCells[1];

    int x = nCells[0];
    int y = nCells[1];
    bool count_x;
    std::vector<std::array<int, 2>> possibilities;
    if(size==1){
        possibilities.push_back({1,1});
    }
    for (int i = 1; i <= size/2; i++) {
        if (size % i == 0) {

            std::cout << i << " x " << size/i << std::endl;
            possibilities.push_back({i, size / i});
        }
    }
    int min = INT_MAX;
    int index = 0;
    for (int i = 0; i < possibilities.size(); i++) {
        if (abs(possibilities[i][0] - possibilities[i][1]) < min) {
            min = abs(possibilities[i][0]-possibilities[i][1]);
            index = i;
        }
    }
    int numberX;
    int numberY;
    if ((possibilities[index][0] > possibilities[index][1] &&nCells[0] > nCells[1]) || (possibilities[index][0] < possibilities[index][1] &&nCells[0] < nCells[1])){
        numberX = possibilities[index][0];
        numberY = possibilities[index][1];
    }else{
        numberX = possibilities[index][1];
        numberY = possibilities[index][0];
    }
    std::cout << "Split for processors " << numberX << " x " << numberY << std::endl;
    
    // Determine if special case (Left border, bottom, ...)
    rankRight = rank + 1;
    rankLeft = rank - 1;
    rankBottom = rank - numberX;
    rankTop = rank + numberX;
    // Rank begins with 0
    if (rank < numberX) { //Bottom Border
        rankBottom = -1;
    }
    if (rank % numberX == numberX - 1) { // Right Border -1 since numberX is a size value
        rankRight = -1;
    }
    if (rank >= size - numberX) { //Top Border
        rankTop = -1;
    }
    if (rank % numberX == 0) { // Left Border
        rankLeft = -1;
    }
    if (rankRight != -1) {
        nCellsLocal[0] = floor(nCells[0] / numberX);
    } else {
        nCellsLocal[0] = nCells[0] - floor(nCells[0] / numberX) * (numberX - 1);
    }
    nodeOffset_[0] = floor(nCells[0] / numberX) * (rank % numberX);
    if (rankTop != -1) {
        nCellsLocal[1] = floor(nCells[1] / numberY);
    } else {
        nCellsLocal[1] = nCells[1] - floor(nCells[1] / numberY) * (numberY - 1);
    }
    nodeOffset_[1] = floor(rank / numberX)  * floor(nCells[1] / numberY) ;
    std::cout << rank << " | nCells " << nCellsLocal[0] << " x " << nCellsLocal[1] << std::endl;

    physicalSize_[0] = dx * nCellsLocal[0]; // physicalSize[0]/numberX;
    physicalSize_[1] = dy * nCellsLocal[1]; // physicalSize[1]/numberY;
    std::cout << "test2" << std::endl;
}

int Partitioning::getRank() {
    return rank;
}

int Partitioning::getRankOfLeftNeighbour() {
    return rankLeft;
}

int Partitioning::getRankOfRightNeighbour() {
    return rankRight;
}

int Partitioning::getRankOfBottomNeighbour() {
    return rankBottom;
}

int Partitioning::getRankOfTopNeighbour() {
    return rankTop;
}

std::array<int, 2> Partitioning::getNCells() {
    return nCellsLocal;
}


int Partitioning::getSize() {
    return size;
}

const std::array<int, 2> Partitioning::nCellsGlobal() {
    return std::array<int,2> {nCellsGlobal_[0], nCellsGlobal_[1]};
    //return nCellsGlobal_;
}

const bool Partitioning::ownPartitionContainsRightBoundary() {
    return getRankOfRightNeighbour() == -1;
}

const bool Partitioning::ownPartitionContainsTopBoundary() {
    return getRankOfTopNeighbour() == -1;
}

int Partitioning::ownRankNo() {
    return rank;
}

std::array<int, 2> Partitioning::nodeOffset() {
    return nodeOffset_;
}

std::array<double, 2> Partitioning::getMeshWidth() {
    return {physicalSizeGlobal_[0]/nCellsGlobal_[0], physicalSizeGlobal_[1]/nCellsGlobal_[1]};
//    return {physicalSize_[0]/(getNCells()[0]), physicalSize_[1]/(getNCells()[1])};
}
