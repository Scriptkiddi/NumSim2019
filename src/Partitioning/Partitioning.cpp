//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include "Partitioning.h"
#include <array>
#include <vector>
#include <iostream>
#include <cmath>


// Partitions are numbered from left to right starting with 0. After a line is finished you start for the next line again from the left side

Partitioning::Partitioning(std::array<int, 2> nCells) {
    nCellsLocal = nCells;

    //Determine partition
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int x = nCells[0];
    int y = nCells[1];
    bool count_x;
    std::vector<std::array<int, 2>> possibilities;
    for (int i = 1; i < size; i++) {
        if (size % i == 0) {
            possibilities.push_back({i, size / i});
        }
    }
    int min = size;
    int index = 0;
    for (int i = 0; i < possibilities.size(); i++) {
        if (abs(possibilities[i][0] - possibilities[i][1]) < min) {
            index = i;
        }
    }
    int numberX = possibilities[index][0];
    int numberY = possibilities[index][1];

    // Determine if special case (Left border, bottom, ...)
    rankRight = rank +1;
    rankLeft = rank - 1;
    rankBottom = rank-numberX;
    rankTop = rank+numberX;
    if(rank < numberX ){ //Bottom Border
        rankBottom = -1;
    }
    if(rank % numberX == numberX-1){ // Right Border -1 since numberX is a size value
        rankRight = -1;
    }
    if (rank >= size-numberX) { //Top Border
        rankTop = -1;
    }
    if(rank % numberX==0){ // Left Border
        rankLeft = -1;
    }

    int localCellsX = 0;
    int localCellsY = 0;

    if(rankRight != -1){
        localCellsX = floor(nCells[0]/numberX);
    }else{
        localCellsX = nCells[0] - floor(nCells[0]/numberX)*(numberX-1);
    }
    if(rankTop != -1){
        localCellsY = floor(nCells[1]/numberY);
    }else{
        localCellsY = nCells[1] - floor(nCells[1]/numberY)*(numberY-1);
    }
    std::cout << rank << "|" << localCellsX << "x" << localCellsY << std::endl;
    std::cout << rank << "|" << rankBottom << "|" << rankTop <<  "|" << rankLeft <<  "|" << rankRight << std::endl;




}


    int Partitioning::getRank(){
        return rank;
    }

    int Partitioning::getRankOfLeftNeighbour(){
        return rankLeft;
    }

    int Partitioning::getRankOfRightNeighbour(){
        return rankRight;
    }

    int Partitioning::getRankOfBottomNeighbour(){
        return rankBottom;
    }

    int Partitioning::getRankOfTopNeighbour(){
        return rankTop;
    }

int Partitioning::getSize() {
    return size;
}
