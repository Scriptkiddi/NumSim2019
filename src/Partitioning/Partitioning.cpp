//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include "Partitioning.h"
#include <array>
#include <vector>
#include <iostream>


// Partitions are numbered from left to right starting with 0. After a line is finished you start for the next line again from the left side

Partitioning::Partitioning(std::array<int, 2> nCells) {
    nCellsLocal = nCells;

    //Determine partition
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    int x = nCells[0];
    int y = nCells[1];
    bool count_x;
    std::vector<std::array<int, 2>> possibilities;
    for (int i = 1; i < commSize; i++) {
        if (commSize % i == 0) {
            possibilities.push_back({i, commSize / i});
        }
    }
    int min = commSize;
    int index = 0;
    for (int i = 0; i < possibilities.size(); i++) {
        if (abs(possibilities[i][0] - possibilities[i][1]) < min) {
            index = i;
        }
    }
    int local_x = possibilities[index][0];
    int local_y = possibilities[index][1];

    std::cout << local_x << "x" << local_y << std::endl;


    //cout << "rank is " << rank << "/" << commSize << endl;

}

