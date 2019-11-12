//
// Created by Fabian on 26.10.2019.
//
#include <iostream>
#include "Computation/Computation.h"
#include "Computation/ComputationParallel.h"
#include <chrono>
#include <mpi.h>

using namespace std::chrono;

int main(int argc, char *argv[]) {
    // Todo parse number of processes
    if (argc == 1 ){
        std::cout << "Please pass settings filename" << std::endl;
        return EXIT_FAILURE;
    }

    //cout << "Running with " << number_of_processes <<  " processes" << endl;

    //Determine parition
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    cout << "rank is " << rank << endl;
    MPI_Finalize();
    return 0;
>>>>>>> cdb7b7649a52749a7826af0c88dc39e686b7fdb6

    ComputationParallel computation;
    computation.initialize(argc, argv);
    auto start = high_resolution_clock::now();
    computation.runSimulation();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop -start);
    cout << "Time in milliseconds: " << duration.count() << endl;

    MPI_Finalize();
    return EXIT_SUCCESS;


}

