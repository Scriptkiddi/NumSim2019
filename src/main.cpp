//
// Created by Fabian on 26.10.2019.
//
#include <iostream>
#include "Computation/Computation.h"
#include "Computation/ComputationParallel.h"
#include <chrono>
#include <mpi.h>
#include <Util/InputParser.h>

using namespace std::chrono;

int main(int argc, char *argv[]) {
    // Todo parse number of processes
    if (argc == 1 ){
        std::cout << "Please pass settings filename" << std::endl;
        return EXIT_FAILURE;
    }

    int number_of_processes = 0;
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-n")){
        number_of_processes = stoi(input.getCmdOption("-n"));
    }else{
        cout << "Please supply the number of Processes to use" << endl;
        return EXIT_FAILURE;
    }
    if(number_of_processes <= 0){
        cout << "Please supply a number of Processes greater or equal to 1" << endl;
        return EXIT_FAILURE;
    }

    cout << "Running with " << number_of_processes <<  " processes" << endl;

    //Determine partition
    MPI_Init(NULL, NULL);


    ComputationParallel computation;
    computation.initialize(argc, argv);
    auto start = high_resolution_clock::now();
    computation.runSimulation();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop -start);
    cout << "Time in miliseconds: " << duration.count() << endl;

    return EXIT_SUCCESS;


}

