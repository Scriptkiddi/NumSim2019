//
// Created by Fabian on 26.10.2019.
//
#include <iostream>
#include "Computation/Computation.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[]) {
    if (argc == 1 ){
        std::cout << "Please pass settings filename" << std::endl;
        return EXIT_FAILURE;
    }
    
    Computation computation;

    computation.initialize(argc, argv);

    auto start = high_resolution_clock::now();
    computation.runSimulation();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop -start);
    cout << "Time in miliseconds: " << duration.count() << endl;

    return EXIT_SUCCESS;


}

