//
// Created by Fabian on 26.10.2019.
//
#include <iostream>
#include "Computation/Computation.h"
#include "Computation/Computation_solid.h"
#include <chrono>

using namespace std::chrono;

int main(int argc, char *argv[]) {
    // Todo parse number of processes
    if (argc == 1 ){
        std::cout << "Please pass settings filename" << std::endl;
        return EXIT_FAILURE;
    }
    
    // TODO WIE AUFRUF VON MAIN_SOLID? NOTWENDIG, DASS EIGENE DATEI ODER HIER EINPFLEGEN?
    
    //Computation computation;
    Computation_solid computation_solid;

    //computation.initialize(argc, argv);
    computation_solid.initialize(argc, argv);
    
    auto start = high_resolution_clock::now();
    //computation.runSimulation();
    computation_solid.runSimulation();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop -start);
    cout << "Time in miliseconds: " << duration.count() << endl;

    return EXIT_SUCCESS;


}

