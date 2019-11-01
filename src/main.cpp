//
// Created by Fabian on 26.10.2019.
//
#include <iostream>
#include "Computation.h"

int main(int argc, char *argv[]) {
    if (argc == 1 ){
        std::cout << "Please pass settings filename" << std::endl;
        return EXIT_FAILURE;
    }

    Computation computation;
    computation.initialize(argc, argv);


    return EXIT_SUCCESS;


}

