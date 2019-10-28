//
// Created by Fabian on 26.10.2019.
//
#include <iostream>
#include "Array2D/FieldVariable.h"

int main() {
//    std::cout << "Hello World" << std::endl;

    FieldVariable a({2, 2}, {0, 0}, {1, 1});
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= 1; j++) {
            a(i,j) = i+j;
            std::cout << a(i, j) << std::endl;
        }
    }
}
