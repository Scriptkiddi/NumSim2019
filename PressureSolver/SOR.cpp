//
// Created by Julia Pelzer on 26.10.2019.
//

#include "SOR.h"

SOR::SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega):
    PressureSolver(discretization,epsilon,maximumNumberOfIterations),
    omega(omega){

}

void SOR::solve() {
    int a = 2;
    while(a <2){
        a++;
        for(int j =0; j<1; j++){
            for(int i=0; i<1; i++){
                //discretization_.get()->p(i,j) = discretization_.get()->p(i,j) + ;
            }
        }
    }
}
