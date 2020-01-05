//
// Created by Julia Pelzer on 04.01.2020.
//

#include "GaussSeidel.h"
#include <cmath>
#include <iostream>
#include "Settings.h"

GaussSeidel::GaussSeidel(std::shared_ptr<Discretization> discretization, std::shared_ptr<Geometry> geometry, double epsilon, int maximumNumberOfIterations, double dt, double heatDiffusivity) :
        TemperatureSolver(discretization, geometry, epsilon, maximumNumberOfIterations, dt, heatDiffusivity) {

}

void GaussSeidel::solve() {
    int iter = 0;
    double eps = 1;
    double dx2 = pow(discretization_.get()->dx(), 2);
    double dy2 = pow(discretization_.get()->dy(), 2);
    while (iter <= maximumNumberOfIterations_ && eps > pow(epsilon_,2)) {
        applyBoundaryValuesTemperature();
        for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
            for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
                if(geometry_.get()->isFluid(i,j)){
                
                //GaussSeidel Temperature
                    discretization_.get()->t(i,j) = 
                        1 / (1/dt_ + 2 * alpha_ * (1/dx2 + 1/dy2)) * 
                        (
                            alpha_ *
                            (
                                 (discretization_.get()->t(i+1,j) +  discretization_.get()->t(i-1,j)) / dx2
                                 + 
                                 (discretization_.get()->t(i,j+1) +  discretization_.get()->t(i,j-1)) / dy2
                            )
                        + discretization_.get()->tOld(i,j) / dt_
                        );

                // GaussSeidel Pressure
                /*
                    discretization_.get()->p(i, j) = 
                        dx2 * dy2 / (2 * (dx2 + dy2)) * 
                        ( (discretization_.get()->p(i - 1, j) + discretization_.get()->p(i + 1, j))
                            / dx2
                            +
                            (discretization_.get()->p(i, j - 1) + discretization_.get()->p(i, j + 1)) 
                            / dy2       //- discretization_.get()->rhs(i, j)
                        );
                        */
                
                }
            }
        }
        //residuum  // maximum of last error and temperature in this timestep minus temperature in the previous timestep
        //TODO right??? DO NEXT
        eps = 0;
        for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
            for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
                if(geometry_.get()->isFluid(i,j)){
                    //eps for GaussSeidel Temperature
                    //TODO not shure if this is correct yet
                    eps = eps + pow(
                        dt_ * alpha_ * 
                            (
                                (discretization_.get()->t(i+1,j) - 2 * discretization_.get()->t(i,j) +  discretization_.get()->t(i-1,j)) / dx2
                                 + 
                                (discretization_.get()->t(i,j+1) - 2 * discretization_.get()->t(i,j) +  discretization_.get()->t(i,j-1)) / dy2
                            )

                    ,2);
                    // eps for GaussSeidel pressure
                    /*
                    eps = eps + pow(
                            discretization_->rhs(i, j)
                            - (discretization_.get()->t(i - 1, j) - 2 * discretization_.get()->t(i, j) +
                            discretization_.get()->t(i + 1, j))
                            / dx2
                            - (discretization_.get()->t(i, j - 1) - 2 * discretization_.get()->t(i, j) +
                            discretization_.get()->t(i, j + 1))
                            / dy2, 2);
                    */
                }
            }
        }
        eps = eps/geometry_.get()->nCellsFluid();
        //eps = eps / ((discretization_.get()->nCells()[0] - 2) * (discretization_.get()->nCells()[1] - 2));
        iter++;
    }
    //std::cout << "pressure solver iterations: " << iter << " eps :" << eps << " epsilon² " << pow(epsilon_,2) <<std::endl;
    applyBoundaryValuesTemperature(); //TODO da war doch mal was vonwegen falsch/richtig mit zweimal die Funktion aufrufen... 
    //TODO ...weiß leider nicht mehr wierum das richtig war. wer es weiß bitte ändern und dann kommentar löschen
}

