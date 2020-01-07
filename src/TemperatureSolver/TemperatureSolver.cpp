//
// Created by Julia Pelzer on 26.10.2019.
//

#include "TemperatureSolver.h"


TemperatureSolver::TemperatureSolver(std::shared_ptr<Discretization> discretization, std::shared_ptr<Geometry> geometry,double epsilon,
                               int maximumNumberOfIterations, double dt, double heatDiffusivity) :
        discretization_(discretization),
        geometry_(geometry),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations),
        dt_(dt),
        alpha_(heatDiffusivity) {

}
void TemperatureSolver::applyBoundaryValuesTemperature() {
    // help variables
    int i_low;
    int i_high;
    int j_low;
    int j_high;

    //outer bounds
    //left bound without corners
    i_low = discretization_.get()->tIBegin() - 1;
    for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
        if (!geometry_.get()->isFluid(i_low, j)) {
            if (geometry_.get()->get_temperature(i_low, j).first == "TN" || geometry_.get()->get_temperature(i_low, j).first == "TPN") {
                discretization_.get()->t(i_low, j) = discretization_.get()->t(i_low + 1, j) -
                                                     discretization_.get()->dx() * alpha_ *
                                                     geometry_.get()->get_temperature(i_low, j).second[0];
            } else { //"TD"
                discretization_.get()->t(i_low, j) = 2 * geometry_.get()->get_temperature(i_low, j).second[0] -
                                                     discretization_.get()->t(i_low + 1, j);
            }
        }
    }

    //right bound without corners
    i_high = discretization_.get()->tIEnd() + 1;
    for (int j = discretization_.get()->tJBegin() - 1; j <= discretization_.get()->tJEnd() + 1; j++) {
        if (!geometry_.get()->isFluid(i_high, j)) {
            if ( geometry_.get()->get_temperature(i_high, j).first == "TN" || geometry_.get()->get_temperature(i_high, j).first == "TPN") {
                discretization_.get()->t(i_high, j) = discretization_.get()->t(i_high - 1, j) -
                                                      discretization_.get()->dx() * alpha_ *
                                                      geometry_.get()->get_temperature(i_high, j).second[0];
            } else {
                discretization_.get()->t(i_high, j) = 2 * geometry_.get()->get_temperature(i_high, j).second[0] -
                                                      discretization_.get()->t(i_high - 1, j);
            }
        }
    }

    //bottom bound with corners
    j_low = discretization_.get()->tJBegin() - 1;
    for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
        if (!geometry_.get()->isFluid(i, j_low)) {
            if (geometry_.get()->get_temperature(i, j_low).first == "TN" || geometry_.get()->get_temperature(i, j_low).first == "TPN") {
                discretization_.get()->t(i, j_low) = discretization_.get()->t(i, j_low + 1) -
                                                     discretization_.get()->dy() * alpha_ *
                                                     geometry_.get()->get_temperature(i, j_low).second[0];
            } else {
                discretization_.get()->t(i, j_low) = 2 * geometry_.get()->get_temperature(i, j_low).second[0] -
                                                     discretization_.get()->t(i, j_low + 1);
            }
        }
    }

    //top bound with corners
    j_high = discretization_.get()->tJEnd() + 1;
    for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
        if (!geometry_.get()->isFluid(i, j_high)) {
            if (geometry_.get()->get_temperature(i, j_high).first == "TN" || geometry_.get()->get_temperature(i, j_high).first == "TPN") {
                discretization_.get()->t(i, j_high) = discretization_.get()->t(i, j_high - 1) -
                                                      discretization_.get()->dy() * alpha_ *
                                                      geometry_.get()->get_temperature(i, j_high).second[0];
            } else {
                discretization_.get()->t(i, j_high) = 2 * geometry_.get()->get_temperature(i, j_high).second[0] -
                                                      discretization_.get()->t(i, j_high - 1);
            }
        }
    }

    //inner cells (obstacles)
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if (!geometry_.get()->isFluid(i, j)) {

                if (geometry_.get()->get_temperature(i,j).first == "TN" || geometry_.get()->get_temperature(i,j).first == "TPN"){
                //Neumann-BC
                    int N = 0; //number of neighbouring fluid cells
                    double tTmp = 0;

                    if (geometry_.get()->isFluid(i - 1,j)){
                        tTmp += discretization_.get()->t(i - 1,j) - discretization_.get()->dx() / alpha_ * geometry_.get()->get_temperature(i,j).second[0];
                        N++;
                    } 
                    if (geometry_.get()->isFluid(i + 1, j)) {
                        tTmp += discretization_.get()->t(i + 1,j) - discretization_.get()->dx() / alpha_ * geometry_.get()->get_temperature(i,j).second[0];
                        N++;
                    }
                    if (geometry_.get()->isFluid(i, j - 1)) {
                        tTmp += discretization_.get()->t(i, j - 1) - discretization_.get()->dx() / alpha_ * geometry_.get()->get_temperature(i,j).second[0];
                        N++;
                    }
                    if (geometry_.get()->isFluid(i, j + 1)) {
                        tTmp += discretization_.get()->t(i ,j + 1) - discretization_.get()->dx() / alpha_ * geometry_.get()->get_temperature(i,j).second[0];
                        N++;
                    }
                     //assign average value
                    if (N > 0) {
                        discretization_.get()->t(i, j) = tTmp / N;
                    }

                } else if(geometry_.get()->get_temperature(i,j).first == "TD" || geometry_.get()->get_temperature(i,j).first == "TPD"){
                //Dirichlet-BC
                    int N = 0; //number of neighbouring fluid cells
                    double tTmp = 0;

                    if (geometry_.get()->isFluid(i - 1, j)) {
                        tTmp += 2 * geometry_.get()->get_temperature(i,j).second[0] - discretization_.get()->t(i - 1, j);
                        N++;
                    }
                    if (geometry_.get()->isFluid(i + 1, j)) {
                        tTmp += 2 * geometry_.get()->get_temperature(i,j).second[0] - discretization_.get()->t(i + 1, j);
                        N++;
                    }
                    if (geometry_.get()->isFluid(i, j - 1)) {
                        tTmp += 2 * geometry_.get()->get_temperature(i,j).second[0] - discretization_.get()->t(i, j - 1);
                        N++;
                    }
                    if (geometry_.get()->isFluid(i, j + 1)) {
                        tTmp += 2 * geometry_.get()->get_temperature(i,j).second[0] - discretization_.get()->t(i, j + 1);
                        N++;
                    }
                    //assign average value
                    if (N > 0) {
                        discretization_.get()->t(i, j) = tTmp / N;
                    }

                } else {
                //Default-BC (homogenous Neumann)
                    int N = 0; //number of neighbouring fluid cells
                    double tTmp = 0;

                    if (geometry_.get()->isFluid(i - 1, j)) {
                        tTmp += discretization_.get()->t(i - 1, j);
                        N++;
                    }
                    if (geometry_.get()->isFluid(i + 1, j)) {
                        tTmp += discretization_.get()->t(i + 1, j);
                        N++;
                    }
                    if (geometry_.get()->isFluid(i, j - 1)) {
                        tTmp += discretization_.get()->t(i, j - 1);
                        N++;
                    }
                    if (geometry_.get()->isFluid(i, j + 1)) {
                        tTmp += discretization_.get()->t(i, j + 1);
                        N++;
                    }
                    //assign average value
                    if (N > 0) {
                        discretization_.get()->t(i, j) = tTmp / N;
                    }
                }
            }
        }
    }
}
