//
// Created by Julia Pelzer on 26.10.2019.
//

#include "PressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, std::shared_ptr<Geometry> geometry,double epsilon,
                               int maximumNumberOfIterations) :
        discretization_(discretization),
        geometry_(geometry),
        epsilon_(epsilon),
        maximumNumberOfIterations_(maximumNumberOfIterations) {

}

void PressureSolver::setBoundaryValues() {

    //outer bounds
    //TODO: only homogeneous Neumann BC for now, change if necessary
    //left bound without corners
    int i_low = discretization_.get()->pIBegin()-1;
    for (int j = discretization_.get()->pJBegin() ; j <= discretization_.get()->pJEnd() ; j++) {
        if(geometry_.get()->pLeftBoundType[j] == "Neumann"){
            discretization_.get()->p(i_low, j) = discretization_.get()->p(i_low + 1, j);
        }else{
            discretization_.get()->p(i_low, j) = 2 * geometry_.get()->pLeftBoundValue[j] - discretization_.get()->p(i_low + 1, j);
        }
    }

    //right bound without corners
    int i_high = discretization_.get()->pIEnd()+1;
    for (int j = discretization_.get()->pJBegin() ; j <= discretization_.get()->pJEnd() ; j++) {
        if(geometry_.get()->pRightBoundType[j] == "Neumann"){
            discretization_.get()->p(i_high, j) = discretization_.get()->p(i_high - 1, j);
        }else{
            discretization_.get()->p(i_high, j) = 2 * geometry_.get()->pRightBoundValue[j] - discretization_.get()->p(i_high + 1, j);
        }
    }

    //bottom bound with corners
    int j_low = discretization_.get()->pJBegin()-1;
    for (int i = discretization_.get()->pIBegin()-1; i <= discretization_.get()->pIEnd()+1; i++) {
        if(geometry_.get()->pBottomBoundType[i] == "Neumann"){
            discretization_.get()->p(i, j_low) = discretization_.get()->p(i, j_low + 1);
        }else{
            discretization_.get()->p(i, j_low) = 2 * geometry_.get()->pBottomBoundValue[i] - discretization_.get()->p(i, j_low + 1);
        }
    }

    //top bound with corners
    int j_high = discretization_.get()->pJEnd()-1;
    for (int i = discretization_.get()->pIBegin()-1; i <= discretization_.get()->pIEnd()+1; i++) {
        if(geometry_.get()->pTopBoundType[i] == "Neumann"){
            discretization_.get()->p(i, j_high) = discretization_.get()->p(i, j_high - 1);
        }else{
            discretization_.get()->p(i, j_high) = 2 * geometry_.get()->pTopBoundValue[i] - discretization_.get()->p(i, j_high - 1);
        }
    }

    //inner cells (obstacles)
    for(int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++){
        for(int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++){
            if(!geometry_.get()->isFluid(i,j)){
                int N = 0;
                double pTmp = 0;

                //check left neighbour
                if(geometry_.get()->isFluid(i - 1, j)){
                    pTmp += discretization_.get()->p(i - 1, j);
                    N++;
                }
                //check right neighbour
                if(geometry_.get()->isFluid(i + 1, j)){
                    pTmp += discretization_.get()->p(i + 1, j);
                    N++;
                }
                //check bottom neighbour
                if(geometry_.get()->isFluid(i, j - 1)){
                    pTmp += discretization_.get()->p(i, j - 1);
                    N++;
                }
                //check top neighbour
                if(geometry_.get()->isFluid(i, j + 1)){
                    pTmp += discretization_.get()->p(i, j + 1);
                    N++;
                }
                //assign average value
                if(N > 0){
                    discretization_.get()->p(i,j) = pTmp / N;
                }
            }
        }
    }

}
