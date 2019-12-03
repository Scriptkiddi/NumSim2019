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
    //left bound with corners
    int i_low = discretization_.get()->pIBegin()-1;
    for (int j = discretization_.get()->pJBegin()-1 ; j <= discretization_.get()->pJEnd()+1 ; j++) {
        if(!geometry_.get()->isFluid(i_low,j) && geometry_.get()->isFluid(i_low + 1,j)){
           /* if(geometry_.get()->pLeftBoundType[j] == "Neumann"){
                discretization_.get()->p(i_low, j) = discretization_.get()->p(i_low + 1, j) - discretization_.get()->dx() * geometry_.get()->pLeftBoundValue[j];
            }else{
                discretization_.get()->p(i_low, j) = 2 * geometry_.get()->pLeftBoundValue[j] - discretization_.get()->p(i_low + 1, j); 
            }*/
            if(geometry_.get()->pressure(i_low,j).first == "PR"){
                discretization_.get()->p(i_low,j) = 2 * geometry_.get()->pressure(i_low,j).second[0] - discretization_.get()->p(i_low+1,j); //2p_in - p(1,j)
            }else{
                discretization_.get()->p(i_low,j) = discretization_.get()->p(i_low+1,j);
            }
        }
    }

    //right bound with corners
    int i_high = discretization_.get()->pIEnd()+1;
    for (int j = discretization_.get()->pJBegin()-1 ; j <= discretization_.get()->pJEnd()+1 ; j++) {
        if(!geometry_.get()->isFluid(i_high,j) && geometry_.get()->isFluid(i_high - 1,j)){
            /*if(geometry_.get()->pRightBoundType[j] == "Neumann"){
                discretization_.get()->p(i_high, j) = discretization_.get()->p(i_high - 1, j) - discretization_.get()->dx() * geometry_.get()->pRightBoundValue[j];
            }else{
                discretization_.get()->p(i_high, j) = 2 * geometry_.get()->pRightBoundValue[j] - discretization_.get()->p(i_high + 1, j); 
            }*/
            if(geometry_.get()->pressure(i_high,j).first == "PR"){
                discretization_.get()->p(i_high,j) = 2 * geometry_.get()->pressure(i_high,j).second[0] - discretization_.get()->p(i_high-1,j); //2p_in - p(1,j)
            }else{
                discretization_.get()->p(i_high,j) = discretization_.get()->p(i_high-1,j);
            }
        }
    }

    //bottom bound without corners
    int j_low = discretization_.get()->pJBegin()-1;
    for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
        if(!geometry_.get()->isFluid(i,j_low) && geometry_.get()->isFluid(i, j_low + 1)){
            /*if(geometry_.get()->pBottomBoundType[i] == "Neumann"){
                discretization_.get()->p(i, j_low) = discretization_.get()->p(i, j_low + 1) - discretization_.get()->dy() * geometry_.get()->pBottomBoundValue[i];
            }else{
                discretization_.get()->p(i, j_low) = 2 * geometry_.get()->pBottomBoundValue[i] - discretization_.get()->p(i, j_low + 1);
            }*/
            if(geometry_.get()->pressure(i,j_low).first == "PR"){
                discretization_.get()->p(i,j_low) = 2 * geometry_.get()->pressure(i,j_low).second[0] - discretization_.get()->p(i,j_low+1); //2p_in - p(1,j)
            }else{
                discretization_.get()->p(i,j_low) = discretization_.get()->p(i,j_low+1);
            }
        }
    }

    //top bound without corners
    int j_high = discretization_.get()->pJEnd()-1;
    for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
        if(!geometry_.get()->isFluid(i,j_high) && geometry_.get()->isFluid(i, j_high - 1)){
            /*if(geometry_.get()->pTopBoundType[i] == "Neumann"){
                discretization_.get()->p(i, j_high) = discretization_.get()->p(i, j_high - 1) - discretization_.get()->dy() * geometry_.get()->pTopBoundValue[i];
            }else{
                discretization_.get()->p(i, j_high) = 2 * geometry_.get()->pTopBoundValue[i] - discretization_.get()->p(i, j_high - 1);
            }*/
            if(geometry_.get()->pressure(i,j_high).first == "PR"){
                discretization_.get()->p(i,j_high) = 2 * geometry_.get()->pressure(i,j_high).second[0] - discretization_.get()->p(i,j_high-1); //2p_in - p(1,j)
            }else{
                discretization_.get()->p(i,j_high) = discretization_.get()->p(i,j_high-1);
            }
        }
    }

    //inner cells (obstacles)
    for(int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++){
        for(int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++){
            if(!geometry_.get()->isFluid(i,j)){
                int N = 0;
                double pTmp = 0;

                //check left and right neighbours
                if(geometry_.get()->isFluid(i - 1, j)){
                    pTmp += discretization_.get()->p(i - 1, j);
                    N++;
                }else if(geometry_.get()->isFluid(i + 1, j)){
                    pTmp += discretization_.get()->p(i + 1, j);
                    N++;
                }
                //check top and bottom neighbours
                if(geometry_.get()->isFluid(i, j - 1)){
                    pTmp += discretization_.get()->p(i, j - 1);
                    N++;
                }else if(geometry_.get()->isFluid(i, j + 1)){
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
