//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Discretization.h"
#include "StaggeredGrid.h"

//todo skript seite 10 ff?

Discretization::Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth) : StaggeredGrid(nCells_, meshWidth_){

}

double Discretization::computeD2uDx2(int i, int j) const {
    return (p(i+1,j)-p(i,j))/dx();
}

double Discretization::computeD2uDy2(int i, int j) const {
    return 0;
}

double Discretization::computeD2vDx2(int i, int j) const {
    return 0;
}

double Discretization::computeD2vDy2(int i, int j) const {
    return 0;
}

double Discretization::computeDpDx(int i, int j) const {
    return 0;
}

double Discretization::computeDpDy(int i, int j) const {
    return 0;
}