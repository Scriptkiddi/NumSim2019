//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_DISCRETIZATION_H
#define CODE_NUMSIM_DISCRETIZATION_H

#include "StaggeredGrid.h"

class Discretization : public StaggeredGrid {
public:
    Discretization (std::array< int, 2 > nCells, std::array< double, 2 > meshWidth);
    virtual double 	computeDu2Dx (int i, int j) const =0;
    virtual double 	computeDv2Dy (int i, int j) const =0;
    virtual double 	computeDuvDx (int i, int j) const =0;
    virtual double 	computeDuvDy (int i, int j) const =0;
    virtual double 	computeD2uDx2 (int i, int j) const; //todo
    virtual double 	computeD2uDy2 (int i, int j) const; //todo
    virtual double 	computeD2vDx2 (int i, int j) const; //todo
    virtual double 	computeD2vDy2 (int i, int j) const; //todo
    virtual double 	computeDpDx (int i, int j) const; //todo
    virtual double 	computeDpDy (int i, int j) const; //todo
    int computeD2TDx2(int i, int j);

    int computeD2TDy2(int i, int j);

    double computeDutDx(int i, int j);

    double computeDvtDy(int i, int j);
};


#endif //CODE_NUMSIM_DISCRETIZATION_H
