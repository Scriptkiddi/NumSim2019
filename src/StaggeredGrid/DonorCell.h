//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_DONORCELL_H
#define CODE_NUMSIM_DONORCELL_H

#include "Discretization.h"


class DonorCell : public Discretization {
private:
    double alpha_;
public:
    DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha);

    virtual double computeDu2Dx(int i, int j) const;

    virtual double computeDv2Dy(int i, int j) const;

    virtual double computeDuvDx(int i, int j) const;

    virtual double computeDuvDy(int i, int j) const;
};


#endif //CODE_NUMSIM_DONORCELL_H
