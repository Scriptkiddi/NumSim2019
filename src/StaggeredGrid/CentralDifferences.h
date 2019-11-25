//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_CENTRALDIFFERENCES_H
#define CODE_NUMSIM_CENTRALDIFFERENCES_H

#include "Discretization.h"


class CentralDifferences : public Discretization {
public:
    CentralDifferences(std::array<int, 2>
                       nCells,
                       std::array<double, 2> meshWidth, double gamma
    );

    virtual double computeDu2Dx(int i, int j) const;

    virtual double computeDv2Dy(int i, int j) const;

    virtual double computeDuvDx(int i, int j) const;

    virtual double computeDuvDy(int i, int j) const;
};


#endif //CODE_NUMSIM_CENTRALDIFFERENCES_H
