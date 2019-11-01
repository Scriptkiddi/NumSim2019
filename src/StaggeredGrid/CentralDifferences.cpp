//
// Created by Julia Pelzer on 26.10.2019.
//

#include <array>
#include "CentralDifferences.h"

CentralDifferences::CentralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth) : Discretization(
        nCells,
        meshWidth) {

}

double CentralDifferences::computeDu2Dx(int i, int j) const {
    return 0;
}

double CentralDifferences::computeDv2Dy(int i, int j) const {
    return 0;
}

double CentralDifferences::computeDuvDx(int i, int j) const {
    return 0;
}

double CentralDifferences::computeDuvDy(int i, int j) const {
    return 0;
}

