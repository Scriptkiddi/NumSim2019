//
// Created by Julia Pelzer on 26.10.2019.
//

#include <array>
#include "CentralDifferences.h"
#include <cmath>

CentralDifferences::CentralDifferences(std::array<int, 2> nCells, std::array<double, 2> meshWidth,
                                       std::shared_ptr<Partitioning> partitioning, FieldVariable u,
                                       FieldVariable v, FieldVariable p, FieldVariable f, FieldVariable g, 
                                       FieldVariable rhs) : Discretization(nCells, meshWidth, partitioning, u, v,
                                                                           p,
                                                                           f,
                                                                           g,
                                                                           rhs) {
}

double CentralDifferences::computeDu2Dx(int i, int j) const {
    return 1 / dx() * (pow((u(i, j) + u(i + 1, j)) / 2, 2) - pow((u(i - 1, j) + u(i, j)) / 2, 2));
}

double CentralDifferences::computeDv2Dy(int i, int j) const {
    return 1 / dy() * (pow((v(i, j) + v(i, j + 1)) / 2, 2) - pow((v(i, j - 1) + v(i, j)) / 2, 2));
}

double CentralDifferences::computeDuvDx(int i, int j) const {
    return 1 / dx()
           * ((u(i, j) + u(i, j + 1)) / 2
              * (v(i, j) + v(i + 1, j)) / 2
              -
              (u(i - 1, j) + u(i - 1, j + 1)) / 2
              * (v(i - 1, j) + v(i, j)) / 2);
}

double CentralDifferences::computeDuvDy(int i, int j) const {
    return 1 / dy() * ((v(i, j) + v(i + 1, j)) / 2 * (u(i, j) + u(i, j + 1)) / 2 -
                       (v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) + u(i, j)) / 2);
}

