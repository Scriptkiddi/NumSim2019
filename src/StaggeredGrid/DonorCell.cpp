//
// Created by Julia Pelzer on 26.10.2019.
//

#include "DonorCell.h"
#include <cstdlib>
#include <cmath>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha, double gamma) : Discretization(nCells,
                                                                                                                meshWidth,
                                                                                                                gamma),
                                                                                                 alpha_(alpha) {

}

double DonorCell::computeDu2Dx(int i, int j) const {

    return 1 / dx() * (pow((u(i, j) + u(i + 1, j)) * 0.5, 2) - pow((u(i - 1, j) + u(i, j)) * 0.5, 2))
           + alpha_ * 1 / dx() * (std::fabs(u(i, j) + u(i + 1, j)) * 0.5 * (u(i, j) - u(i + 1, j)) * 0.5 -
                                  std::fabs(u(i - 1, j) + u(i, j)) * 0.5 * (u(i - 1, j) - u(i, j)) * 0.5);
}

double DonorCell::computeDv2Dy(int i, int j) const {
    return 1 / dy() * (pow((v(i, j) + v(i, j + 1)) * 0.5, 2) - pow((v(i, j - 1) + v(i, j)) * 0.5, 2))
           + alpha_ * 1 / dy() * (std::fabs(v(i, j) + v(i, j + 1)) * 0.5 * (v(i, j) - v(i, j + 1)) *0.5 -
                                  std::fabs(v(i, j - 1) + v(i, j)) * 0.5 * (v(i, j - 1) - v(i, j)) * 0.5);
}

double DonorCell::computeDuvDx(int i, int j) const {
    return 1 / dx() * ((u(i, j) + u(i, j + 1)) / 2 * (v(i, j) + v(i + 1, j)) / 2 -
                       (u(i - 1, j) + u(i - 1, j + 1)) / 2 * (v(i - 1, j) + v(i, j)) / 2)
           + alpha_ * 1 / dx() * (std::fabs(u(i, j) + u(i, j + 1)) / 2 * (v(i, j) - v(i + 1, j)) / 2 -
                                  std::fabs(u(i - 1, j) + u(i - 1, j + 1)) / 2 * (v(i - 1, j) - v(i, j)) / 2);
}
//todo Problem mit Konvertierung von double zu int?? andere Funktion als abs?? oder doch nicht?

double DonorCell::computeDuvDy(int i, int j) const {
    return 1 / dy() * ((v(i, j) + v(i + 1, j)) / 2 * (u(i, j) + u(i, j + 1)) / 2 -
                       (v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) + u(i, j)) / 2)
           + alpha_ * 1 / dy() * (std::fabs(v(i, j) + v(i + 1, j)) / 2 * (u(i, j) - u(i, j + 1)) / 2 -
                                  std::fabs(v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) - u(i, j)) / 2);
}
