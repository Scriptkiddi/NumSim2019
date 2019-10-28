//
// Created by Julia Pelzer on 26.10.2019.
//

#include "DonorCell.h"
#include <cstdlib>
#include <cmath>

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha) : Discretization(nCells,
                                                                                                                meshWidth),
                                                                                                 alpha_(alpha) {

}

double DonorCell::computeDu2Dx(int i, int j) const {

    return 1 / dx() * (sqrt((u_(i, j) + u(i + 1, j)) / 2) - sqrt((u(i - 1, j) + u(, j)) / 2))
           + alpha_ * 1 / dx() * (std::abs(u(i, j) + u(i + 1, j)) / 2 * (u(i, j) - u(i + 1, j)) / 2 -
                                  std::abs(u(i - 1, j) + u(i, j)) / 2 * (u(i - 1, j) - u(i, j)) / 2);
}

double DonorCell::computeDv2Dy(int i, int j) const {
    return 1 / dy() * (sqrt((v(i, j) + v(i + 1, j)) / 2) - sqrt((v(i - 1, j) + u(, j)) / 2))
           + alpha_ * 1 / dy() * (std::abs(v(i, j) + v(i + 1, j)) / 2 * (v(i, j) - v(i + 1, j)) / 2 -
                                  std::abs(v(i - 1, j) + v(i, j)) / 2 * (v(i - 1, j) - v(i, j)) / 2);
}

double DonorCell::computeDuvDx(int i, int j) const {
    return 1 / dx() * ((u(i, j) + u(i + 1, j)) / 2 * (v(i, j) + v(i, j + 1)) / 2 -
                       (u(i, j - 1) + u(i + 1, j - 1)) / 2 * (v(i, j - 1) + v(i, j)) / 2)
           + alpha_ * 1 / dx() * (abs(u(i, j) + u(i + 1, j)) / 2 * (v(i, j) - v(i, j + 1)) / 2 -
                                abs(u(i, j - 1) + u(i + 1, j - 1)) / 2 * (v(i, j - 1) - v(i, j)) / 2);
}

double DonorCell::computeDuvDy(int i, int j) const {
    return 1 / dy() * ((v(i, j) + v(i + 1, j)) / 2 * (u(i, j) + u(i, j + 1)) / 2 -
                       (v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) + u(i, j)) / 2)
           + alpha_ * 1 / dy() * (abs(v(i, j) + v(i + 1, j)) / 2 * (u(i, j) - u(i, j + 1)) / 2 -
                                abs(v(i, j - 1) + v(i + 1, j - 1)) / 2 * (u(i, j - 1) - u(i, j)) / 2);
}
