//
// Created by Julia Pelzer on 26.10.2019.
//

#include "DonorCell.h"

DonorCell::DonorCell(std::array<int, 2> nCells, std::array<double, 2> meshWidth, double alpha) {

}

double DonorCell::computeDu2Dx(int i, int j) const {
    return 0;
}

double DonorCell::computeDv2Dy(int i, int j) const {
    return 0;
}

double DonorCell::computeDuvDx(int i, int j) const {
    return 0;
}

double DonorCell::computeDuvDy(int i, int j) const {
    return 0;
}
