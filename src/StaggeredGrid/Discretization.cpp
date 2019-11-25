//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Discretization.h"
#include "StaggeredGrid.h"
#include <math.h>

//todo skript seite 10 ff?

Discretization::Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth) : StaggeredGrid(nCells,
                                                                                                           meshWidth) {

}

double Discretization::computeD2uDx2(int i, int j) const {
    return (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / pow(dx(), 2);
}

double Discretization::computeD2uDy2(int i, int j) const {
    return (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) / pow(dy(), 2);
}

double Discretization::computeD2vDx2(int i, int j) const {
    return (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) / pow(dx(), 2);
}

double Discretization::computeD2vDy2(int i, int j) const {
    return (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1)) / pow(dy(), 2);
}

double Discretization::computeDpDx(int i, int j) const {
    return (p(i + 1, j) - p(i, j)) / dx();
}

double Discretization::computeDpDy(int i, int j) const {
    return (p(i, j + 1) - p(i, j)) / dy();
}

int Discretization::computeD2TDx2(int i, int j) {
    return (t(i + 1, j) - 2 * t(i, j) + t(i - 1, j)) / pow(dx(), 2);
}

int Discretization::computeD2TDy2(int i, int j) {
    return (t(i, j + 1) - 2 * t(i, j) + t(i, j - 1)) / pow(dy(), 2);
}

double Discretization::computeDutDx(int i, int j) {
    // todo are we using gamma?
    return 1 / dx() *
           ((u(i, j) * (t(i + 1, j) + t(i, j)) / 2 - u(i - 1, j) * (t(i, j) + t(i - 1, j)) / 2) +
            abs(u(i, j)) * (t(i, j) - t(i + 1, j)) / 2 - abs(u(i - 1, j)) * (t(i - 1, j) - t(i, j)) / 2);

}

double Discretization::computeDvtDy(int i, int j) {
    // todo are we using gamma?
    return 1 / dy() *
           ((v(i, j) * (t(i + 1, j) + t(i, j)) / 2 - v(i - 1, j) * (t(i, j) + t(i - 1, j)) / 2) +
            abs(v(i, j)) * (t(i, j) - t(i + 1, j)) / 2 - abs(v(i - 1, j)) * (t(i - 1, j) - t(i, j)) / 2);

}
