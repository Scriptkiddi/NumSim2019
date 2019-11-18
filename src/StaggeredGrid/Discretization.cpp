//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Discretization.h"
#include "StaggeredGrid.h"
#include <math.h>

Discretization::Discretization(std::array<int, 2> nCells, std::array<double, 2> meshWidth,std::shared_ptr<Partitioning> partitioning, FieldVariable u,
                               FieldVariable v, FieldVariable p, FieldVariable f, FieldVariable g, FieldVariable rhs)
        : StaggeredGrid(nCells,
                        meshWidth,
                        partitioning,
                        u,
                        v,
                        p,
                        f,
                        g,
                        rhs) {

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