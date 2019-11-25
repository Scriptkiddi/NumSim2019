//
// Created by Julia Pelzer on 26.10.2019.
//

#include "StaggeredGrid.h"
#include "../Array2D/Array2D.h"
#include "../Array2D/FieldVariable.h"

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCellsBoundary, std::array<double, 2> meshWidth) :
        meshWidth_(meshWidth),
        nCells_(nCellsBoundary),
        u_(nCellsBoundary, {0 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        v_(nCellsBoundary, {-0.5 * meshWidth[0], 0 * meshWidth[1]}, meshWidth),
        p_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        rhs_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        f_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        g_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
    t_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth){
    }

    const std::array<double, 2> StaggeredGrid::meshWidth() const {
    return meshWidth_;
}

const std::array<int, 2> StaggeredGrid::nCells() const {
    return nCells_;
}

const FieldVariable &StaggeredGrid::u() const {
    return u_;
}

const FieldVariable &StaggeredGrid::v() const {
    return v_;
}

const FieldVariable &StaggeredGrid::p() const {
    return p_;
}

const FieldVariable &StaggeredGrid::t() const {
    return t_;
}

double StaggeredGrid::u(int i, int j) const {
    return u_(i, j);
}

double StaggeredGrid::v(int i, int j) const {
    return v_(i, j);
}

double &StaggeredGrid::u(int i, int j) {
    return u_(i, j);
}

double &StaggeredGrid::v(int i, int j) {
    return v_(i, j);
}

double StaggeredGrid::p(int i, int j) const {
    return p_(i, j);
}

double &StaggeredGrid::p(int i, int j) {
    return p_(i, j);
}

double &StaggeredGrid::rhs(int i, int j) {
    return rhs_(i, j);
}

double &StaggeredGrid::f(int i, int j) {
    return f_(i, j);
}

double &StaggeredGrid::g(int i, int j) {
    return g_(i, j);
}

double &StaggeredGrid::t(int i, int j) {
    return t_(i, j);
}

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}
// ncells beinhaltet ghost cells
// https://numsim-exercises.readthedocs.io/en/latest/exercise1/hints.html
// 1 weil uIBegin gibt uns die erste nicht ghost cell
int StaggeredGrid::uIBegin() const {
    return 1;
}

int StaggeredGrid::uIEnd() const { //one after last valid index for u in x direction
    // minus 1 für size to index und nochmal -1 weil der wert der letzten richtigen zelle auf 0 gesetzt wird anstelle der in der ghost cell
    return nCells_[0] - 3;
}

int StaggeredGrid::uJBegin() const {
    return 1;
}

int StaggeredGrid::uJEnd() const { //one after last valid index for u in y direction
    return nCells_[1]-2;
}

int StaggeredGrid::vIBegin() const {
    return 1;
}

int StaggeredGrid::vIEnd() const {
    return nCells_[0]-2;
}

int StaggeredGrid::vJBegin() const {
    return 1;
}

int StaggeredGrid::vJEnd() const {
    return nCells_[1] - 3;
}

int StaggeredGrid::pIBegin() const {
    return 1;
}

int StaggeredGrid::pIEnd() const {
    return nCells_[0]-2;
}

int StaggeredGrid::pJBegin() const {
    return 1;
}

int StaggeredGrid::pJEnd() const {
    return nCells_[1]-2;
}

int StaggeredGrid::tIBegin() const {
    return 1;
}

int StaggeredGrid::tIEnd() const {
    return nCells_[0]-2;
}

int StaggeredGrid::tJBegin() const {
    return 1;
}

int StaggeredGrid::tJEnd() const {
    return nCells_[1]-2;
}
