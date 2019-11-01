//
// Created by Julia Pelzer on 26.10.2019.
//

#include "StaggeredGrid.h"
#include "../Array2D/Array2D.h"
#include "../Array2D/FieldVariable.h"

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCells, std::array<double, 2> meshWidth) :
        meshWidth_(meshWidth),
        nCells_(nCells),
        u_(nCells, {1 * meshWidth[0], 0.5 * meshWidth[1]}, meshWidth),
        v_(nCells, {0.5 * meshWidth[0], 1 * meshWidth[1]}, meshWidth),
        p_(nCells, {0.5 * meshWidth[0], 0.5 * meshWidth[1]}, meshWidth),
        rhs_(nCells, {0.5 * meshWidth[0], 0.5 * meshWidth[1]}, meshWidth),
        f_(nCells, {0.5 * meshWidth[0], 0.5 * meshWidth[1]}, meshWidth),
        g_(nCells, {0.5 * meshWidth[0], 0.5 * meshWidth[1]}, meshWidth) {
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

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}

int StaggeredGrid::uIBegin() const {
    return 0; //todo evtl anpassen mit ghost cells
}

int StaggeredGrid::uIEnd() const { //one after last valid index for u in x direction
    return nCells_[0] - 1; //todo anfassen, falls wir Zahl der Zellen (N+1 vs N+2) ändern
}

int StaggeredGrid::uJBegin() const {
    return 0;
}

int StaggeredGrid::uJEnd() const { //one after last valid index for u in y direction
    return nCells_[1] - 1; //todo anfassen, falls wir Zahl der Zellen (N+1 vs N+2) ändern
}

int StaggeredGrid::vIBegin() const {
    return 0;
}

int StaggeredGrid::vIEnd() const {
    return nCells_[0] - 1;
}

int StaggeredGrid::vJBegin() const {
    return 0;
}

int StaggeredGrid::vJEnd() const {
    return nCells_[1] - 1;
}

int StaggeredGrid::pIBegin() const {
    return 0;
}

int StaggeredGrid::pIEnd() const {
    return nCells_[0] - 1;
}

int StaggeredGrid::pJBegin() const {
    return 0;
}

int StaggeredGrid::pJEnd() const {
    return nCells_[1] - 1;
}