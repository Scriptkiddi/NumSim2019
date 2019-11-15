//
// Created by Julia Pelzer on 26.10.2019.
//

#include "StaggeredGrid.h"
#include "../Array2D/Array2D.h"
#include "../Array2D/FieldVariable.h"

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCellsBoundary, std::array<double, 2> meshWidth, std::shared_ptr<Partitioning> partitioning, FieldVariable u,
                             FieldVariable v, FieldVariable p, FieldVariable f, FieldVariable g, FieldVariable rhs) :
        meshWidth_(meshWidth),
        partition_(partitioning),
        nCells_(nCellsBoundary), u_(u), v_(v), p_(p), rhs_(rhs), g_(g), f_(f) {
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
    if(partition_.get()->getRankOfLeftNeighbour() == -1) {
        return 1;
    }else{
        return 2;
    }
}

int StaggeredGrid::uIEnd() const {
    if(partition_.get()->getRankOfRightNeighbour() == -1) {
        return nCells_[0] - 1;
    }else{
        return nCells_[0] - 2;
    }
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
    if(partition_.get()->getRankOfBottomNeighbour() == -1) {
        return 1;
    }else{
        return 2;
    }
}

int StaggeredGrid::vJEnd() const {
    if(partition_.get()->getRankOfTopNeighbour() == -1) {
        return nCells_[1] - 1;
    }else{
        return nCells_[1] - 2;
    }
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