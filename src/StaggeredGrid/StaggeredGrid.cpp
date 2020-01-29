//
// Created by Julia Pelzer on 26.10.2019.
//

#include "StaggeredGrid.h"
#include "../Array2D/Array2D.h"
#include "../Array2D/FieldVariable.h"
#include "../Array2D/FieldVector.h"
#include <cmath>

StaggeredGrid::StaggeredGrid(std::array<int, 2> nCellsBoundary, int nVelo, std::array<double, 2> meshWidth) :
        meshWidth_(meshWidth),
        nCells_(nCellsBoundary),
        nVelo_(nVelo),
        u_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        v_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        p_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        rho_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        t_(nCellsBoundary, {-0.5 * meshWidth[0], -0.5 * meshWidth[1]}, meshWidth),
        f_(nCellsBoundary, nVelo, meshWidth),
        feq_(nCellsBoundary, nVelo, meshWidth),
        ftmp_(nCellsBoundary, nVelo, meshWidth){
            w_.resize(nVelo,0.0);
            c_.resize(nVelo,{0.0,0.0});
            e_.resize(nVelo,{0.0,0.0});
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

double StaggeredGrid::rho(int i, int j) const {
    return rho_(i, j);
}

double &StaggeredGrid::rho(int i, int j) {
    return rho_(i, j);
}

double &StaggeredGrid::t(int i, int j) {
    return t_(i, j);
}

double StaggeredGrid::f(int i, int j, int k) const {
    return f_(i, j, k);
}

double &StaggeredGrid::f(int i, int j, int k) {
    return f_(i, j, k);
}

double StaggeredGrid::feq(int i, int j, int k) const {
    return feq_(i, j, k);
}

double &StaggeredGrid::feq(int i, int j, int k) {
    return feq_(i, j, k);
}

double StaggeredGrid::ftmp(int i, int j, int k) const {
    return ftmp_(i, j, k);
}

double &StaggeredGrid::ftmp(int i, int j, int k) {
    return ftmp_(i, j, k);
}

double &StaggeredGrid::w(int k){
    return w_[k];
}

double StaggeredGrid::w(int k) const{
    return w_[k];
}

double &StaggeredGrid::c(int k, int l){
    return c_[k][l];
}

double StaggeredGrid::c(int k, int l) const{
    return c_[k][l];
}

double &StaggeredGrid::e(int k, int l){
    return e_[k][l];
}

double StaggeredGrid::e(int k, int l) const{
    return e_[k][l];
}

double StaggeredGrid::dx() const {
    return meshWidth_[0];
}

double StaggeredGrid::dy() const {
    return meshWidth_[1];
}

int StaggeredGrid::iBegin() const {
    return 1;
}

int StaggeredGrid::iEnd() const {
    return nCells_[0]-2;
}

int StaggeredGrid::jBegin() const {
    return 1;
}

int StaggeredGrid::jEnd() const {
    return nCells_[1]-2;
}

int StaggeredGrid::kBegin() const {
    return 0;
}

int StaggeredGrid::kEnd() const {
    return nVelo_ - 1;
}