//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_STAGGEREDGRID_H
#define CODE_NUMSIM_STAGGEREDGRID_H


#include "../Array2D/Array2D.h"
#include "../Array2D/FieldVariable.h"
#include "../Array2D/FieldVector.h"

class StaggeredGrid {
protected:
    const std::array<int, 2> nCells_;
    const std::array<double, 2> meshWidth_;
    int nVelo_;
    FieldVariable u_;
    FieldVariable v_;
    FieldVariable p_;
    FieldVariable rho_;
    FieldVariable t_;
    FieldVector f_;
    FieldVector feq_;
    FieldVector ftmp_;
    std::vector<double> w_;
    std::vector<std::array<double,2>> c_;
    std::vector<std::array<double,2>> e_;

    //todo: wie werden u,v,p usw. gesetzt?

public:
    StaggeredGrid(std::array<int, 2> nCellsBoundary, int nVelo, std::array<double, 2> meshWidth); //todo am Ende wieder nach oben nach PROTECTED verschieben

    const std::array<double, 2> meshWidth() const;

    const std::array<int, 2> nCells() const;

    const FieldVariable &u() const;

    const FieldVariable &v() const;

    const FieldVariable &p() const;

    const FieldVariable &t() const;

    double u(int i, int j) const;

    double &u(int i, int j);

    double v(int i, int j) const;

    double &v(int i, int j);

    double p(int i, int j) const;

    double &p(int i, int j);

    double rho(int i, int j) const;

    double &rho(int i, int j);

    double t(int i, int j) const;

    double &t(int i, int j);

    double f(int i, int j, int k) const;

    double &f(int i, int j, int k);

    double feq(int i, int j, int k) const;

    double &feq(int i, int j, int k);

    double ftmp(int i, int j, int k) const;

    double &ftmp(int i, int j, int k);

    double w (int k) const;

    double &w (int k);

    double c (int k, int l) const;

    double &c (int k, int l);

    double e (int k, int l) const;

    double &e (int k, int l);

    double dx() const;

    double dy() const;

    int iBegin() const;

    int iEnd() const;

    int jBegin() const;

    int jEnd() const;

    int kBegin() const;

    int kEnd() const;
};


#endif //CODE_NUMSIM_STAGGEREDGRID_H
