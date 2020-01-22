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
    FieldVariable u_;
    FieldVariable v_;
    FieldVariable p_;
    FieldVariable t_;

    FieldVector f_;
    FieldVector feq_;
    FieldVector ftmp_;

    //todo: wie werden u,v,p usw. gesetzt?

public:
    StaggeredGrid(std::array<int, 2> nCellsBoundary, std::array<double, 2> meshWidth, int nVelo); //todo am Ende wieder nach oben nach PROTECTED verschieben

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

    double t(int i, int j) const;

    double &t(int i, int j);

    double dx() const;

    double dy() const;

    int uIBegin() const;

    int uIEnd() const;

    int uJBegin() const;

    int uJEnd() const;

    int vIBegin() const;

    int vIEnd() const;

    int vJBegin() const;

    int vJEnd() const;

    int pIBegin() const;

    int pIEnd() const;

    int pJBegin() const;

    int pJEnd() const;

    int tIBegin() const;

    int tIEnd() const;

    int tJBegin() const;

    int tJEnd() const;
};


#endif //CODE_NUMSIM_STAGGEREDGRID_H
