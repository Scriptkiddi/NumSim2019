//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_FIELDVARIABLE_H
#define CODE_NUMSIM_FIELDVARIABLE_H

#include "Array2D.h"

class FieldVariable : public Array2D {
private:
    const std::array<double, 2> origin_;
    const std::array<double, 2> meshWidth_;
public:
    FieldVariable(std::array<int, 2> size, std::array<double, 2> origin,
                  std::array<double, 2> meshWidth);

    double interpolateAt(double x, double y) const;

    void setToZero();
};

#endif //CODE_NUMSIM_FIELDVARIABLE_H
