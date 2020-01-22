//
// Created by Julia Pelzer on 26.10.2019.
//

#ifndef CODE_NUMSIM_FIELDVARIABLE_H
#define CODE_NUMSIM_FIELDVARIABLE_H

#include "Array2D.h"

class FieldVector{
private:
    std::vector<double> data_;
    std::array<int,3> size_;
    const std::array<double, 2> meshWidth_;
public:
    explicit FieldVector(std::array<int,3> size);
    std::array<int,3> size() const;
    double &operator()(int i, int j, int k); //i, j = cell indizes, k direction of velocity (0,...,8)
    double operator()(int i, int j, int k) const;

    FieldVector(std::array<int, 3> size, std::array<double, 2> origin,
                  std::array<double, 2> meshWidth);
};

#endif //CODE_NUMSIM_FIELDVARIABLE_H
