//
// Created by Julia Pelzer on 26.10.2019.
//

#include "FieldVariable.h"
#include <cmath>
#include <cassert>


const std::array<double, 2> origin_{};
const std::array<double, 2> meshWidth_{};

FieldVariable::FieldVariable(std::array<int, 2> size, std::array<double, 2> origin,
                             std::array<double, 2> meshWidth) :
        Array2D(size), origin_(origin), meshWidth_(meshWidth) {}

double FieldVariable::interpolateAt(double x, double y) const {

    assert(x>=0 && x<= meshWidth_[0]*(size_[0]-1) && y>=0 && y<=meshWidth_[1]*(size_[1]-1));

    int i = floor((x - origin_[0]) / meshWidth_[0]);
    double alpha = (x - origin_[0] - i * meshWidth_[0]) / meshWidth_[0];

    int j = floor((y - origin_[1]) / meshWidth_[1]);
    double beta = (y - origin_[1] - j * meshWidth_[1]) / meshWidth_[1];

    const int index = j * size_[0] + i;


    return ((1 - beta) * ((1 - alpha) * data_[index] + alpha * data_[index + 1]) +
            beta * ((1 - alpha) * data_[index + size_[0]] + alpha * data_[index + size_[0] + 1]));

}