//
// Created by Julia Pelzer on 26.10.2019.
//

#include "FieldVariable.h"
#include "math.h"


const std::array<double, 2> origin_{};
const std::array<double, 2> meshWidth_{};

FieldVariable::FieldVariable(std::array<int, 2> size, std::array<double, 2> origin,
                             std::array<double, 2> meshWidth) :
        Array2D(size), origin_(origin), meshWidth_(meshWidth) {}

double FieldVariable::interpolateAt(double x, double y) const {
    int i = floor((x - origin_[1]) / meshWidth_[1]);
    double x_rest = x - origin_[1] - i * meshWidth_[1];

    int j = floor((y - origin_[2]) / meshWidth_[2]);
    double y_rest = y - origin_[2] - j * meshWidth_[2];

    double alpha = x_rest / meshWidth_[1];
    double beta = y_rest / meshWidth_[2];

    double value = 0; // to be continued

    return value;
}

