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
    assert(x >= 0 && x <= meshWidth_[0] * (size_[0]) && y >= 0 && y <= meshWidth_[1] * (size_[1]));
    int i = std::floor((x - origin_[0])/ meshWidth_[0] );
    int j = std::floor((y - origin_[1])/ meshWidth_[1] );
    if (i == size()[0] - 1) {
        i=i-1;
    }
    if (j == size()[1] - 1) {
        j=j-1;
    }

    double blX = i * meshWidth_[0] + origin_[0];
    double blY = j * meshWidth_[1] + origin_[1];

    double blDistanceX = x - blX;
    double blDistanceY = y - blY;

    double alphaX = blDistanceX / meshWidth_[0];
    double alphaY = blDistanceY / meshWidth_[1];

    return (1 - alphaX) * (1- alphaY) * operator()(i, j)
           + alphaX * (1- alphaY) * operator()(i + 1, j)
           + (1 - alphaX) * alphaY * operator()(i, j + 1)
           + alphaX * alphaY * operator()(i + 1, j + 1);

}

void FieldVariable::setToZero() {


}
