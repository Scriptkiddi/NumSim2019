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
    // calculate corspeonding cell (bottom left point)
    // cartesian (0,0) -> index (0,0) etc.
    int i = std::floor((x - origin_[0])/ meshWidth_[0] );
    int j = std::floor((y - origin_[1])/ meshWidth_[1] );

    // the right and upper boundaries have to be corrected, as those would use cells outside our grid
    if (i == size()[0] - 1) {
        i=i-1;
    }

    if (j == size()[1] - 1) {
        j=j-1;
    }

    // calculate cartesian coordinates of the bottom left point
    double blX = i * meshWidth_[0] + origin_[0];
    double blY = j * meshWidth_[1] + origin_[1];

    // calculate distances for bilinear interpolation
    // distance to bottom left point
    double blDistanceX = x - blX;
    double blDistanceY = y - blY;

    // calculate coefficient for bilinear interpolation
    double alphaX = blDistanceX / meshWidth_[0];
    double alphaY = blDistanceY / meshWidth_[1];

    // calculate bilinear interpolation for cartesian (x,y)
    return (1 - alphaX) * (1- alphaY) * (*this)(i, j)
           + alphaX * (1- alphaY) * (*this)(i + 1, j)
           + (1 - alphaX) * alphaY * (*this)(i, j + 1)
           + alphaX * alphaY * (*this)(i + 1, j + 1);

    //assert(x >= 0 && x <= meshWidth_[0] * (size_[0]) && y >= 0 && y <= meshWidth_[1] * (size_[1]));

    //int i = floor((x - origin_[0]) / meshWidth_[0]);
    //int j = floor((y - origin_[1]) / meshWidth_[1]);
    //if(i == size()[0]){
    //    i = i -1;
    //}
    //if(j == size()[0]){
    //    j = j -1;
    //}
    //double alpha = (x - origin_[0] - i * meshWidth_[0]) / meshWidth_[0];
    //double beta = (y - origin_[1] - j * meshWidth_[1]) / meshWidth_[1];
    //const int index = j * size_[0] + i;

    //return ((1 - beta) * ((1 - alpha) * operator()(i,j) + alpha *  operator()(i+1,j)) +
    //        beta * ((1 - alpha) * operator()(i,j+1)+ alpha * operator()(i+1,j+1)));

}