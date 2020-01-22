//
// Created by Julia Pelzer on 26.10.2019.
//

#include "FieldVariable.h"
#include <cmath>
#include <cassert>

const std::array<double, 2> meshWidth_{};

FieldVector::FieldVector(std::array<int, 3> size,
                             std::array<double, 2> meshWidth) :
        size_(size), meshWidth_(meshWidth) {
    // allocate data, initialize to 0
    data_.resize(size_[0]*size_[1]*size_[2], 0.0);
}

//! get the size
std::array<int,3> FieldVector::size() const
{
    return size_;
}

double &FieldVector::operator()(int i, int j, int k)
{
    const int index = j*size_[0]*size_[2] + i*size[2] + k;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(0 <= k && k < size_[2]);
    assert(j*size_[0]*size_[2] + i*size[2] + k < (int)data_.size());

    return data_[index];
}

double FieldVector::operator()(int i, int j, int k) const
{
    const int index = j*size_[0]*size_[2] + i*size[2] + k;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(0 <= k && k < size_[2]);
    assert(j*size_[0]*size_[2] + i*size[2] + k < (int)data_.size());

    return data_[index];
}

