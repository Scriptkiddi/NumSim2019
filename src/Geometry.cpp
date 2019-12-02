
#include <assert.h>
#include <iostream>
#include "Geometry.h"

Geometry::Geometry(std::array<int, 2> size): size_(size){
    data_.resize(size[0] * size[1]);
}


std::pair<std::basic_string<char>, std::vector<double>> Geometry::operator()(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) data_.size());

    return data_[index];
}

//! get the size
std::array<int, 2> Geometry::size() const {
    return size_;
}

