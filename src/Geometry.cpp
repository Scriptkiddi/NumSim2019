
#include <assert.h>
#include <iostream>
#include "Geometry.h"
#include "Util/GeometryParser.h"
#include <memory>

Geometry::Geometry(std::array<int, 2> size) : size_(size) {
    pressure_.resize(size[0] * size[1]);
    velocity_.resize(size[0] * size[1]);
    temperature_.resize(size[0] * size[1]);
    state_.resize(size[0] * size[1]);
    cout << size[0] << endl;
    cout << size[1] << endl;
}

std::pair<std::basic_string<char>, std::vector<double>> Geometry::pressure(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) pressure_.size());

    return pressure_[index];
}

std::pair<std::basic_string<char>, std::vector<double>> Geometry::temperature(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) temperature_.size());

    return temperature_[index];
}

std::pair<std::basic_string<char>, std::vector<double>> Geometry::velocity(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) velocity_.size());

    return velocity_[index];
}

std::pair<std::basic_string<char>, std::vector<double>> Geometry::state(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) state_.size());

    return state_[index];
}

//! get the size
std::array<int, 2> Geometry::size() const {
    return size_;
}

bool Geometry::isFluid(int i, int j) {
/*
    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    //assert(j * size_[0] + i < (int) data_.size());

    return isFluid_[nCellsX * j + i];
*/
    if(this->state(i,j).first == "S"){
        return false;
    }else{
        return true;
    }

}

int Geometry::nCellsFluid() {
    return nCellsFluid_;
}

void Geometry::countFluidCells() {
    for (int j = 0; j < nCellsY; j++) {
        for (int i = 0; i < nCellsX; i++) {
            if (isFluid(i, j)) {
                nCellsFluid_++;
            }
        }
    }
}