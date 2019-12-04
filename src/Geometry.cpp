
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
}

std::pair<std::string, std::vector<double>> Geometry::get_pressure(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) pressure_.size());

    return pressure_[index];
}
void Geometry::set_pressure(int i, int j, std::pair<std::string, std::vector<double>> value) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) pressure_.size());

    pressure_[index] = value;
}


std::pair<std::string, std::vector<double>> Geometry::get_velocity(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) velocity_.size());

    return velocity_[index];
}
void Geometry::set_velocity(int i, int j, std::pair<std::string, std::vector<double>> value) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) velocity_.size());
    cout << value.second.size() << endl;
    cout << i << "x" << j << endl;

    velocity_[index] = value;
    cout << velocity_[index].second.size() << endl;
    cout << "--" << endl;
}
std::pair<std::string, std::vector<double>> Geometry::get_temperature(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) temperature_.size());

    return temperature_[index];
}
void Geometry::set_temperature(int i, int j, std::pair<std::string, std::vector<double>> value) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) temperature_.size());

    temperature_[index] = value;
}
std::pair<std::string, std::vector<double>> Geometry::get_state(int i, int j) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) state_.size());

    return state_[index];
}
void Geometry::set_state(int i, int j, std::pair<std::string, std::vector<double>> value) {
    const int index = j * size_[0] + i;

    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) state_.size());

    state_[index] = value;
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
    if(this->get_state(i,j).first == "F"){
        return true;
    }else{
        return false;
    }

}

int Geometry::nCellsFluid() {
    return nCellsFluid_;
}

void Geometry::countFluidCells() {

    for (int j = 0; j < size_[1]; j++) {
        for (int i = 0; i < size_[0]; i++) {
            if (isFluid(i, j)) {
                nCellsFluid_++;
            }
        }
    }
}