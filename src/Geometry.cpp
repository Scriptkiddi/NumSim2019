
#include <assert.h>
#include <iostream>
#include "Geometry.h"
#include "Util/GeometryParser.h"
#include <memory>

Geometry::Geometry(std::array<int, 2> size): size_(size){
    data_.resize(size[0] * size[1]);
}

//TODO delete? replace with velocity, temperature, pressure
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

bool Geometry::isFluid(int i, int j){
/*
    // assert that indices are in range
    assert(0 <= i && i < size_[0]);
    assert(0 <= j && j < size_[1]);
    assert(j * size_[0] + i < (int) data_.size());

    return isFluid_[nCellsX * j + i];
*/
    if(this.state(i,j)=="S"){
        return false;
    }else{
        return true;
    }

}

int Geometry::nCellsFluid(){
    return nCellsFluid_;
}

void Geometry::countFluidCells(){
    for(int j = 0; j <= nCellsY; j++){
        for(int i = 0; i <= nCellsX; i++){
            if(isFluid(i,j)){
                nCellsFluid_ ++;
            }
        }
    }
}