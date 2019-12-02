
#ifndef CODE_NUMSIM_GEOMETRY_H
#define CODE_NUMSIM_GEOMETRY_H

#include<vector>
#include<array>
#include "Array2D/Array2D.h"

class Geometry{

public:
    //! constructor
    explicit Geometry(std::array<int,2> size);

    //! get the size
    std::array<int,2> size() const;
    bool isFluid(int x, int y);

    void countFluidCells();

    int nCellsX;
    int nCellsY;

    int nCellsFluid;

    std::vector< std::pair <std::string,std::vector<double>>> data_;  //< storage array values, in row-major order

    //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    std::pair<std::basic_string<char>, std::vector<double>> operator()(int i, int j);

    //! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
    std::pair<std::basic_string<char>, std::vector<double>> operator()(int i, int j) const;
    const std::array<int,2> size_ = {0,0};    //< width, height of the domain


private:
    std::vector<bool> isFluid_;
};


#endif
