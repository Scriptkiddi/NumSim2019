
#ifndef CODE_NUMSIM_GEOMETRY_H
#define CODE_NUMSIM_GEOMETRY_H

#include <vector>
#include <array>
#include "Array2D/Array2D.h"
#include <memory>

class Geometry{

public:
    //! constructor
    Geometry(std::array<int,2> size);

    //! get the size
    std::array<int,2> size() const;
    bool isFluid(int x, int y);

    void countFluidCells();

    int nCellsFluid();    

    //attributes 

    std::vector< std::pair <std::string,std::vector<double>>> temperature_;  //< storage array values, in row-major order
    std::vector< std::pair <std::string,std::vector<double>>> velocity_;  //< storage array values, in row-major order
    std::vector< std::pair <std::string,std::vector<double>>> pressure_;  //< storage array values, in row-major order
    std::vector< std::pair <std::string,std::vector<double>>> state_;  //< storage array values, in row-major order

    //! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
    std::pair<std::string, std::vector<double>> get_velocity(int i, int j); // "PR", ""
    void set_velocity(int i, int j, std::pair<std::string, std::vector<double>> value); // "PR", ""
    std::pair<std::string, std::vector<double>> get_temperature(int i, int j); // "PR", ""
    void set_temperature(int i, int j, std::pair<std::string, std::vector<double>> value); // "PR", ""
    std::pair<std::string, std::vector<double>> get_pressure(int i, int j); // "PR", ""
    void set_pressure(int i, int j, std::pair<std::string, std::vector<double>> value); // "PR", ""
    std::pair<std::string, std::vector<double>> get_state(int i, int j); // // "F" (fluid), "S" (solid)
    void set_state(int i, int j, std::pair<std::string, std::vector<double>> value); //"F" (fluid), "S" (solid)

    const std::array<int,2> size_ = {0,0};    //< width, height of the domain
    
private:
    std::vector<bool> isFluid_;

    int nCellsX;
    int nCellsY;

    int nCellsFluid_;

};


#endif
