
#ifndef CODE_NUMSIM_GEOMETRY_H
#define CODE_NUMSIM_GEOMETRY_H

#include<vector>
#include<array>
#include "Array2D/Array2D.h"

struct Geometry{

    Geometry(int nX, int nY);

    int nCellsX;

    int nCellsY;

    bool isFluid(int x, int y);
    
    int nCellsFluid;

    std::vector<std::string> VeloTopBoundType;

    std::vector<std::string> VeloBottomBoundType;

    std::vector<std::string> VeloLeftBoundType;

    std::vector<std::string> VeloRightBoundType;

    std::vector<std::string> pTopBoundType;

    std::vector<std::string> pBottomBoundType;

    std::vector<std::string> pLeftBoundType;

    std::vector<std::string> pRightBoundType;

    std::vector<std::string> tTopBoundType;

    std::vector<std::string> tBottomBoundType;

    std::vector<std::string> tLeftBoundType;

    std::vector<std::string> tRightBoundType;

    std::vector<double> uTopBoundValue;

    std::vector<double> uBottomBoundValue;

    std::vector<double> uLeftBoundValue;

    std::vector<double> uRightBoundValue;

    std::vector<double> vTopBoundValue;

    std::vector<double> vBottomBoundValue;

    std::vector<double> vLeftBoundValue;

    std::vector<double> vRightBoundValue;

    std::vector<double> pTopBoundValue;

    std::vector<double> pBottomBoundValue;

    std::vector<double> pLeftBoundValue;

    std::vector<double> pRightBoundValue;

    std::vector<double> tTopBoundValue;

    std::vector<double> tBottomBoundValue;

    std::vector<double> tLeftBoundValue;

    std::vector<double> tRightBoundValue;

    private:
        std::vector<bool> isFluid_;
};


#endif
