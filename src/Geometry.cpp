
#include "Geometry.h"

Geometry::Geometry(int nX, int nY): nCellsX(nX), nCellsY(nY){

    // allocate data, initialize to 0
    int nCellsFluid = nX * nY;

    isFluid_.resize(nX*nY);

    VeloTopBoundType.resize(nX, "");

    VeloBottomBoundType.resize(nX, "");

    VeloLeftBoundType.resize(nY, "");

    VeloRightBoundType.resize(nY, "");

    pTopBoundType.resize(nX, "");

    pBottomBoundType.resize(nX, "");

    pLeftBoundType.resize(nY, "");

    pRightBoundType.resize(nY, "");

    tTopBoundType.resize(nX, "");

    tBottomBoundType.resize(nX, "");

    tLeftBoundType.resize(nY, "");

    tRightBoundType.resize(nY, "");

    uTopBoundValue.resize(nX, 0.0);

    uBottomBoundValue.resize(nX, 0.0);

    uLeftBoundValue.resize(nY, 0.0);

    uRightBoundValue.resize(nY, 0.0);

    vTopBoundValue.resize(nX, 0.0);

    vBottomBoundValue.resize(nX, 0.0);

    vLeftBoundValue.resize(nY, 0.0);

    vRightBoundValue.resize(nY, 0.0);

    pTopBoundValue.resize(nX, 0.0);

    pBottomBoundValue.resize(nX, 0.0);

    pLeftBoundValue.resize(nY, 0.0);

    pRightBoundValue.resize(nY, 0.0);

    tTopBoundValue.resize(nX, 0.0);

    tBottomBoundValue.resize(nX, 0.0);

    tLeftBoundValue.resize(nY, 0.0);

    tRightBoundValue.resize(nY, 0.0);
}

bool Geometry::isFluid(int x, int y){
    return isFluid_[nCellsX * y + x];
}