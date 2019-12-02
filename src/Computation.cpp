//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Computation.h"
#include "Settings.h"
#include <cmath>
#include <iostream>
#include <StaggeredGrid/CentralDifferences.h>
#include <StaggeredGrid/DonorCell.h>
#include <PressureSolver/SOR.h>
#include <PressureSolver/GaussSeidel.h>

using namespace std;

void Computation::initialize(int argc, char **argv) {
    std::cout << "Running with" << argv[0] << std::endl;
    Settings settings;
    settings_ = settings;
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();
    array<int, 2> nCellsBoundary = {settings_.nCells[0] + 2, settings_.nCells[1] + 2}; // Mit Ghost cells

    Geometry geometry(nCellsBoundary[0], nCellsBoundary[1]);
    geometry_ = geometry;

    //initialize temperature boundary values //TEMPORARILY
    t_h = 1;
    t_c = 0;

    //initialize initial values
    uInit = settings_.uInit;
    vInit = settings_.vInit;
    pInit = settings_.pInit;
    tInit = settings_.tInit;

    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] /
                    (nCellsBoundary[0]);
    meshWidth_[1] = settings_.physicalSize[1] / (nCellsBoundary[1]);

    //initialize discretization
    if (!settings_.useDonorCell) {
        CentralDifferences grid(nCellsBoundary, meshWidth_, settings_.gamma);
        discretization_ = make_shared<CentralDifferences>(grid);
    } else {
        DonorCell grid(nCellsBoundary, meshWidth_, settings_.alpha, settings_.gamma);
        discretization_ = make_shared<DonorCell>(grid);
    }

    //initialize explicit pressureSolver
    if (settings_.pressureSolver == "SOR") {
        SOR pSolver(discretization_, std::make_shared<Geometry>(geometry_),settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
        pressureSolver_ = make_unique<SOR>(pSolver);
    } else {
        //GaussSeidel pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
        //pressureSolver_ = make_unique<PressureSolver>(pSolver);
        std::cout << "Please select SOR-solver" << std::endl;
    }
    //initialize outputWriters
    OutputWriterText outText(discretization_);
    outputWriterText_ = make_unique<OutputWriterText>(outText);

    OutputWriterParaview outPara(discretization_);
    outputWriterParaview_ = make_unique<OutputWriterParaview>(outPara);
}

void Computation::runSimulation() {
    double t = 0;
    applyInitialConditions();
    applyBoundaryValuesVelocities();
    applyBoundaryValuesTemperature();
    while (t < settings_.endTime) {
        computeTimeStepWidth();
        computeTemperature();
        applyBoundaryValuesTemperature();
        PreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        applyBoundaryValuesVelocities();
        t += dt_;
        outputWriterParaview_.get()->writeFile(t);
        outputWriterText_.get()->writeFile(t);
        cout << "current time: " << t << " dt: " << dt_ << " pressure solver iterations: " << endl;
    }
}

void Computation::computeTimeStepWidth() {
    // Take min time step stemming from speed considerations
    double uMaximum = discretization_.get()->u(discretization_.get()->uIBegin(), discretization_.get()->uJEnd());
    for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
        for (int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd() + 1; i++) {
            if (uMaximum < fabs(discretization_.get()->u(i, j))) {
                uMaximum = fabs(discretization_.get()->u(i, j));
            }
        }
    }
    double vMaximum = discretization_.get()->v(discretization_.get()->vIBegin(), discretization_.get()->vJEnd());
    for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
        for (int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++) {
            if (vMaximum < fabs(discretization_.get()->v(i, j))) {
                vMaximum = fabs(discretization_.get()->v(i, j));
            }
        }
    }
    double condition_diffusion = pow(discretization_.get()->dx(), 2) * pow(discretization_.get()->dy(), 2) /
                                 (pow(discretization_.get()->dx(), 2) + pow(discretization_.get()->dy(), 2)) *
                                 settings_.re /
                                 2;
    double condition_convection1 = discretization_.get()->dx() / uMaximum;
    double condition_convection2 = discretization_.get()->dy() / vMaximum;

    // Min time requirements for temperature
    double condition_temp = settings_.re * settings_.prandtl / 2 *
                            pow((1 / pow(discretization_.get()->dx(), 2) + 1 / pow(discretization_.get()->dy(), 2)),
                                -1);

    dt_ = min(condition_convection1, condition_convection2);
    dt_ = min(condition_diffusion, dt_);
    dt_ = min(condition_temp, dt_);
    dt_ = min(settings_.maximumDt, dt_) * settings_.tau;
}

void Computation::applyBoundaryValuesVelocities() {
    // u
    // TODO: add F
    //left bound
    int i_low = discretization_.get()->uIBegin() - 1;
    for(int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++){
        if(!geometry_.isFluid(i_low,j) && geometry_.isFluid(i_low + 1, j)){
            if(geometry_.VeloLeftBoundType[j] == "Neumann"){
                discretization_.get()->u(i_low,j) = discretization_.get()->u(i_low + 1, j) - discretization_.get()->dx() * geometry_.uLeftBoundValue[j];
                discretization_.get()->f(i_low,j) = discretization_.get()->u(i_low,j);
            }else{
                discretization_.get()->u(i_low,j) = geometry_.uLeftBoundValue[j];
                discretization_.get()->f(i_low,j) = discretization_.get()->u(i_low,j);
            }
        }
    }
    //right bound
    int i_high = discretization_.get()->uIEnd() + 1;
    for(int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++){
        if(!geometry_.isFluid(i_high,j) && geometry_.isFluid(i_high - 1, j)){
            if(geometry_.VeloRightBoundType[j] == "Neumann"){
                discretization_.get()->u(i_high,j) = discretization_.get()->u(i_high - 1, j) - discretization_.get()->dx() * geometry_.uRightBoundValue[j];
                discretization_.get()->f(i_high,j) = discretization_.get()->u(i_high,j);
            }else{
                discretization_.get()->u(i_high,j) = geometry_.uRightBoundValue[j];
                discretization_.get()->f(i_high,j) = discretization_.get()->u(i_high,j);
            }
        }
    }
    //bottom bound
    int j_low = discretization_.get()->uJBegin() - 1;
    for(int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd() + 1; i++){
        if(!geometry_.isFluid(i,j_low) && geometry_.isFluid(i,j_low + 1)){
            if(geometry_.VeloRightBoundType[i] == "Neumann"){
                double uOld = discretization_.get()->u(i,j_low);
                discretization_.get()->u(i,j_low) = discretization_.get()->u(i, j_low + 1) - discretization_.get()->dy() * geometry_.uBottomBoundValue[i];
                discretization_.get()->f(i,j_low) = 2 * discretization_.get()->u(i,j_low) - uOld;
            }else{
                discretization_.get()->u(i,j_low) = 2 * geometry_.uBottomBoundValue[i] - discretization_.get()->u(i,j_low + 1);
                discretization_.get()->u(i,j_low) = discretization_.get()->u(i,j_low);
            }
        }
    }
    //top bound
    int j_high = discretization_.get()->uJEnd() - 1;
    for(int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd() + 1; i++){
        if(!geometry_.isFluid(i,j_low) && geometry_.isFluid(i,j_low + 1)){
            if(geometry_.VeloRightBoundType[i] == "Neumann"){
                double uOld = discretization_.get()->u(i,j_high);
                discretization_.get()->u(i,j_high) = discretization_.get()->u(i, j_high + 1) - discretization_.get()->dy() * geometry_.uTopBoundValue[i];
                discretization_.get()->f(i,j_high) = 2 * discretization_.get()->u(i,j_high) - uOld;
            }else{
                discretization_.get()->u(i,j_high) = 2 * geometry_.uTopBoundValue[i] - discretization_.get()->u(i,j_high - 1);
                discretization_.get()->f(i,j_high) = discretization_.get()->u(i,j_high);
            }
        }
    }
    // inner cells
    for(int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++){
        for(int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); i++){
            //checking right neighbour
            if(!geometry_.isFluid(i,j) && !geometry_.isFluid(i + 1,j)){
                //check top and bottom neighbours
                if(geometry_.isFluid(i,j + 1)){
                    discretization_.get()->u(i,j) = - discretization_.get()->u(i, j + 1);
                    discretization_.get()->f(i,j) = discretization_.get()->u(i,j);
                }else if(geometry_.isFluid(i, j - 1)){
                    discretization_.get()->u(i,j) = - discretization_.get()->u(i, j - 1);
                    discretization_.get()->f(i,j) = discretization_.get()->u(i,j);
                }
            }
        }
    }

    // v
    // TODO: add G
    //left bound
    int i_low = discretization_.get()->vIBegin() - 1;
    for(int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++){
        if(!geometry_.isFluid(i_low,j) && geometry_.isFluid(i_low + 1, j)){
            if(geometry_.VeloLeftBoundType[j] == "Neumann"){
                double vOld;
                discretization_.get()->v(i_low,j) = discretization_.get()->v(i_low + 1, j) - discretization_.get()->dx() * geometry_.vLeftBoundValue[j];
                discretization_.get()->g(i_low,j) = 2 * discretization_.get()->v(i_low,j) - vOld;
            }else{
                discretization_.get()->v(i_low,j) = 2 * geometry_.vLeftBoundValue[j] - discretization_.get()->v(i_low + 1,j);
                discretization_.get()->g(i_low,j) = discretization_.get()->v(i_low,j);
            }
        }
    }
    //right bound
    int i_high = discretization_.get()->vIEnd() + 1;
    for(int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++){
        if(!geometry_.isFluid(i_high,j) && geometry_.isFluid(i_high - 1, j)){
            if(geometry_.VeloRightBoundType[j] == "Neumann"){
                double vOld = discretization_.get()->v(i_high,j);
                discretization_.get()->v(i_high,j) = discretization_.get()->v(i_high - 1, j) - discretization_.get()->dx() * geometry_.vRightBoundValue[j];
                discretization_.get()->g(i_high,j) = 2 * discretization_.get()->v(i_high,j) - vOld;
            }else{
                discretization_.get()->v(i_high,j) = 2 * geometry_.vRightBoundValue[j] - discretization_.get()->v(i_high - 1,j);
                discretization_.get()->g(i_high,j) = discretization_.get()->v(i_high,j);
            }
        }
    }
    //bottom bound
    int j_low = discretization_.get()->vJBegin() - 1;
    for(int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++){
        if(!geometry_.isFluid(i,j_low) && geometry_.isFluid(i,j_low + 1)){
            if(geometry_.VeloRightBoundType[i] == "Neumann"){
                discretization_.get()->v(i,j_low) = discretization_.get()->v(i, j_low + 1) - discretization_.get()->dy() * geometry_.vBottomBoundValue[i];
                discretization_.get()->g(i,j_low) = discretization_.get()->v(i,j_low);
            }else{
                discretization_.get()->v(i,j_low + 1) = geometry_.vBottomBoundValue[i];
                discretization_.get()->g(i,j_low) = discretization_.get()->v(i,j_low);
            }
        }
    }
    //top bound
    int j_high = discretization_.get()->vJEnd() - 1;
    for(int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++){
        if(!geometry_.isFluid(i,j_high) && geometry_.isFluid(i,j_high + 1)){
            if(geometry_.VeloRightBoundType[i] == "Neumann"){
                discretization_.get()->v(i,j_high) = discretization_.get()->v(i, j_high + 1) - discretization_.get()->dy() * geometry_.vTopBoundValue[i];
                discretization_.get()->g(i,j_high) = discretization_.get()->v(i,j_high);
            }else{
                discretization_.get()->v(i,j_high) = geometry_.vBottomBoundValue[i];
                discretization_.get()->g(i,j_high) = discretization_.get()->v(i,j_high);
            }
        }
    }
    //inner cells
    for(int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++){
        for(int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); i++){
            if(!geometry_.isFluid(i,j) && !geometry_.isFluid(i, j + 1)){
                if(geometry_.isFluid(i + 1, j)){
                    discretization_.get()->v(i,j) = - discretization_.get()->v(i + 1, j);
                    discretization_.get()->g(i,j) = discretization_.get()->v(i,j);
                }else if(geometry_.isFluid(i - 1, j)){
                    discretization_.get()->v(i,j) = - discretization_.get()->v(i - 1, j);
                    discretization_.get()->g(i,j) = discretization_.get()->v(i,j);
                }
            }
        }
    }
}

void Computation::applyBoundaryValuesTemperature(){
    //outer bounds
    //left bound without corners
    int i_low = discretization_.get()->tIBegin()-1;
    for (int j = discretization_.get()->tJBegin() ; j <= discretization_.get()->tJEnd() ; j++) {
        if(geometry_.tLeftBoundType[j] == "Neumann"){
            discretization_.get()->t(i_low, j) = discretization_.get()->t(i_low + 1, j) - discretization_.get()->dx() * geometry_.tLeftBoundValue[j];
        }else{
            discretization_.get()->t(i_low, j) = 2 * geometry_.tLeftBoundValue[j] - discretization_.get()->t(i_low + 1, j); 
        }
    }

    //right bound without corners
    int i_high = discretization_.get()->tIEnd()+1;
    for (int j = discretization_.get()->tJBegin() ; j <= discretization_.get()->tJEnd() ; j++) {
        if(geometry_.tRightBoundType[j] == "Neumann"){
            discretization_.get()->t(i_high, j) = discretization_.get()->t(i_high - 1, j) - discretization_.get()->dx() * geometry_.tRightBoundValue[j];
        }else{
            discretization_.get()->t(i_high, j) = 2 * geometry_.tRightBoundValue[j] - discretization_.get()->t(i_high + 1, j); 
        }
    }

    //bottom bound with corners
    int j_low = discretization_.get()->tJBegin()-1;
    for (int i = discretization_.get()->tIBegin()-1; i <= discretization_.get()->tIEnd()+1; i++) {
        if(geometry_.tBottomBoundType[i] == "Neumann"){
            discretization_.get()->t(i, j_low) = discretization_.get()->t(i, j_low + 1) - discretization_.get()->dy() * geometry_.tBottomBoundValue[i];
        }else{
            discretization_.get()->t(i, j_low) = 2 * geometry_.tBottomBoundValue[i] - discretization_.get()->t(i, j_low + 1);
        }
    }

    //top bound with corners
    int j_high = discretization_.get()->tJEnd()-1;
    for (int i = discretization_.get()->tIBegin()-1; i <= discretization_.get()->tIEnd()+1; i++) {
        if(geometry_.tTopBoundType[i] == "Neumann"){
            discretization_.get()->t(i, j_high) = discretization_.get()->t(i, j_high - 1) - discretization_.get()->dy() * geometry_.tTopBoundValue[i];
        }else{
            discretization_.get()->t(i, j_high) = 2 * geometry_.tTopBoundValue[i] - discretization_.get()->t(i, j_high - 1);
        }
    }

    //inner cells (obstacles)
    for(int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++){
        for(int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++){
            if(!geometry_.isFluid(i,j)){
                int N = 0;
                double tTmp = 0;

                //check left and right neighbours
                if(geometry_.isFluid(i - 1, j)){
                    tTmp += discretization_.get()->t(i - 1, j);
                    N++;
                }else if(geometry_.isFluid(i + 1, j)){
                    tTmp += discretization_.get()->t(i + 1, j);
                    N++;
                }
                //check top and bottom neighbours
                if(geometry_.isFluid(i, j - 1)){
                    tTmp += discretization_.get()->t(i, j - 1);
                    N++;
                }else if(geometry_.isFluid(i, j + 1)){
                    tTmp += discretization_.get()->t(i, j + 1);
                    N++;
                }
                //assign average value
                if(N > 0){
                    discretization_.get()->t(i,j) = tTmp / N;
                }
            }
        }
    }
}

void Computation::PreliminaryVelocities() {
    for (int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++) {
        for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
            if(geometry_.isFluid(i,j) && geometry_.isFluid(i + 1, j)){
                discretization_.get()->f(i, j) =
                        (discretization_.get()->u(i, j) +
                        dt_ * (1 / settings_.re * (discretization_.get()->computeD2uDx2(i, j) +
                                                    discretization_.get()->computeD2uDy2(i, j)) -
                                discretization_.get()->computeDu2Dx(i, j) -
                                discretization_.get()->computeDuvDy(i, j) + settings_.g[0]))
                        - dt_ * settings_.beta * settings_.g[0] *
                        (discretization_.get()->t(i, j) + discretization_.get()->t(i + 1, j)) / 2;
            }
        }
    }

    for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
        for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
            if(geometry_.isFluid(i,j) && geometry_.isFluid(i, j + 1)){
                discretization_.get()->g(i, j) =
                        (discretization_.get()->v(i, j) +
                        dt_ * (1 / settings_.re * (discretization_.get()->computeD2vDy2(i, j) +
                                                    discretization_.get()->computeD2vDx2(i, j)) -
                                discretization_.get()->computeDv2Dy(i, j) -
                                discretization_.get()->computeDuvDx(i, j) + settings_.g[1]))
                        - dt_ * settings_.beta * settings_.g[1] *
                        (discretization_.get()->t(i, j) + discretization_.get()->t(i, j + 1)) / 2;
            }
        }
    }
}

void Computation::computeRightHandSide() {
    for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
        for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
            if(geometry_.isFluid(i,j)){
                discretization_.get()->rhs(i, j) =
                        1 / dt_ * (((discretization_.get()->f(i, j) - discretization_.get()->f(i - 1, j)) /
                                    discretization_.get()->dx()) +
                                ((discretization_.get()->g(i, j) - discretization_.get()->g(i, j - 1)) /
                                    discretization_.get()->dy()));
            }
        }
    }
}

void Computation::computePressure() {

    pressureSolver_->solve();
}

void Computation::computeVelocities() {
    for (int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++) {
        for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
            if(geometry_.isFluid(i,j) && geometry_.isFluid(i + 1, j)){
                discretization_.get()->u(i, j) =
                        discretization_.get()->f(i, j) - dt_ * discretization_.get()->computeDpDx(i, j);
            }
        }
    }

    for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            if(geometry_.isFluid(i,j) && geometry_.isFluid(i, j + 1)){
                discretization_.get()->v(i, j) =
                        discretization_.get()->g(i, j) - dt_ * discretization_.get()->computeDpDy(i, j);
            }
        }
    }
}

void Computation::computeTemperature() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if(geometry_.isFluid(i,j)){
                discretization_.get()->t(i, j) = discretization_.get()->t(i, j) +
                                                dt_ * (1 / settings_.re * 1 / settings_.prandtl * (
                                                        discretization_.get()->computeD2TDx2(i, j)
                                                        +
                                                        discretization_.get()->computeD2TDy2(i, j)
                                                )
                                                        - discretization_.get()->computeDuTDx(i, j)
                                                        - discretization_.get()->computeDvTDy(i, j)
                                                );
            }
        }
    }
}

void Computation::applyInitialConditions(){
    for (int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++) {
        for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
            if(geometry_.isFluid(i,j) && geometry_.isFluid(i + 1, j)){
                discretization_.get()->u(i, j) = uInit;
            }else{
                discretization_.get()->u(i, j) = std::nan(""); 
            }
        }
    }

    for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            if(geometry_.isFluid(i,j) && geometry_.isFluid(i, j + 1)){
                discretization_.get()->v(i, j) = vInit;
            }else{
                discretization_.get()->v(i,j) = std::nan("");
            }
        }
    }

    for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
        for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
            if(geometry_.isFluid(i,j)){
                discretization_.get()->p(i, j) = pInit;
            }else{
                discretization_.get()->p(i,j) = std::nan("");
            }
        }
    }

    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if(geometry_.isFluid(i,j)){
                discretization_.get()->t(i, j) = tInit;
            }else{
                discretization_.get()->t(i,j) = std::nan("");
            }
        }
    }
}