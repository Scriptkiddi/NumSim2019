//
// Created by Julia Pelzer on 26.10.2019.
//

#include <memory>
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

    geometry_ = settings_.geometry;

    //initialize initial values
    uInit = settings_.uInit;
    vInit = settings_.vInit;
    pInit = settings_.pInit;
    tInit = settings_.tInit;

    //initialize meshWidth
    meshWidth_[0] = settings_.physicalSize[0] / (settings_.nCells[0]);
    meshWidth_[1] = settings_.physicalSize[1] / (settings_.nCells[1]);

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
        SOR pSolver(discretization_, geometry_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
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
    double counter = 0;
    while (t < settings_.endTime) {
        computeTimeStepWidth();
        computeTemperature();
        applyBoundaryValuesTemperature();
        PreliminaryVelocities();
        computeTemperature();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        applyBoundaryValuesVelocities();
        t += dt_;
        counter += dt_;
        if (counter >= settings_.outputFileEveryDt){
            counter = 0;
            cout << "Writing paraview timestep: "<< t << endl;
            outputWriterParaview_.get()->writeFile(t);
            cout << "current time: " << t << " dt: " << dt_ << " pressure solver iterations: " << endl;h
            //outputWriterText_.get()->writeFile(t);
        }

    }

}

void Computation::computeTimeStepWidth() {
    double uMaximum = 0;
    for (int j = discretization_.get()->uJBegin()-1; j <= discretization_.get()->uJEnd()+1; j++) {
        for (int i = discretization_.get()->uIBegin()-1; i <= discretization_.get()->uIEnd()+1; i++) {
            if (uMaximum < abs(discretization_.get()->u(i, j))) {
                uMaximum = abs(discretization_.get()->u(i, j));
            }
        }
    }
    double vMaximum = 0;
    for (int j = discretization_.get()->vJBegin()-1; j <= discretization_.get()->vJEnd()+1; j++) {
        for (int i = discretization_.get()->vIBegin()-1; i <= discretization_.get()->vIEnd()+1; i++) {
            if (vMaximum < abs(discretization_.get()->v(i, j))) {
                vMaximum = abs(discretization_.get()->v(i, j));
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

//TODO wann Corners definieren!!! (siehe Skript)
//TODO laufen die Schleifen lang genug/ kurz genug/ richtig?

void Computation::applyBoundaryValuesVelocities() {
// type  options : "Neumann" (Abl = 0) oder "Dirichlet" (Wert = 0)
// type: input options : "SLW" (slip wall), "NSW" (no slip wall), "IN" (in wall), "OUT" (out wall), "PR" (pressure),
// "TD"( temperature Dirichlet), "TN" (temperature Neumann), "F" (fluid cell), "S" (solid cell)
// NoSlip ist komplett dirichlet, slip ist dirichlet in normal- und neumann in tangentialrichtung

    // help variables
    int i_low;
    int i_high;
    int j_low;
    int j_high;

    // u
    //left bound
    i_low = discretization_.get()->uIBegin() - 1;
    for(int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++){
        if(!geometry_.get()->isFluid(i_low,j) && geometry_.get()->isFluid(i_low + 1, j)){
            if(geometry_.get()->get_velocity(i_low,j).first == "NSW"){
                discretization_.get()->u(i_low,j) = 0;
                discretization_.get()->f(i_low,j) = discretization_.get()->u(i_low,j);
            }else if(geometry_.get()->get_velocity(i_low,j).first == "SLW"){
                discretization_.get()->u(i_low,j) = 0;
                discretization_.get()->f(i_low,j) = discretization_.get()->u(i_low,j);
            }else if(geometry_.get()->get_velocity(i_low,j).first == "IN"){
                discretization_.get()->u(i_low,j) = geometry_.get()->get_velocity(i_low,j).second[0];
                discretization_.get()->f(i_low,j) = discretization_.get()->u(i_low,j);
            }else if(geometry_.get()->get_velocity(i_low,j).first == "OUT"){
                //richtige Stelle für uOld und Neumann etc?
                //passt das so mit u(n), u(n+1) statt u(n-1), u(n)??? Indexfehler im Skript oder liegen wir falsch?
                double uOld = discretization_.get()->u(i_low,j);
                discretization_.get()->u(i_low,j) = discretization_.get()->u(i_low+1,j);
                discretization_.get()->f(i_low,j) = 2 * discretization_.get()->u(i_low,j) - uOld;
            }
        }
    }
    //right bound
    i_high = discretization_.get()->uIEnd() + 1;
    for(int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++){
        if(!geometry_.get()->isFluid(i_high,j) && geometry_.get()->isFluid(i_high - 1, j)){
            if(geometry_.get()->get_velocity(i_high,j).first == "NSW"){
                discretization_.get()->u(i_high,j) = 0;
                discretization_.get()->f(i_high,j) = discretization_.get()->u(i_high,j);
            }else if(geometry_.get()->get_velocity(i_high,j).first == "SLW"){
                discretization_.get()->u(i_high,j) = 0;
                discretization_.get()->f(i_high,j) = discretization_.get()->u(i_high,j);
            }else if(geometry_.get()->get_velocity(i_high,j).first == "IN"){
                discretization_.get()->u(i_high,j) = geometry_.get()->get_velocity(i_high,j).second[0];
                discretization_.get()->f(i_high,j) = discretization_.get()->u(i_high,j);
            }else if(geometry_.get()->get_velocity(i_high,j).first == "OUT"){
                double uOld = discretization_.get()->u(i_high,j);
                discretization_.get()->u(i_high,j) = discretization_.get()->u(i_high - 1,j);
                discretization_.get()->f(i_high,j) = 2 * discretization_.get()->u(i_high,j) - uOld;
            }
        }
    }
    //bottom bound
    j_low = discretization_.get()->uJBegin() - 1; //Schleife kürzer
    for(int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++){
        if(!geometry_.get()->isFluid(i,j_low) && geometry_.get()->isFluid(i,j_low + 1)){
            if(geometry_.get()->get_velocity(i,j_low).first == "NSW"){
                discretization_.get()->u(i,j_low) = - discretization_.get()->u(i,j_low + 1); //should be -u(1,j)
                discretization_.get()->f(i,j_low) = discretization_.get()->u(i,j_low);
            }else if(geometry_.get()->get_velocity(i,j_low).first == "SLW"){
                double uOld = discretization_.get()->u(i,j_low);
                discretization_.get()->u(i,j_low) = discretization_.get()->u(i,j_low + 1); //u(1,j)
                discretization_.get()->f(i,j_low) = 2 * discretization_.get()->u(i,j_low) - uOld;
            }else if(geometry_.get()->get_velocity(i,j_low).first == "IN"){
                discretization_.get()->u(i,j_low) = 2* geometry_.get()->get_velocity(i,j_low).second[0] - discretization_.get()->u(i,j_low+1);
                discretization_.get()->f(i,j_low) = discretization_.get()->u(i,j_low);
            }else if(geometry_.get()->get_velocity(i,j_low).first == "OUT"){
                double uOld = discretization_.get()->u(i,j_low);
                discretization_.get()->u(i,j_low) = discretization_.get()->u(i,j_low + 1);
                discretization_.get()->f(i,j_low) = 2 * discretization_.get()->u(i,j_low)- uOld;
            }
        }
    }
    //top bound
    j_high = discretization_.get()->uJEnd() + 1; //Schleife kürzer
    for(int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++){
        if(!geometry_.get()->isFluid(i,j_high) && geometry_.get()->isFluid(i,j_high - 1)){
            if(geometry_.get()->get_velocity(i,j_high).first == "NSW"){
                discretization_.get()->u(i,j_high) = - discretization_.get()->u(i,j_high - 1); //should be -u(1,j)
                discretization_.get()->f(i,j_high) = discretization_.get()->u(i,j_high);
            }else if(geometry_.get()->get_velocity(i,j_high).first == "SLW"){
                double uOld = discretization_.get()->u(i,j_high);
                discretization_.get()->u(i,j_high) = discretization_.get()->u(i,j_high - 1); //u(1,j)
                discretization_.get()->f(i,j_high) = 2 * discretization_.get()->u(i,j_high) - uOld;
            }else if(geometry_.get()->get_velocity(i,j_high).first == "IN"){
                discretization_.get()->u(i,j_high) = 2* geometry_.get()->get_velocity(i,j_high).second[0] - discretization_.get()->u(i,j_high-1);
                discretization_.get()->f(i,j_high) = discretization_.get()->u(i,j_high);
            }else if(geometry_.get()->get_velocity(i,j_high).first == "OUT"){
                double uOld = discretization_.get()->u(i,j_high);
                discretization_.get()->u(i,j_high) = discretization_.get()->u(i,j_high - 1);
                discretization_.get()->f(i,j_high) = 2 * discretization_.get()->u(i,j_high) - uOld;
            }
        }
    }
    // inner cells
    for(int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++){
        for(int j = discretization_.get()->uJBegin(); j <= discretization_.get()->uJEnd(); j++){
            if(!geometry_.get()->isFluid(i,j)){
                if(!geometry_.get()->isFluid(i + 1,j)){
                    discretization_.get()->u(i,j) = 0;
                    discretization_.get()->f(i,j) = 0;
                }else if(geometry_.get()->isFluid(i,j + 1)){
                    discretization_.get()->u(i,j) = - discretization_.get()->u(i, j + 1);
                    discretization_.get()->f(i,j) = discretization_.get()->u(i,j);
                }else if(geometry_.get()->isFluid(i, j - 1)){
                    discretization_.get()->u(i,j) = - discretization_.get()->u(i, j - 1);
                    discretization_.get()->f(i,j) = discretization_.get()->u(i,j);
                }
            }else if(!geometry_.get()->isFluid(i + 1, j)){
                discretization_.get()->u(i,j) = 0;
                discretization_.get()->f(i,j) = 0;
            }
        }
    }



    // v
    //left bound
    i_low = discretization_.get()->vIBegin() - 1; //Schleife kürzer
    for(int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++){
        if(!geometry_.get()->isFluid(i_low, j) && geometry_.get()->isFluid(i_low + 1, j)){
            if(geometry_.get()->get_velocity(i_low,j).first == "NSW"){
                discretization_.get()->v(i_low,j) = - discretization_.get()->v(i_low + 1,j); //should be -v(1,j)
                discretization_.get()->g(i_low,j) = discretization_.get()->v(i_low,j);
            }else if(geometry_.get()->get_velocity(i_low,j).first == "SLW"){
                double vOld = discretization_.get()->v(i_low,j);
                discretization_.get()->v(i_low,j) = discretization_.get()->v(i_low + 1,j); //v(1,j)
                discretization_.get()->g(i_low,j) = 2 * discretization_.get()->v(i_low,j) - vOld;
            }else if(geometry_.get()->get_velocity(i_low,j).first == "IN"){
                discretization_.get()->v(i_low,j) = 2* geometry_.get()->get_velocity(i_low,j).second[1] - discretization_.get()->v(i_low+1,j); //2v_in(0,j*h)-v(1,j)
                discretization_.get()->g(i_low,j) = discretization_.get()->v(i_low,j);
            }else if(geometry_.get()->get_velocity(i_low,j).first == "OUT"){
                double vOld = discretization_.get()->v(i_low,j);
                discretization_.get()->v(i_low,j) = discretization_.get()->v(i_low + 1,j);
                discretization_.get()->g(i_low,j) = 2 * discretization_.get()->v(i_low,j) - vOld;
            }
        }
    }
    //right bound
    i_high = discretization_.get()->vIEnd() + 1;//Schleife kürzer
    for(int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++){
        if(!geometry_.get()->isFluid(i_high,j) && geometry_.get()->isFluid(i_high - 1, j)){
            /*if(geometry_.get()->vRightBoundType[j] == "Neumann"){
                double vOld = discretization_.get()->v(i_high,j);
                discretization_.get()->v(i_high,j) = discretization_.get()->v(i_high - 1, j) - discretization_.get()->dx() * geometry_.get()->vRightBoundValue[j];
                discretization_.get()->g(i_high,j) = 2 * discretization_.get()->v(i_high,j) - vOld;
            */
            if(geometry_.get()->get_velocity(i_high,j).first == "NSW"){
                discretization_.get()->v(i_high,j) = - discretization_.get()->v(i_high - 1,j); //should be -v(end-1,j)
                discretization_.get()->g(i_high,j) = discretization_.get()->v(i_high,j);
            }else if(geometry_.get()->get_velocity(i_high,j).first == "SLW"){
                double vOld = discretization_.get()->v(i_high,j);
                discretization_.get()->v(i_high,j) = discretization_.get()->v(i_high - 1,j); //v(end-1,j)
                discretization_.get()->g(i_high,j) = 2 * discretization_.get()->v(i_high,j) - vOld;
            }else if(geometry_.get()->get_velocity(i_high,j).first == "IN"){
                discretization_.get()->v(i_high,j) = 2* geometry_.get()->get_velocity(i_high,j).second[1] - discretization_.get()->v(i_high-1,j); //2v_in(0,j*h)-v(1,j)
                discretization_.get()->g(i_high,j) = discretization_.get()->v(i_high,j);
            }else if(geometry_.get()->get_velocity(i_high,j).first == "OUT"){
                double vOld = discretization_.get()->v(i_high,j);
                discretization_.get()->v(i_high,j) = discretization_.get()->v(i_high - 1,j);
                discretization_.get()->g(i_high,j) = 2 * discretization_.get()->v(i_high,j) - vOld;
            }
        }
    }
    //bottom bound
    j_low = discretization_.get()->vJBegin() - 1;
    for(int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++){
        if(!geometry_.get()->isFluid(i,j_low) && geometry_.get()->isFluid(i,j_low + 1)){
            if(geometry_.get()->get_velocity(i,j_low).first == "NSW"){
                discretization_.get()->v(i,j_low) = 0;
                discretization_.get()->g(i,j_low) = discretization_.get()->v(i,j_low);
            }else if(geometry_.get()->get_velocity(i,j_low).first == "SLW"){
                discretization_.get()->v(i,j_low) = 0;
                discretization_.get()->g(i,j_low) = discretization_.get()->v(i,j_low);
            }else if(geometry_.get()->get_velocity(i,j_low).first == "IN"){
                discretization_.get()->v(i,j_low) = geometry_.get()->get_velocity(i,j_low).second[1]; //v_in(0,j*h-h/2) //TODO trafo von li/re und u/v? allg. dafuq?
                discretization_.get()->g(i,j_low) = discretization_.get()->v(i,j_low);
            }else if(geometry_.get()->get_velocity(i,j_low).first == "OUT"){
                double vOld = discretization_.get()->v(i,j_low);
                discretization_.get()->v(i,j_low) = discretization_.get()->v(i,j_low+1);
                discretization_.get()->g(i,j_low) = 2 * discretization_.get()->v(i,j_low) - vOld;
            }
        }
    }
    //top bound
    j_high = discretization_.get()->vJEnd() + 1;
    for(int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++){
        if(!geometry_.get()->isFluid(i,j_high) && geometry_.get()->isFluid(i,j_high - 1 )){
            if(geometry_.get()->get_velocity(i,j_high).first == "NSW"){
                discretization_.get()->v(i,j_high) = 0;
                discretization_.get()->g(i,j_high) = discretization_.get()->v(i,j_high);
            }else if(geometry_.get()->get_velocity(i,j_high).first == "SLW"){
                discretization_.get()->v(i,j_high) = 0;
                discretization_.get()->g(i,j_high) = discretization_.get()->v(i,j_high);
            }else if(geometry_.get()->get_velocity(i,j_high).first == "IN"){
                discretization_.get()->v(i,j_high) = geometry_.get()->get_velocity(i,j_high).second[1]; //v_in(0,j*h-h/2) //TODO trafo von li/re und u/v? allg. dafuq?
                discretization_.get()->g(i,j_high) = discretization_.get()->v(i,j_high);
            }else if(geometry_.get()->get_velocity(i,j_high).first == "OUT"){
                double vOld = discretization_.get()->v(i,j_high);
                discretization_.get()->v(i,j_high) = discretization_.get()->v(i,j_high - 1);
                discretization_.get()->g(i,j_high) = 2 * discretization_.get()->v(i,j_high) - vOld;
            }
        }
    }
    //inner cells //TODO korrekt?
    for(int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++){
        for(int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++){
            if(!geometry_.get()->isFluid(i,j)){
                if(!geometry_.get()->isFluid(i, j + 1)){
                    discretization_.get()->v(i,j) = 0;
                    discretization_.get()->g(i,j) = 0;
                }else if(geometry_.get()->isFluid(i + 1, j)){
                    discretization_.get()->v(i,j) = - discretization_.get()->v(i + 1, j);
                    discretization_.get()->g(i,j) = discretization_.get()->v(i,j);
                }else if(geometry_.get()->isFluid(i - 1, j)){
                    discretization_.get()->v(i,j) = - discretization_.get()->v(i - 1, j);
                    discretization_.get()->g(i,j) = discretization_.get()->v(i,j);
                }
            }else if(!geometry_.get()->isFluid(i, j + 1)){
                discretization_.get()->v(i,j) = 0;
                discretization_.get()->g(i,j) = 0;
            }
        }
    }
}

void Computation::applyBoundaryValuesTemperature(){
    // TODO Warum setzen wir hier auf einmal wieder Ecken über oben und unten? Hatten wir da nicht mal gesagt, dass wir das über links und rechts machen? Oder dnekt ihr, das hängt davon ab, wo der Fluss ist?

    // help variables
    int i_low;
    int i_high;
    int j_low;
    int j_high;

    //outer bounds
    //left bound without corners
    i_low = discretization_.get()->tIBegin()-1;
    for (int j = discretization_.get()->tJBegin() ; j <= discretization_.get()->tJEnd() ; j++) {
        if(!geometry_.get()->isFluid(i_low, j) && geometry_.get()->isFluid(i_low + 1, j)){
            if(geometry_.get()->get_temperature(i_low,j).first == "TN"){
                discretization_.get()->t(i_low, j) = discretization_.get()->t(i_low + 1, j) - discretization_.get()->dx() * geometry_.get()->get_temperature(i_low,j).second[0];
            }else{ //"TD"
                discretization_.get()->t(i_low, j) = 2 * geometry_.get()->get_temperature(i_low,j).second[0] - discretization_.get()->t(i_low + 1, j);
            }
        }
    }

    //right bound without corners
    i_high = discretization_.get()->tIEnd()+1;
    for (int j = discretization_.get()->tJBegin() ; j <= discretization_.get()->tJEnd() ; j++) {
        if(!geometry_.get()->isFluid(i_high, j) && geometry_.get()->isFluid(i_high - 1, j)){
            if(geometry_.get()->get_temperature(i_high,j).first == "TN"){
                discretization_.get()->t(i_high, j) = discretization_.get()->t(i_high - 1, j) - discretization_.get()->dx() * geometry_.get()->get_temperature(i_high,j).second[0];
            }else{
                discretization_.get()->t(i_high, j) = 2 * geometry_.get()->get_temperature(i_high,j).second[0] - discretization_.get()->t(i_high + 1, j);
            }
        }
    }

    //bottom bound with corners
    j_low = discretization_.get()->tJBegin()-1;
    for (int i = discretization_.get()->tIBegin()-1; i <= discretization_.get()->tIEnd()+1; i++) {
        if(!geometry_.get()->isFluid(i, j_low) && geometry_.get()->isFluid(i, j_low + 1)){
            if(geometry_.get()->get_temperature(i,j_low).first == "TN"){
                discretization_.get()->t(i, j_low) = discretization_.get()->t(i, j_low + 1) - discretization_.get()->dy() * geometry_.get()->get_temperature(i,j_low).second[0];
            }else{
                discretization_.get()->t(i, j_low) = 2 * geometry_.get()->get_temperature(i,j_low).second[0] - discretization_.get()->t(i, j_low + 1);
            }
        }
    }

    //top bound with corners
    j_high = discretization_.get()->tJEnd()-1;
    for (int i = discretization_.get()->tIBegin()-1; i <= discretization_.get()->tIEnd()+1; i++) {
        if(!geometry_.get()->isFluid(i, j_high) && geometry_.get()->isFluid(i, j_high - 1)){
            if(geometry_.get()->get_temperature(i,j_high).first == "TN"){
                discretization_.get()->t(i, j_high) = discretization_.get()->t(i, j_high - 1) - discretization_.get()->dy() * geometry_.get()->get_temperature(i,j_high).second[0];
            }else{
                discretization_.get()->t(i, j_high) = 2 * geometry_.get()->get_temperature(i,j_high).second[0] - discretization_.get()->t(i, j_high - 1);
            }
        }
    }

    //inner cells (obstacles)
    for(int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++){
        for(int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++){
            if(!geometry_.get()->isFluid(i,j)){
                int N = 0; //number of neighbouring fluid cells
                double tTmp = 0;

                //check left and right neighbours
                if(geometry_.get()->isFluid(i - 1, j)){
                    tTmp += discretization_.get()->t(i - 1, j);
                    N++;
                }else if(geometry_.get()->isFluid(i + 1, j)){
                    tTmp += discretization_.get()->t(i + 1, j);
                    N++;
                }
                //check top and bottom neighbours
                if(geometry_.get()->isFluid(i, j - 1)){
                    tTmp += discretization_.get()->t(i, j - 1);
                    N++;
                }else if(geometry_.get()->isFluid(i, j + 1)){
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
            if(geometry_.get()->isFluid(i,j) && geometry_.get()->isFluid(i + 1, j)){
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

    for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            if(geometry_.get()->isFluid(i,j) && geometry_.get()->isFluid(i, j + 1)){
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
            if(geometry_.get()->isFluid(i,j)){
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
            if(geometry_.get()->isFluid(i,j) && geometry_.get()->isFluid(i + 1, j)){
                discretization_.get()->u(i, j) =
                        discretization_.get()->f(i, j) - dt_ * discretization_.get()->computeDpDx(i, j);
            }
        }
    }

    for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            if(geometry_.get()->isFluid(i,j) && geometry_.get()->isFluid(i, j + 1)){
                discretization_.get()->v(i, j) =
                        discretization_.get()->g(i, j) - dt_ * discretization_.get()->computeDpDy(i, j);
            }
        }
    }
}

void Computation::computeTemperature() {
    for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) {
        for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) {
            if(geometry_.get()->isFluid(i,j)){
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
            if(geometry_.get()->isFluid(i,j) && geometry_.get()->isFluid(i + 1, j)){
                discretization_.get()->u(i, j) = uInit;
            }else{
                discretization_.get()->u(i, j) = std::nan("");
            }
        }
    }

    for (int j = discretization_.get()->vJBegin(); j <= discretization_.get()->vJEnd(); j++) {
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            if(geometry_.get()->isFluid(i,j) && geometry_.get()->isFluid(i, j + 1)){
                discretization_.get()->v(i, j) = vInit;
            }else{
                discretization_.get()->v(i,j) = std::nan("");
            }
        }
    }

    for (int j = discretization_.get()->pJBegin(); j <= discretization_.get()->pJEnd(); j++) {
        for (int i = discretization_.get()->pIBegin(); i <= discretization_.get()->pIEnd(); i++) {
            if(geometry_.get()->isFluid(i,j)){
                discretization_.get()->p(i, j) = pInit;
            }else{
                discretization_.get()->p(i,j) = std::nan("");
            }
        }
    }
for (int j = discretization_.get()->tJBegin(); j <= discretization_.get()->tJEnd(); j++) { for (int i = discretization_.get()->tIBegin(); i <= discretization_.get()->tIEnd(); i++) { if(geometry_.get()->isFluid(i,j)){
                discretization_.get()->t(i, j) = tInit;
            }else{
                discretization_.get()->t(i,j) = std::nan("");
            }
        }
    }
}
