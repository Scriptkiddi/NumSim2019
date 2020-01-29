//
// Created by Fabian on 26.10.2019.
//
using namespace std;

#ifndef NUMSIM2019_SETTINGS_H
#define NUMSIM2019_SETTINGS_H

#include <string>
#include <array>
#include <memory>
#include "Geometry.h"

struct Settings {
    void loadFromFile(std::string filename);

    void printSettings();

    std::array<int, 2> nCells;

    std::array<double, 2> physicalSize;

    std::array<double, 2> DirichletBCTop = {0.0,0.0};

    std::array<double, 2> DirichletBCRight = {0.0,0.0};

    std::array<double, 2> DirichletBCBottom = {0.0,0.0};

    std::array<double, 2> DirichletBCLeft = {0.0,0.0};

    double endTime = 10.0;

    double timeStepRelaxation = 0.5;

    double maximumDt = 0.1;
    
    double fInit = 0; //TODO vielleicht in ein Feld ändern

    int nVelo = 9;
    
    double cs; //Speed of sound in [m/s]

    double rhoInit; //Density

    double viscosity = 0.0000185; //dynamic viscosity

    double M = 0.0289645; //molar mass

    double T = 293.15; //Temperature

    double R = 8.31446261; //universal Gas constant

    double gamma = 1.4; //adiabatic coefficient
    
    std::shared_ptr<Geometry> geometry;

    double outputFileEveryDt = 0;
};


#endif //NUMSIM2019_SETTINGS_H
