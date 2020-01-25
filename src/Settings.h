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

    double re = 1000;

    double endTime = 10.0;

    double tau = 0.5;

    double maximumDt = 0.1;

    std::array<double, 2> g{0., 0.};

    std::array<double, 2> dirichletBcBottom;

    std::array<double, 2> dirichletBcTop;

    std::array<double, 2> dirichletBcLeft;

    std::array<double, 2> dirichletBcRight;

    double prandtl;

    double tInit = 0;
    double fInit = 0; //TODO vielleicht in ein Feld ändern

    int nVelo = 9;

    std::shared_ptr<Geometry> geometry;

    double outputFileEveryDt = 0;
};


#endif //NUMSIM2019_SETTINGS_H
