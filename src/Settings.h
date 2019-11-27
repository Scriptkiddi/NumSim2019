//
// Created by Fabian on 26.10.2019.
//
using namespace std;

#ifndef NUMSIM2019_SETTINGS_H
#define NUMSIM2019_SETTINGS_H

#include <string>
#include <array>

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

    bool useDonorCell = false;

    double alpha = 0.5;

    std::array<double, 2> dirichletBcBottom;

    std::array<double, 2> dirichletBcTop;

    std::array<double, 2> dirichletBcLeft;

    std::array<double, 2> dirichletBcRight;

    std::string pressureSolver = "SOR";

    double omega = 1.0;

    double epsilon = 1e-5;

    int maximumNumberOfIterations = 1e5;

    double prandtl;

    double beta;

    double gamma = 1;

    double uInit;
    double vInit;
    double pInit;
    double tInit;
};


#endif //NUMSIM2019_SETTINGS_H
