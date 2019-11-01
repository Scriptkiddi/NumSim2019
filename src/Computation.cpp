//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Computation.h"
#include "Settings.h"
#include <cmath>

using namespace std;

void Computation::intitialize(int argc, char **argv) {
    cout << "Running with" << settingsFilename << endl;
    Settings settings;
    settings.loadFromFile(argv[1]);
    settings_(settings);
    StaggeredGrid grid({2, 2}, {1, 1}); // einmal anlegen und füllen, dann nur noch überschreiben
}

void Computation::runSimulation() {
    // initialize(); //todo fill
    double t = 0;
    while (t < settings_.endTime) {
        computeTimeStepWidth();
        applyBoundaryValues();
        PreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();
        //outputwriter oÄ aufrufen
        t += dt_;
    }
}

void Computation::computeTimeStepWidth() {
    double uMaximum = grid.u(0, 0);
    for (int j = grid.uJBegin; j <= grid.uJEnd; j++) {
        for (int i = grid.uIBegin; i <= grid.uIEnd; i++) {
            if (uMaximum < grid.u(i, j)) {
                uMaximum = grid.u(i, j);
            }
        }
    }
    double vMaximum = grid.u(0, 0);
    for (int j = grid.uJBegin; j <= grid.uJEnd; j++) {
        for (int i = grid.uIBegin; i <= grid.uIEnd; i++) {
            if (vMaximum < grid.u(i, j)) {
                vMaximum = grid.u(i, j);
            }
        }
    }
    double condition_diffusion = pow(grid.dx()) * pow(grid.dy()) / (pow(grid.dx()) + pow(grid.dy())) * Re_ / 2 *
                                 0.99; //*0.99, damit echt kleiner
    double condition_convection1 = grid.dx() / uMaximum;
    double condition_convection2 = grid.dy() / vMaximum;

    dt_ = std::min(condition_diffusion, condition_convection1, condition_convection2);
    dt_ = std::max(settings_.maximumDt, dt_); //todo korrekt?
}

void Computation::applyBoundaryValues() {


}

void Computation::PreliminaryVelocities() {
    for (int j = grid.uJBegin; j <= grid.uJEnd; j++) {
        for (int i = grid.uIBegin; i <= grid.uIEnd; i++) {
            grid.f(i, j) = grid.u(i,j) + dt_ * (1/Re_ * (discretization_.computeD2uDx2(i,j) + discretization_.computeD2uDy2(i,j)) - discretization_.computeDu2Dx(i,j) - discretization_.computeDuvDy(i,j) + settings_.g[0]);
        }
    }

    for (int j = grid.vJBegin; j <= grid.vJEnd; j++) {
        for (int i = grid.vIBegin; i <= grid.vIEnd; i++) {
            grid.g(i, j) = grid.v(i,j) + dt_ * (1/Re_ * (discretization_.computeD2vDy2(i,j) + discretization_.computeD2vDx2(i,j)) - discretization_.computeDv2Dy(i,j) - discretization_.computeDuvDx(i,j) + settings_.g[1]);
        }
    }
}

void Computation::computeRightHandSide() {
    for (int j = grid.uJBegin; j <= grid.uJEnd; j++) {
        for (int i = grid.uIBegin; i <= grid.uIEnd; i++) {
            grid.rhs(i, j) = 1 / dt_ * ((grid.f(i, j) - grid.f(i - 1, j)) / grid.dx() +
                                        (grid.g(i, j) - grid.g(i, j - 1)) / grid.dy());
        }
    }
}

void Computation::computePressure() {
    while (it < itmax && residuum > epsilon) {
        //todo compute sor(..);
    }
}

void Computation::computeVelocities() {
    for (int j = grid.uJBegin; j <= grid.uJEnd; j++) {
        for (int i = grid.uIBegin; i <= grid.uIEnd; i++) {
            grid.u(i, j) = grid.f(i, j) - dt_ *;// todo dp/dx (i,j) von Schritt n+1
        }
    }

    for (int j = grid.vJBegin; j <= grid.vJEnd; j++) {
        for (int i = grid.vIBegin; i <= grid.vIEnd; i++) {
            grid.v(i, j) = grid.g(i, j) - dt_ *;//todo  dp/dx (i,j) von Schritt n+1
        }
    }
}
