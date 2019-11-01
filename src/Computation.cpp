//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Computation.h"
#include "Settings.h"
#include <cmath>

using namespace std;

void Computation::initialize(int argc, char **argv) {
    cout << "Running with" << argv[0] << endl;
    Settings settings;
    settings_ = settings;
    settings.loadFromFile(argv[1]);
    settings_.printSettings();
    //StaggeredGrid grid({2, 2}, {1, 1}); // einmal anlegen und füllen, dann nur noch überschreiben
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
}

void Computation::applyBoundaryValues() {


}

void Computation::PreliminaryVelocities() {
}

void Computation::computeRightHandSide() {
}

void Computation::computePressure() {
}

void Computation::computeVelocities() {
}
