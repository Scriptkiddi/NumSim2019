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
    settings.loadFromFile(argv[1]);
    settings_.printSettings();

    //initialize meshWidth
    meshWidth_[1] = settings_.physicalSize[1] /
                   (settings_.nCells[1] - 2); //todo stimmt das mit der Anzahl der Zellen? (mit 2 Ghostcells)
    meshWidth_[2] = settings_.physicalSize[2] / (settings_.nCells[2] - 2);

    //initialize discretization
    if (!settings_.useDonorCell) {
        CentralDifferences grid(settings_.nCells, meshWidth_);
        discretization_ = make_shared<CentralDifferences>(grid);
    } else {
        DonorCell grid(settings_.nCells, meshWidth_, settings_.alpha);
        discretization_ = make_shared<DonorCell>(grid);
    }

    //initialize explicit pressureSolver
    if (settings_.pressureSolver == "SOR") {
        SOR pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
        pressureSolver_ = make_unique<SOR>(pSolver);
    } else {
        //GaussSeidel pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
        //pressureSolver_ = make_unique<PressureSolver>(pSolver);
        std::cout << "Please select SOR-solver" << std::endl;
    }
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
    double uMaximum = discretization_->u(discretization_->uIBegin(), discretization_->uJEnd());
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
        for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) {
            if (uMaximum < discretization_->u(i, j)) {
                uMaximum = discretization_->u(i, j);
            }
        }
    }
    double vMaximum = discretization_->v(discretization_->vIBegin(), discretization_->vJBegin());
    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++) {
        for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
            if (vMaximum < discretization_->v(i, j)) {
                vMaximum = discretization_->v(i, j);
            }
        }
    }
    double condition_diffusion = pow(discretization_->dx(), 2) * pow(discretization_->dy(), 2) /
                                 (pow(discretization_->dx(), 2) + pow(discretization_->dy(), 2)) * settings_.re /
                                 2;
    double condition_convection1 = discretization_->dx() / uMaximum;
    double condition_convection2 = discretization_->dy() / vMaximum;

    //dt_ = std::min(condition_diffusion, condition_convection1, condition_convection2);
    //dt_ = std::min(settings_.maximumDt, dt_) * .9;//*0.9, damit echt kleiner
}

void Computation::applyBoundaryValues() {
    //todo: Zellen in den Ecken werden noch ignoriert, evtl. nach Bedarf ergänzen

    // u
    //unterer Rand
    int j = discretization_->uJBegin();
    for (int i = discretization_->uIBegin() + 1; i <= discretization_->uIEnd() - 1; i++) {
        discretization_->u(i, j) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i, j + 1);
    }
    //rechter und linker Rand
    int i_low = discretization_->uIBegin();
    int i_high = discretization_->uIEnd();
    for (j = discretization_->uJBegin() + 1; j <= discretization_->uJEnd() - 1; j++) {
        discretization_->u(i_low, j) = settings_.dirichletBcLeft[0];
        discretization_->u(i_high, j) = settings_.dirichletBcRight[0];
    }
    // oberer Rand
    j = discretization_->uJEnd();
    for (int i = discretization_->uIBegin() + 1; i <= discretization_->uIEnd() - 1; i++) {
        discretization_->u(i, j) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i, j - 1);
    }

    // v
    //unterer Rand
    j = discretization_->vJBegin();
    for (int i = discretization_->vIBegin() + 1; i <= discretization_->vIEnd() - 1; i++) {
        discretization_->v(i, j) = settings_.dirichletBcBottom[1] - discretization_->v(i, j + 1);
    }
    //rechter und linker Rand
    i_low = discretization_->vIBegin();
    i_high = discretization_->vIEnd();
    for (j = discretization_->vJBegin() + 1; j <= discretization_->vJEnd() - 1; j++) {
        discretization_->v(i_low, j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(i_low + 1, j);
        discretization_->v(i_high, j) =
                2 * settings_.dirichletBcRight[1] - discretization_->v(i_high - 1, j);
    }
    // oberer Rand
    j = discretization_->vJEnd();
    for (int i = discretization_->vIBegin() + 1; i <= discretization_->vIEnd() - 1; i++) {
        discretization_->v(i, j) = 2 * settings_.dirichletBcTop[1] - discretization_->v(i, j - 1);
    }

}

void Computation::PreliminaryVelocities() {
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
        for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) {
            discretization_->f(i, j) =
                    discretization_->u(i, j) + dt_ * (1 / settings_.re * (discretization_->computeD2uDx2(i, j) +
                                                                          discretization_->computeD2uDy2(i, j)) -
                                                      discretization_->computeDu2Dx(i, j) -
                                                      discretization_->computeDuvDy(i, j) + settings_.g[0]);
        }
    }

    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++) {
        for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
            discretization_->g(i, j) =
                    discretization_->v(i, j) + dt_ * (1 / settings_.re * (discretization_->computeD2vDy2(i, j) +
                                                                          discretization_->computeD2vDx2(i, j)) -
                                                      discretization_->computeDv2Dy(i, j) -
                                                      discretization_->computeDuvDx(i, j) + settings_.g[1]);
        }
    }
}

void Computation::computeRightHandSide() {
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
        for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) {
            discretization_->rhs(i, j) =
                    1 / dt_ * ((discretization_->f(i, j) - discretization_->f(i - 1, j)) / discretization_->dx() +
                               (discretization_->g(i, j) - discretization_->g(i, j - 1)) / discretization_->dy());
        }
    }
}

void Computation::computePressure() {

    pressureSolver_->solve();
}

void Computation::computeVelocities() {
    for (int j = discretization_->uJBegin(); j <= discretization_->uJEnd(); j++) {
        for (int i = discretization_->uIBegin(); i <= discretization_->uIEnd(); i++) {
            discretization_->u(i, j) = discretization_->f(i, j) - dt_ *discretization_->computeDpDx(i,j);
        }
    }

    for (int j = discretization_->vJBegin(); j <= discretization_->vJEnd(); j++) {
        for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
            discretization_->v(i, j) = discretization_->g(i, j) - dt_ *discretization_->computeDpDy(i,j);
        }
    }
}
