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
    settings_(settings);
    StaggeredGrid grid({2, 2}, {1, 1}); // einmal anlegen und füllen, dann nur noch überschreiben

    //initialize meshWidth
    std::array<int, 2> meshWidth;
    meshWidth[1] = settings_.physicalSize[1]/(settings_.nCells[1]-2); //todo stimmt das mit der Anzahl der Zellen? (mit 2 Ghostcells)
    meshWidth[2] = settings_.physicalSize[2]/(settings_.nCells[2]-2);
    meshWidth_ = meshWidth;

    //initialize discretization
    if (settings_.useDonorCell == "false") {
        CentralDifferences grid(settings_.nCells, meshWidth_);
    } else {
        DonorCell grid(settings_.nCells, meshWidth_, settings_.alpha);
    }
    discretization_ = grid;

    //initialize explicit pressureSolver
    if (settings_.pressureSolver == "SOR") {
        SOR pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else {
        GaussSeidel pSolver(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    }
    pressureSolver_ = pSolver;
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
    dt_ = std::min(settings_.maximumDt, dt_);
}

void Computation::applyBoundaryValues() {
    //todo: Zellen in den Ecken werden noch ignoriert, evtl. nach Bedarf ergänzen

    // u
    //unterer Rand
    int j = discretization_.get()->uJBegin();
    for (int i = discretization_.get()->uIBegin() + 1; i <= discretization_.get()->uIEnd() - 1; i++) {
        discretization_.get()->u(i, j) = 2 * settings_.dirichletBcBottom[0] - discretization_.get()->u(i, j + 1);
    }
    //rechter und linker Rand
    int i_low = discretization_.get()->uIBegin();
    int i_high = discretization_.get()->uIEnd();
    for (j = discretization_.get()->uJBegin() + 1; j <= discretization_.get()->uJEnd() - 1; j++) {
        discretization_->u(i_low, j) = settings_.dirichletBcLeft[0];
        discretization_.get()->u(i_high, j) = settings_.dirichletBcRight[0];
    }
    // oberer Rand
    j = discretization_.get()->uJEnd();
    for (int i = discretization_.get()->uIBegin() + 1; i <= discretization_.get()->uIEnd() - 1; i++) {
        discretization_.get()->u(i, j) = 2 * settings_.dirichletBcTop[0] - discretization_.get()->u(i, j - 1);
    }

    // v
    //unterer Rand
    int j = discretization_.get()->vJBegin();
    for (int i = discretization_.get()->vIBegin() + 1; i <= discretization_.get()->vIEnd() - 1; i++) {
        discretization_.get()->v(i, j) = settings_.dirichletBcBottom[1] - discretization_.get()->v(i, j + 1);
    }
    //rechter und linker Rand
    int i_low = discretization_.get()->vIBegin();
    int i_high = discretization_.get()->vIEnd();
    for (j = discretization_.get()->vJBegin() + 1; j <= discretization_.get()->vJEnd() - 1; j++) {
        discretization_.get()->v(i_low, j) = 2 * settings_.dirichletBcLeft[1] - discretization_.get()->v(i_low + 1, j);
        discretization_.get()->v(i_high, j) =
                2 * settings_.dirichletBcRight[1] - discretization_.get()->v(i_high - 1, j);
    }
    // oberer Rand
    j = discretization_.get()->vJEnd();
    for (int i = discretization_.get()->vIBegin() + 1; i <= discretization_.get()->vIEnd() - 1; i++) {
        discretization_.get()->v(i, j) = 2 * settings_.dirichletBcTop[1] - discretization_.get()->v(i, j - 1);
    }

}

void Computation::PreliminaryVelocities() {
    for (int j = grid.uJBegin; j <= grid.uJEnd; j++) {
        for (int i = grid.uIBegin; i <= grid.uIEnd; i++) {
            grid.f(i, j) = grid.u(i, j) + dt_ * (1 / settings_.re * (discretization_.computeD2uDx2(i, j) +
                                                                     discretization_.computeD2uDy2(i, j)) -
                                                 discretization_.computeDu2Dx(i, j) -
                                                 discretization_.computeDuvDy(i, j) + settings_.g[0]);
        }
    }

    for (int j = grid.vJBegin; j <= grid.vJEnd; j++) {
        for (int i = grid.vIBegin; i <= grid.vIEnd; i++) {
            grid.g(i, j) = grid.v(i, j) + dt_ * (1 / settings_.re * (discretization_.computeD2vDy2(i, j) +
                                                                     discretization_.computeD2vDx2(i, j)) -
                                                 discretization_.computeDv2Dy(i, j) -
                                                 discretization_.computeDuvDx(i, j) + settings_.g[1]);
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

    pressureSolver_.solve();
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
