//
// Created by fritz on 11/12/19.
//

#include <mpi.h>
#include <StaggeredGrid/CentralDifferences.h>
#include <StaggeredGrid/DonorCell.h>
#include <Communication.h>
#include "ComputationParallel.h"
#include <PressureSolver/GaussSeidel.h>
#include <PressureSolver/SOR.h>


ComputationParallel::ComputationParallel(std::string settingsFilename) : Computation(settingsFilename),
                                                                         partitioning_(
                                                                                 Partitioning(settings_.nCells, settings_.physicalSize)) {

}

void ComputationParallel::initialize(int argc, char **argv) {
    array<int, 2> nCellsBoundary = {partitioning_.getNCells()[0] + 2,
                                    partitioning_.getNCells()[1] + 2};

    //initialize meshWidth
    meshWidth_ = partitioning_.getMeshWidth();

    //init grid
    int nX = partitioning_.getNCells()[0];
    int nY = partitioning_.getNCells()[1];

    int nUX;
    int nVX;
    int nPX;
    int nPY;
    int nUY;
    int nVY;
    double offsetX = 0;
    double offsetY = 0;

    if (partitioning_.getRankOfLeftNeighbour() == -1) {
        nUX = nX + 1;
        nVX = nX + 2;
        nPX = nX + 2;
    }else{
        nUX = nX + 2;
        nVX = nX + 3;
        nPX = nX + 3;
        offsetX = -meshWidth_[0];
    }
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        nUY = nY + 2;
        nVY = nY + 1;
        nPY = nY + 2;
    }
    else{
        nUY = nY + 3;
        nVY = nY + 2;
        nPY = nY + 3;
        offsetY = -meshWidth_[1];
    }
    
    FieldVariable u = FieldVariable({nUX, nUY}, {0 * meshWidth_[0] + offsetX, -0.5 * meshWidth_[1] + offsetY}, meshWidth_);
    FieldVariable f = FieldVariable({nUX, nUY}, {-0.5 * meshWidth_[0]+ offsetX, -0.5 * meshWidth_[1] + offsetY}, meshWidth_);
    FieldVariable v = FieldVariable({nVX, nVY}, {-0.5 * meshWidth_[0]+ offsetX, 0 * meshWidth_[1] + offsetY}, meshWidth_);
    FieldVariable g = FieldVariable({nVX, nVY}, {-0.5 * meshWidth_[0]+ offsetX, -0.5 * meshWidth_[1] + offsetY}, meshWidth_);
    FieldVariable p = FieldVariable({nPX, nPY}, {-0.5 * meshWidth_[0]+ offsetX, -0.5 * meshWidth_[1] + offsetY}, meshWidth_);
    FieldVariable rhs = FieldVariable({nPX, nPY}, {-0.5 * meshWidth_[0]+ offsetX, -0.5 * meshWidth_[1] + offsetY}, meshWidth_);

    std::cout << partitioning_.getRank() <<"| u(" << u.size()[0] << "," << u.size()[1] << ")" << std::endl;
    std::cout << partitioning_.getRank() <<"| v(" << v.size()[0] << "," << v.size()[1] << ")" << std::endl;
    std::cout << partitioning_.getRank() <<"| f(" << f.size()[0] << "," << f.size()[1] << ")" << std::endl;
    std::cout << partitioning_.getRank() <<"| g(" << g.size()[0] << "," << g.size()[1] << ")" << std::endl;
    std::cout << partitioning_.getRank() <<"| p(" << p.size()[0] << "," << p.size()[1] << ")" << std::endl;
    std::cout << partitioning_.getRank() <<"| r(" << rhs.size()[0] << "," << rhs.size()[1] << ")" << std::endl;

    //initialize discretization
    if (!settings_.useDonorCell) {
        CentralDifferences grid(nCellsBoundary, meshWidth_, make_shared<Partitioning>(partitioning_), u, v, p, f, g,
                                rhs);
        discretization_ = make_shared<CentralDifferences>(grid);
        std::cout << "use Central Differences" << std::endl;
    } else {
        DonorCell grid(nCellsBoundary, meshWidth_, settings_.alpha, make_shared<Partitioning>(partitioning_), u, v, p,
                       f, g, rhs);
        discretization_ = make_shared<DonorCell>(grid);
        std::cout << "use Donor Cell" << std::endl;
    }
    communication_ = make_shared<Communication>(make_shared<Partitioning>(partitioning_), discretization_);

    //initialize explicit pressureSolver
    if (settings_.pressureSolver == "SOR"){
        SOR pSolver(discretization_, communication_, partitioning_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
       pressureSolver_ = make_unique<SOR>(pSolver);
       std::cout << "Pressuresolver = SOR" << std::endl;
    } else {
        GaussSeidel pSolver(discretization_, communication_, partitioning_, settings_.epsilon, settings_.maximumNumberOfIterations);
        pressureSolver_ = make_unique<GaussSeidel>(pSolver);
        std::cout << "Pressuresolver = GaussSeidel" << std::endl;
    }

    //initialize outputWriters
    OutputWriterTextParallel outText(discretization_, partitioning_);
    outputWriterText_ = make_unique<OutputWriterTextParallel>(outText);

    OutputWriterParaviewParallel outputWriterParaviewParallel(discretization_, partitioning_);
    outputWriterParaview_ = make_unique<OutputWriterParaviewParallel>(outputWriterParaviewParallel);
}

void ComputationParallel::runSimulation() {
    double t = 0;

    applyBoundaryValues();

    cout << "Starting Simulation" << endl;
    outputWriterParaview_.get()->writeFile(t);
    outputWriterText_.get()->writeFile(t);
    double counter = 0;
    while (t < settings_.endTime) {
        computeTimeStepWidth();
        applyBoundaryValues();
        PreliminaryVelocities();
        //cout << "Communicating f&g" << endl;
        communication_.get()->communicate(discretization_.get()->f(), "f");
        communication_.get()->communicate(discretization_.get()->g(), "g");
        computeRightHandSide();
        //cout << "Computing p" << endl;
        computePressure();
        //cout << "Computing u&v" << endl;
        computeVelocities();
        //cout << "Communicating u&v" << endl;
        communication_.get()->communicate(discretization_.get()->u(), "u");
        communication_.get()->communicate(discretization_.get()->v(), "v");
        t += dt_;
        counter += dt_;
        if (counter >= 1){
            counter = 0;
            cout << "Writing paraview timestep: "<< t << endl;
            outputWriterParaview_.get()->writeFile(t);
        //    outputWriterText_.get()->writeFile(t);
        }
        //cout << partitioning_.getRank() << "|current time: " << t << " dt: " << dt_ << " pressure solver iterations: "
        //     << endl;

    }


}


void ComputationParallel::computeTimeStepWidth() {
    double uMaximum = 0;
    for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
        for (int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd() + 1; i++) {
            if (uMaximum < abs(discretization_.get()->u(i, j))) {
                uMaximum = abs(discretization_.get()->u(i, j));
            }
        }
    }
    double vMaximum = 0;
    for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
        for (int i = discretization_.get()->vIBegin() - 1; i <= discretization_.get()->vIEnd() + 1; i++) {
            if (vMaximum < abs(discretization_.get()->v(i, j))) {
                vMaximum = abs(discretization_.get()->v(i, j));
            }
        }
    }
    double condition_diffusion = pow(discretization_.get()->dx(), 2) * pow(discretization_.get()->dy(), 2) /
                                 (pow(discretization_.get()->dx(), 2) + pow(discretization_.get()->dy(), 2)) *
                                 settings_.re / 2;
    double condition_convection1 = discretization_.get()->dx() / abs(uMaximum);
    double condition_convection2 = discretization_.get()->dy() / abs(vMaximum);
    /*
    std::cout << partitioning_.getRank() << "|condition_convection1 " << condition_convection1 << std::endl;
    std::cout << partitioning_.getRank() << "|condition_convection2 " << condition_convection2 << std::endl;
    std::cout << partitioning_.getRank() << "|condition_diffusion " << condition_diffusion << std::endl;
    std::cout << partitioning_.getRank() << "|settings_.maximumDt " << settings_.maximumDt << std::endl;
    */
    dt_ = min(condition_convection1, condition_convection2);
    //std::cout << partitioning_.getRank() << "|dt_" << dt_ << std::endl;
    dt_ = min(condition_diffusion, dt_);
    //std::cout << partitioning_.getRank() << "|dt_" << dt_ << std::endl;
    dt_ = min(settings_.maximumDt, dt_) * settings_.tau;
    //std::cout << partitioning_.getRank() << "|dt_" << dt_ << std::endl;

    MPI_Allreduce(&dt_, &dtAll_, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt_ = dtAll_;
    //std::cout << partitioning_.getRank() << "| finales dt_" << dt_ << std::endl;
//TODO Warum ist dt nicht konstant bei 0.05??

//TODO pressure ist viel zu klein!
//TODO jetzt: v zurückverfolgen, da tritt der Fehler offensichtlich auf.
}

void ComputationParallel::applyBoundaryValues() {
    // U
    //Linker Rand
    if (partitioning_.getRankOfLeftNeighbour() == -1) {
        int i_low = discretization_.get()->uIBegin() - 1;
        for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
            discretization_.get()->u(i_low, j) = settings_.dirichletBcLeft[0];
            discretization_.get()->f(i_low, j) = discretization_.get()->u(i_low, j);
        }
    }

    // Rechter Rand
    if (partitioning_.getRankOfRightNeighbour() == -1) {
        int i_high = discretization_.get()->uIEnd() + 1;
        for (int j = discretization_.get()->uJBegin() - 1; j <= discretization_.get()->uJEnd() + 1; j++) {
            discretization_.get()->u(i_high, j) = settings_.dirichletBcRight[0];
            discretization_.get()->f(i_high, j) = discretization_.get()->u(i_high, j);
        }
    }

    //unterer Rand
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        int j_low = discretization_.get()->uJBegin() - 1;
        for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
            discretization_.get()->u(i, j_low) =
                    2 * settings_.dirichletBcBottom[0] - discretization_.get()->u(i, j_low + 1);
            discretization_.get()->f(i, j_low) = discretization_.get()->u(i, j_low);
        }
    }
    // oberer Rand
    if (partitioning_.getRankOfTopNeighbour() == -1) {
        if (partitioning_.getRankOfLeftNeighbour() == -1 && partitioning_.getRankOfRightNeighbour() == -1) {
            int j_high = discretization_.get()->uJEnd() + 1;
            for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd(); i++) {
                discretization_.get()->u(i, j_high) =
                        2 * settings_.dirichletBcTop[0] - discretization_.get()->u(i, j_high - 1);
                discretization_.get()->f(i, j_high) = discretization_.get()->u(i, j_high);
            }
        } else if (partitioning_.getRankOfLeftNeighbour() == -1) {
            int j_high = discretization_.get()->uJEnd() + 1;
            for (int i = discretization_.get()->uIBegin(); i <= discretization_.get()->uIEnd() + 1; i++) {
                discretization_.get()->u(i, j_high) =
                        2 * settings_.dirichletBcTop[0] - discretization_.get()->u(i, j_high - 1);
                discretization_.get()->f(i, j_high) = discretization_.get()->u(i, j_high);
            }
        } else if (partitioning_.getRankOfRightNeighbour() == -1) {
            int j_high = discretization_.get()->uJEnd() + 1;
            for (int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd(); i++) {
                discretization_.get()->u(i, j_high) =
                        2 * settings_.dirichletBcTop[0] - discretization_.get()->u(i, j_high - 1);
                discretization_.get()->f(i, j_high) = discretization_.get()->u(i, j_high);
            }
        } else {
            int j_high = discretization_.get()->uJEnd() + 1;
            for (int i = discretization_.get()->uIBegin() - 1; i <= discretization_.get()->uIEnd() + 1; i++) {
                discretization_.get()->u(i, j_high) =
                        2 * settings_.dirichletBcTop[0] - discretization_.get()->u(i, j_high - 1);
                discretization_.get()->f(i, j_high) = discretization_.get()->u(i, j_high);
            }
        }
    }

    // v
    //unterer Rand
    // oberer Rand
    if (partitioning_.getRankOfTopNeighbour() == -1) {
        int j_high = discretization_.get()->vJEnd() + 1;
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            discretization_.get()->v(i, j_high) = settings_.dirichletBcTop[1];
            discretization_.get()->g(i, j_high) = discretization_.get()->v(i, j_high);
        }
    }
    if (partitioning_.getRankOfBottomNeighbour() == -1) {
        int j_low = discretization_.get()->vJBegin() - 1;
        for (int i = discretization_.get()->vIBegin(); i <= discretization_.get()->vIEnd(); i++) {
            discretization_.get()->v(i, j_low) = settings_.dirichletBcBottom[1];
            discretization_.get()->g(i, j_low) = discretization_.get()->v(i, j_low);
        }
    }

    //rechter und linker Rand
    if (partitioning_.getRankOfRightNeighbour() == -1) {
        int i_high = discretization_.get()->vIEnd() + 1;
        for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
            discretization_.get()->v(i_high, j) =
                    2 * settings_.dirichletBcRight[1] - discretization_.get()->v(i_high - 1, j);
            discretization_.get()->g(i_high, j) = discretization_.get()->v(i_high, j);
        }
    }
    if (partitioning_.getRankOfLeftNeighbour() == -1) {
        int i_low = discretization_.get()->vIBegin() - 1;
        for (int j = discretization_.get()->vJBegin() - 1; j <= discretization_.get()->vJEnd() + 1; j++) {
            discretization_.get()->v(i_low, j) =
                    2 * settings_.dirichletBcLeft[1] - discretization_.get()->v(i_low + 1, j);
            discretization_.get()->g(i_low, j) = discretization_.get()->v(i_low, j);

        }
    }
}
