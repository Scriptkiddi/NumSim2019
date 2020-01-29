//
// Created by Fabian on 26.10.2019.
//

#include "Settings.h"
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>
#include <Util/GeometryParser.h>
#include "Util/Utils.h"
#include <cmath>

void Settings::loadFromFile(std::string filename) {
    std::ifstream file(filename.c_str(), std::ios::in);

    if (!file.is_open()) {
        std::cout << "Can not open file " << filename << std::endl;
        exit(1);
        return;
    }

    while (!file.eof()) {
        std::string line;
        getline(file, line);
        line = Utils::trim(line);
        if (line[0] == '#') {
            continue;
        }
        if (line.find_first_of('=') == std::string::npos) {
            continue;
        }
        vector <string> values;
        Utils::split(line, values, '=');
        string parameterName = Utils::trim(values[0]);
        Utils::split(line, values, '=');
        string parameterValueString = Utils::trim(values[1]);
        if (line.find('#')) {
            vector <string> value_clean;
            Utils::split(values[1], value_clean, '#');
            parameterValueString = Utils::trim(value_clean[0]);
        }
        const char *parameterValue = parameterValueString.c_str();
        if (parameterName == "physicalSizeX") {
            this->physicalSize[0] = atof(parameterValue);
        }
        else if ( parameterName == "physicalSizeY") {
            this->physicalSize[1] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCTopX") {
            this->DirichletBCTop[0] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCTopY") {
            this->DirichletBCTop[1] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCBottomX") {
            this->DirichletBCBottom[0] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCBottomY") {
            this->DirichletBCBottom[1] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCLeftX") {
            this->DirichletBCLeft[0] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCLeftY") {
            this->DirichletBCLeft[1] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCRghtX") {
            this->DirichletBCRight[0] = atof(parameterValue);
        }
        else if ( parameterName == "DirichletBCRightY") {
            this->DirichletBCRight[1] = atof(parameterValue);
        }
        else if ( parameterName == "endTime") {
            this->endTime = atof(parameterValue);
        }
        else if ( parameterName == "nCellsX") {
            this->nCells[0] = atoi(parameterValue);
        }
        else if ( parameterName == "nCellsY") {
            this->nCells[1] = atoi(parameterValue);
        }
        else if ( parameterName == "timeStepRelaxation") {
            this->timeStepRelaxation = atof(parameterValue);
        }
        else if ( parameterName == "maximumDt") {
            this->maximumDt = atof(parameterValue);
        }
        else if ( parameterName == "fInit") {
            this->fInit = atof(parameterValue);
        }
        else if ( parameterName == "Temperature") {
            this->T = atof(parameterValue);
        }
        else if ( parameterName == "molarMass") {
            this->M = atof(parameterValue);
        }
        else if ( parameterName == "adiabaticIndex") {
            this->gamma = atof(parameterValue);
        }
        else if ( parameterName == "initialDensity") {
            this->rhoInit = atof(parameterValue);
        }
        else if ( parameterName == "viscosity") {
            this->viscosity = atof(parameterValue);
        }
        else if ( parameterName == "nVelo") {
            this->nVelo = atof(parameterValue);
        }
        else if (parameterName == "outputFileEveryDt"){
            this->outputFileEveryDt = atof(parameterValue);
        }
    }

    cs = sqrt(gamma * R * T / M);
    fInit = rhoInit / nVelo;
};

void Settings::printSettings() {
    std::cout << "Settings: " << std::endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x "
              << nCells[1] << std::endl
              << "  endTime: " << endTime << " s, timeStepRelaxation: " << timeStepRelaxation
              << ", maximum dt: " << maximumDt << std::endl
              << "  number of velocities per cell: " << nVelo << ", initial f: " << fInit << ", initial T: " << T << std::endl;

}

