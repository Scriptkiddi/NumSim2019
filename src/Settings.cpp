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
        else if ( parameterName == "endTime") {
            this->endTime = atof(parameterValue);
        }
        else if ( parameterName == "re") {
            this->re = atof(parameterValue);
        }
        else if ( parameterName == "gX") {
            this->g[0] = atof(parameterValue);
        }
        else if ( parameterName == "gY") {
            this->g[1] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletBottomX") {
            this->dirichletBcBottom[0] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletBottomY") {
            this->dirichletBcBottom[1] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletTopX") {
            this->dirichletBcTop[0] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletTopY") {
            this->dirichletBcTop[1] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletLeftX") {
            this->dirichletBcLeft[0] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletLeftY") {
            this->dirichletBcLeft[1] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletRightX") {
            this->dirichletBcRight[0] = atof(parameterValue);
        }
        else if ( parameterName == "dirichletRightY") {
            this->dirichletBcRight[1] = atof(parameterValue);
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
        else if ( parameterName == "prandtl") {
            this->prandtl = atof(parameterValue);
        }
        else if ( parameterName == "fInit") {
            this->fInit = atof(parameterValue);
        }
        else if ( parameterName == "tInit") {
            this->tInit = atof(parameterValue);
        }
        else if ( parameterName == "nVelo") {
            this->nVelo = atof(parameterValue);
        }
        else if ( parameterName == "geometryFile"){
            cout << "create parser" << endl;
            cout << parameterValue << endl;
            GeometryParser parser = GeometryParser();
            //if( strcmp(parameterValue, "free")){
            //    this->geometry = nullptr;
            //}else{
                cout << "start parsing" << endl;
                this->geometry = parser.parseGeometryFile(parameterValue, this);
            //}
        } else if (parameterName == "outputFileEveryDt"){
            this->outputFileEveryDt = atof(parameterValue);
        }
    }

    
    cs = sqrt(gamma * R * T / M);
};

void Settings::printSettings() {
    std::cout << "Settings: " << std::endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x "
              << nCells[1] << std::endl
              << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), timeStepRelaxation: " << timeStepRelaxation
              << ", maximum dt: " << maximumDt << std::endl
              << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"
              << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"
              << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
              << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
              << ", number of velocities per cell: " << nVelo << ", initial f: " << fInit << ", initial T: " << tInit << std::endl;

}

