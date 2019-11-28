//
// Created by Fabian on 26.10.2019.
//

#include "Settings.h"
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>
#include "Util/Utils.h"

void Settings::loadFromFile(std::string filename) {
    std::ifstream file(filename.c_str(), std::ios::in);

    if (!file.is_open()) {
        std::cout << "Can not open file" << filename << std::endl;
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
        else if ( parameterName == "useDonorCell") {
            this->useDonorCell = strcmp(parameterValue, "true") == 0;
        }
        else if ( parameterName == "alpha") {
            this->alpha = atof(parameterValue);
        }
        else if ( parameterName == "tau") {
            this->tau = atof(parameterValue);
        }
        else if ( parameterName == "maximumDt") {
            this->maximumDt = atof(parameterValue);
        }
        else if ( parameterName == "pressureSolver") {
            this->pressureSolver = parameterValue;
        }
        else if ( parameterName == "omega") {
            this->omega = atof(parameterValue);
        }
        else if ( parameterName == "epsilon") {
            this->epsilon = atof(parameterValue);
        }
        else if ( parameterName == "maximumNumberOfIterations") {
            this->maximumNumberOfIterations = atof(parameterValue);
        }
        else if ( parameterName == "prandtl") {
            this->prandtl = atof(parameterValue);
        }
        else if ( parameterName == "beta") {
            this->beta = atof(parameterValue);
        }
        else if ( parameterName == "uInit") {
            this->uInit = atof(parameterValue);
        }
        else if ( parameterName == "vInit") {
            this->vInit = atof(parameterValue);
        }
        else if ( parameterName == "pInit") {
            this->pInit = atof(parameterValue);
        }
        else if ( parameterName == "tInit") {
            this->tInit = atof(parameterValue);
        }
    }
};

void Settings::printSettings() {
    std::cout << "Settings: " << std::endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x "
              << nCells[1] << std::endl
              << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau
              << ", maximum dt: " << maximumDt << std::endl
              << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"
              << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"
              << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
              << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
              << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
              << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon
              << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl
              << ", initial u: " << uInit << ", initial v: " << vInit << ", initial p: " << pInit << ", initial T: " << tInit << std::endl;             ;

}

std::string trim(const std::string &str,
                 const std::string &whitespace = " \t") {
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}
