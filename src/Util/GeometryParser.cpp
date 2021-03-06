//
// Created by fritz on 12/2/19.
//

#include <fstream>
#include <iostream>
#include <Settings.h>
#include <sstream>
#include <cstring>
#include <assert.h>
#include "GeometryParser.h"
#include "Utils.h"

shared_ptr<Geometry> GeometryParser::parseGeometryFile(std::string filename, Settings *settings) {
    std::ifstream file(filename.c_str(), std::ios::in);

    if (!file.is_open()) {
        std::cout << "Can not open file " << filename << std::endl;
        exit(1);
    }


    bool meshFound = false;
    int meshLineNumber = 1;
    while (!file.eof()) {
        std::string line;
        getline(file, line);
        line = Utils::trim(line);
        if (line[0] == '#') {
            continue;
        }
        if (line.find_first_of('=') == std::string::npos && !meshFound) {
            continue;
        }
        std::vector<std::string> values;
        Utils::split(line, values, '=');
        std::string parameterName = Utils::trim(values[0]);
        Utils::split(line, values, '=');
        std::string parameterValueString = Utils::trim(values[1]);
        if (line.find('#')) {
            std::vector<std::string> value_clean;
            Utils::split(values[1], value_clean, '#');
            parameterValueString = Utils::trim(value_clean[0]);
        }
        const char *parameterValue = parameterValueString.c_str();
        if (parameterName == "physicalSizeX") {
            settings->physicalSize[0] = atof(parameterValue);
        } else if (parameterName == "physicalSizeY") {
            settings->physicalSize[1] = atof(parameterValue);
        } else if (parameterName == "xOrigin") {
            settings->origin[0] = atof(parameterValue);
        } else if (parameterName == "yOrigin") {
            settings->origin[1] = atof(parameterValue);
        } else if (parameterName == "nCellsX") {
            settings->nCells[0] = atoi(parameterValue);
        } else if (parameterName == "nCellsY") {
            settings->nCells[1] = atoi(parameterValue);
        } else if (parameterName == "Mesh") {
            meshFound = true;
            geometry_ = make_shared<Geometry>(Geometry({settings->nCells[0]+2, settings->nCells[1]+2}));
            continue;
        }
        if (meshFound) {
            //cout << "Create geometry" << settings->nCells[0]  << "x" << settings->nCells[1] << endl;
            parseMeshLine(line, settings->nCells[1]+2-meshLineNumber);
            meshLineNumber++;

        }
    }
    geometry_.get()->countFluidCells();

    //cout << "testing" << endl;
    //cout << "out: " << geometry_.get()->get_velocity(0,1).first << endl;
    //cout << geometry_.get()->get_velocity(1,0).first << endl;
    return geometry_;
}

void GeometryParser::parseMeshLine(string line, int lineNumber) {
    std::stringstream lineStream(line);
    std::string cell;

    int columnNumber = 0;
    while (std::getline(lineStream, cell, ',')) {
        parseMeshCell(cell, columnNumber, lineNumber);
        columnNumber++;
    }
}

void GeometryParser::parseMeshCell(string basicString, int columnNumber, int lineNumber) {
    vector <string> boundaries;
    Utils::split(basicString, boundaries, ';');
    for (int i = 0; i < boundaries.size() ; ++i) {
        vector<string> values;

        if(boundaries[i].find(':')  == std::string::npos){
            values.push_back(boundaries[i]);
        }else{
            Utils::split(boundaries[i], values, ':');
        }
        //cout << "line:"<< lineNumber << " col: " << columnNumber << "val: " << values[0] << endl;

        geometry_.get()->set_state(columnNumber, lineNumber, {"S", {0}});

        if (values[0] == "NSW"){ // Noslip wall

            geometry_.get()->set_velocity(columnNumber, lineNumber, {"NSW", {0}});
        }
        else if (values[0] == "SLW"){ // Slip wall
            geometry_.get()->set_velocity(columnNumber, lineNumber, {"SLW", {0}});
        }
        else if (values[0] == "IN"){ // IN wall
            assert(values.size() == 3);
            double u = stod(values[1]);
            double v = stod(values[2]);
            geometry_.get()->set_velocity(columnNumber, lineNumber, {"IN", {u, v}});
        }
        else if (values[0] == "OUT"){ // Out
            geometry_.get()->set_velocity(columnNumber, lineNumber, {"OUT", {0}});
        }
        else if (values[0] == "PR"){ // Pressure
            double pr = stod(values[1]);
            geometry_.get()->set_pressure(columnNumber, lineNumber, {"PR", {pr}});

        }
        else if (values[0] == "TD"){ // T Dirichlet
            double td = stod(values[1]);
            geometry_.get()->set_temperature(columnNumber, lineNumber, {"TD", {td}});

        }
        else if (values[0] == "TPD"){ // preCICE T Dirichlet
            geometry_.get()->set_temperature(columnNumber, lineNumber, {"TPD", {0}});

        }
        else if (values[0] == "TPN"){ // preCICE T Neumann
            geometry_.get()->set_temperature(columnNumber, lineNumber, {"TPN", {0}});

        }
        else if (values[0] == "TN"){ // T Neumann
            double tn = stod(values[1]);
            geometry_.get()->set_temperature(columnNumber, lineNumber, {"TN", {tn}});
        }
        else if (values[0] == "F"){ // Fluid cell
            geometry_.get()->set_state(columnNumber, lineNumber, {"F", {0}});
        }
    }

}
