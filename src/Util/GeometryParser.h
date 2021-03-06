//
// Created by fritz on 12/2/19.
//

#ifndef NUMSIM2019_GEOMETRYPARSER_H
#define NUMSIM2019_GEOMETRYPARSER_H


#include <Geometry.h>
#include <Settings.h>

class GeometryParser {
public:


    std::shared_ptr<Geometry> parseGeometryFile(std::string filePath, Settings *settings);

    void parseMeshLine(string line, int lineNumber);

    void parseMeshCell(string basicString, int lineNumber);

private:
    std::shared_ptr<Geometry> geometry_;

    void parseMeshCell(string basicString, int columnNumber, int lineNumber);
};


#endif //NUMSIM2019_GEOMETRYPARSER_H
