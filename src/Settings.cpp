//
// Created by Fabian on 26.10.2019.
//

#include "Settings.h"
#include <fstream>
#include <iostream>

void loadFromFile(std::string settingsFilename){
    std::ifstream file(filename.c_str(), std::ios::in);

    if (!file.is_open()){
        std::cout << "Can not open file" << settingsFilename << std::endl;
    }

    while(!file.eof()){
        std::string line;
        getline(file, line);
        std::cout << line << std::endl;
    }
};
