//
// Created by fritz on 11/1/19.
//

#ifndef NUMSIM_UTILS_H
#define NUMSIM_UTILS_H

#include <string>
#include <vector>


class Utils {
public:
    static std::string trim(const std::string &str, const std::string &whitespace = " \t");

    static void split(const std::string &str, std::vector<std::string> &cont, char delim);

};


#endif //NUMSIM_UTILS_H
