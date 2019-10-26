//
// Created by Fabian on 26.10.2019.
//

#ifndef NUMSIM2019_RESULT_H
#define NUMSIM2019_RESULT_H

#include "Array2D.h"

class Result {
private:
    int timestep;
    double time;
    Array2D u;
    Array2D v;
    Array2D p;

public:
    Result(double time, int Nx, int Ny, int timestep);

    double get(char name, int i, int j);

    Array2D get_all(char name);

    void set(char name, int i, int j, double value);

    void set_time(double time);

    double get_time();
};


#endif //NUMSIM2019_RESULT_H
