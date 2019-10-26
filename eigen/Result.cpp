//
// Created by Fabian on 26.10.2019.
//

#include "Result.h"

Result(double time, int Nx, int Ny, int timestep) {
    this.time = time;
    this.timestep = timestep;
    Array2D
    this.u({Nx, Ny});
    Array2D
    this.v({Nx, Ny});
    Array2D
    this.p({Nx, Ny});
}

double get(char name, int i, int j) {
    assert(name == u || name == v || name == p);
    return this.name(i, j);
}

Array2D get_all(char name) {
    return this.name;
}

void set(char name, int i, int j, double value) {
    this.name(i, j) = value;
}

void set_time(double time) {
    this.time = time;
}

