//
// Created by Julia Pelzer on 26.10.2019.
//

#include "Computation.h"
#include "Settings.h"

using namespace std;

void Computation::intitialize(int argc, char **argv) {
    cout << "Running with" << settingsFilename << endl;
    Settings settings;
    settings.loadFromFile(argv[1]);
    settings_(settings);
}
