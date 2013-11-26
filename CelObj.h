#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

#ifndef CELOBJ_H
#define CELOBJ_H

class CelObj {

public:
    string name;
    double mass;

    rowvec position;
    rowvec velocity;

    ofstream* outfile;

    CelObj();
    CelObj(int);
    CelObj(string, double, rowvec, rowvec);
    ~CelObj();

    rowvec getForce(CelObj);

    double getPotentialEnergy (CelObj);
    double getKineticEnergy   ();

    void makeOutfile(string);
    void closeOutfile();
    void writeData();
};

#endif
