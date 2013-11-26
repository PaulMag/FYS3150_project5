#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

#ifndef SOLSYS_H
#define SOLSYS_H

#include "CelObj.h"

#include <sstream>
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

class SolSys {

public:
    int N; // no of bodies
    int dim; // dimensionality
    double G; // gravitational constant

    vector<CelObj> bodies;

    ofstream* outfile;

    SolSys();
    SolSys(int);
    ~SolSys();

    void addCelObj(CelObj);
    void addCelObj(string, double, rowvec, rowvec);
    void addCelObj(string, double, double, double, double,
                                   double, double, double);
    void addCelObj(string, double, double, double,
                                   double, double);

    void setPositions (mat);
    void setVelocities(mat);

    mat getPositions ();
    mat getVelocities();

    rowvec getCenterOfMass();
    void   setCenterOfMass();

    rowvec getTotalMomentum();
    void   setTotalMomentum(CelObj*);
    void   setTotalMomentum();

    cube findForces();
    mat  findAccels();

    double getPotentialEnergy ();
    double getPotentialEnergy (int);
    vec getPotentialEnergies  ();
    double getKineticEnergy   ();
    double getKineticEnergy   (int);
    vec getKineticEnergies    ();
    double getEquilibriumEnergy();

    void makeDataFiles();
    void makeDataFiles(string);

    void rungeKutta4();
    void rungeKutta4(double h);
    void leapFrog   (double h);

    void moveSystem(double, int,    string);
    void moveSystem(double, double, string);
    void moveSystem(double, int);
    void moveSystem(double, double);
};

#endif
