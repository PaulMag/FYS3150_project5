#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

#include <sstream>

using namespace std;
using namespace arma;

#include "CelObj.h"

CelObj:: CelObj() {
    /* Default dimensionality is 2. */
    position = zeros<rowvec>(2);
    velocity = zeros<rowvec>(2);
}

CelObj:: CelObj(int dim) {
    position = zeros<rowvec>(dim);
    velocity = zeros<rowvec>(dim);
}

CelObj:: CelObj(string n, double m, rowvec x, rowvec v) {
    name = n;
    mass = m;
    position = x;
    velocity = v;
}

CelObj:: ~CelObj() {
    /* Destructor is called all the time for some reason.
    if (outfileOpen) { // avoid closing if it was never opened
        outfile->close();
        cout << name << ".dat was closed." << endl;
    }
    */
}

rowvec CelObj:: getForce(CelObj other) {
    /* Find the force that act upon this body from another CelObj, NOT included
     * the gravity constant G.
     */
    rowvec r = other.position - position;
    return mass * other.mass * r / pow(norm(r,2), 3);
}

double CelObj:: getKineticEnergy() {
    return 0.5 * mass * pow(norm(velocity,2), 2); // 1/2 m*v^2
}

double CelObj:: getPotentialEnergy(CelObj other) {
    /* Find the potential energy on this body from the gravitational field of
     * another CelObj, NOT included the gravity constant G.
     */
    rowvec r = other.position - position;
    return - mass * other.mass / norm(r,2);
}

void CelObj:: makeOutfile(string location) {
    outfile  = new ofstream(location.c_str());
}

void CelObj:: closeOutfile() {
    outfile->close();
}

void CelObj:: writeData() {
    /* Currently only 2D writing is available. */
    //*outfile << position(0) << "," << position(1) << "," << position(2) << endl;
    *outfile << position(0) << "," << position(1) <<  endl;
}
