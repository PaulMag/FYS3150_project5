#include <iostream>
#include <cmath>
#include <fstream>
#include <armadillo>

#include <ctime>

#include "CelObj.h"
#include "SolSys.h"

#include "gaussiandeviate.h"

using namespace std;
using namespace arma;

/* Stuff to remember:
 *
 * mat A = zeros<mat>(n, n);
 * vec a(n);
 */

int main(int argc, char* argv[]) {

    clock_t start, finish;
    start = clock();

    SolSys mycluster = SolSys(3);

    int N; // number of stars
    istringstream NS(argv[1]);
    NS >> N;

    double totalMass = 100 * 10.0;

    double mu = totalMass / N; // [M_sun] // avg mass
    double stdev = mu / 10.;   // [M_sun] // standard deviation in mass


    double pi = 3.14159265359;
    double R0 = 20.;   // [l.y.] // radius of sphere
    double G_tau = pi*pi * R0*R0*R0 / (8 * N * mu);
    mycluster.G = G_tau; // [M_sun * tau_crunch^2 / (l.y)^3] = [M_sun / (l.y)^3]

    mycluster.setEps(0.15);

    time_t seed = time(0);
    long seed1 = - 1055  - seed;
    long seed2 = - 13937 - seed;

    double m;
    rowvec x;
    rowvec v = zeros<rowvec>(3);

    for (int i=0; i<N;) {
        x << ran2(&seed1) << ran2(&seed1) << ran2(&seed1);
        x -= 0.5; // center on origo
        if (norm(x,2) <= 0.5) {
            /* Stars are only allowed to be within this sphere. */
            m = gaussian_deviate(&seed2) * stdev + mu;
            mycluster.addCelObj("", m, x*2*R0, v);
            i++; // star was successfully addded
        }
    }

    /*
    mycluster.setCenterOfMass();  // Place CM in origo
    mycluster.setTotalMomentum(); // Set v_CM = 0
    */

    istringstream timeS(argv[2]);
    istringstream hS   (argv[3]);
    double time, h;
    timeS >> time;
    hS    >> h;
    string location = argv[4];

    mycluster.moveSystem(time, h, location);
    finish = clock();

    double compTime = double(finish - start) / CLOCKS_PER_SEC;
    cout << "Total computation time:  " << compTime  << " seconds." << endl;
}
