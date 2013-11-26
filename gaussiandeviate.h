#include <cmath>

#ifndef GAUSSIANDEVIATE_H
#define GAUSSIANDEVIATE_H

#include <sstream>
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

double ran2(long *);
double gaussian_deviate(long *);

#endif
