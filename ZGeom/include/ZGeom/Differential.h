#ifndef ZGEOM_DIFFERENTIAL_H
#define ZGEOM_DIFFERENTIAL_H

#include <functional>
#include "EigenSystem.h"
#include "VecN.h"

namespace ZGeom {

double calHeatKernel(const EigenSystem& es, int x, int y, double t);
double calKernel(const EigenSystem& es, std::function<double(double)> transferFunc, int x, int y);

}   // end of namespace


#endif