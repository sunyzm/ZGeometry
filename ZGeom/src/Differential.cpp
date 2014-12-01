#include "Differential.h"
#include <cmath>

namespace ZGeom {

double calHeatKernel(const EigenSystem& es, int x, int y, double t)
{
    double result(0);
    for (int k = 0; k < es.eigVecCount(); ++k) {
        double *phi = es.getEigVec(k).c_ptr();
        result += std::exp(-es.getEigVal(k)*t) * phi[x] * phi[y];
    }
    return result;
}

double calKernel(const EigenSystem& es, std::function<double(double)> transferFunc, int x, int y)
{
    double result(0);
    for (int k = 0; k < es.eigVecCount(); ++k) {
        double lambda = es.getEigVal(k);
        double *phi = es.getEigVec(k).c_ptr();
        double term = transferFunc(lambda) * phi[x] * phi[y];
        result += term;
    }
    return result;
}

}   // end of namespace