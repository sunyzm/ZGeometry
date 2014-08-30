#include "VecN.h"
#include <numeric>
#include <mkl.h>

template<>
inline double ZGeom::VecN<double>::dot(const VecN<double>& v2) const
{
    assert(mDim == v2.mDim);
    return cblas_ddot(mDim, mVec, 1, v2.mVec, 1);
    //return std::inner_product(mVec, mVec + mDim, v2.mVec, 0.);
}