#include "VecN.h"
#define USE_MKL

namespace ZGeom
{
#ifdef USE_MKL
#include <mkl.h>
    /* specialization for operator *= */
    template<>
    const VecN<double>& VecN<double>::operator *= (double scale)
    {
        cblas_dscal(mDim, scale, mVec, 1);
        return *this;
    }
    template<>
    const VecN<float>& VecN<float>::operator *= (float scale)
    {
        cblas_sscal(mDim, scale, mVec, 1);
        return *this;
    }

    /* specialization for dot */
    template<>
    float VecN<float>::dot(const VecN<float>& v2) const
    {
        assert(mDim == v2.mDim);
        return cblas_sdot(mDim, mVec, 1, v2.mVec, 1);
    }
    template<>
    double VecN<double>::dot(const VecN<double>& v2) const
    {
        assert(mDim == v2.mDim);
        return cblas_ddot(mDim, mVec, 1, v2.mVec, 1);
    }

#endif
}