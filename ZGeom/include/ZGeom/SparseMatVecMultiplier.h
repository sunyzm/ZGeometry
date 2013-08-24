#ifndef ZGEOM_SPARSE_MAT_VEC_MULTIPLIER_H
#define ZGEOM_SPARSE_MAT_VEC_MULTIPLIER_H

#include <mkl.h>
#include "MatVecFunctor.h"
#include "SparseMatrix.h"

namespace ZGeom
{

    class SparseMatVecMultiplier : public MatVecFunctor
    {
    public:
        virtual void operator() (double *in, double* out);    
        SparseMatVecMultiplier(const SparseMatrix<double>& mat, bool isSymmetric = false);
        SparseMatVecMultiplier(int order, int nonzeroCount, int* rows, int *cols, double *vals, bool isSymmetric = false);
        ~SparseMatVecMultiplier();

    private:
        bool mIsSymmetric;
        MKL_INT* mRowInd;
        MKL_INT* mColInd;
        double* mVal;
        MKL_INT mNonzeroCount;
    };
}

#endif