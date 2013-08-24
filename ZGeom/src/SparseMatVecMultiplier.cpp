#include "SparseMatVecMultiplier.h"
#include <algorithm>
#include <mkl.h>

namespace ZGeom
{

    SparseMatVecMultiplier::SparseMatVecMultiplier(int order, int nonzeroCount, int* rows, int *cols, double *vals, bool isSymmetric)
    {
        mIsSymmetric = isSymmetric;
        mOrder = order;

        mRowInd = new MKL_INT[nonzeroCount];
        mColInd = new MKL_INT[nonzeroCount];
        mVal    = new double[nonzeroCount];

        mNonzeroCount = 0;
        for (int i = 0; i < nonzeroCount; ++i) {
            if (isSymmetric && rows[i] > cols[i]) continue;
            mRowInd[mNonzeroCount] = rows[i];
            mColInd[mNonzeroCount] = cols[i];
            mVal[mNonzeroCount]    = vals[i];
            mNonzeroCount++;
        }
    }

    SparseMatVecMultiplier::SparseMatVecMultiplier(const SparseMatrix<double>& mat, bool isSymmetric)
    {
        mIsSymmetric = isSymmetric;
        mOrder = mat.rowCount();

        std::vector<MKL_INT> vRowInd, vColInd;
        std::vector<double> vVal;
        MatrixForm matForm = (isSymmetric ? MAT_UPPER : MAT_FULL);
        mat.convertToCOO(vRowInd, vColInd, vVal, matForm);

        mNonzeroCount = vVal.size();
        mRowInd = new MKL_INT[mNonzeroCount];
        mColInd = new MKL_INT[mNonzeroCount];
        mVal    = new double[mNonzeroCount];

        for (int k = 0; k < mNonzeroCount; ++k ) {
            mRowInd[k] = vRowInd[k];
            mColInd[k] = vColInd[k];
            mVal[k]    = vVal[k];
        }
    }

    SparseMatVecMultiplier::~SparseMatVecMultiplier()
    {
        delete []mRowInd;
        delete []mColInd;
        delete []mVal;
    }

    void SparseMatVecMultiplier::operator() (double *in, double* out)
    {
        MKL_INT order = mOrder;

        if (mIsSymmetric) {
            char uplo = 'U';
            mkl_dcoosymv(&uplo, &order, mVal, mRowInd, mColInd, &mNonzeroCount, in, out);
        } else {
            char transa = 'N';
            mkl_dcoogemv(&transa, &order, mVal, mRowInd, mColInd, &mNonzeroCount, in, out);
        }
    }

} // end of namespace mm_laplace

