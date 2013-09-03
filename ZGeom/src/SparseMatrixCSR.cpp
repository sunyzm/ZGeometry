#include "SparseMatrixCSR.h"
#include <cassert>
#include <mkl.h>

namespace ZGeom
{
    template<>
    void ZGeom::addCSR<double,MKL_INT>(const SparseMatrixCSR<double, MKL_INT>& m1, const SparseMatrixCSR<double, MKL_INT>& m2, double beta, SparseMatrixCSR<double, MKL_INT>& m3)
    {
        assert(m1.rowCount() == m2.rowCount());    

        MKL_INT rowCount = m1.rowCount();
        MKL_INT nzmax = m1.nonzeroCount() + m2.nonzeroCount();       //maximum count of non-zero elements

        m3.mVal = new double[nzmax];
        m3.mColIdx = new MKL_INT[nzmax];
        m3.mRowPtr = new MKL_INT[rowCount+1]; 
        m3.mRowCount = rowCount;

        char trans = 'N';
        MKL_INT request = 0;
        MKL_INT sort = 3;
        MKL_INT m = rowCount, n = rowCount;
        MKL_INT info;

        mkl_dcsradd(&trans, &request, &sort, &m, &n, m1.nzVal(), m1.colIdx(), m1.rowPtr(), &beta, 
            m2.nzVal(), m2.colIdx(), m2.rowPtr(), m3.mVal, m3.mColIdx, m3.mRowPtr, &nzmax, &info); 

        m3.mNonzeroCount = m3.mRowPtr[rowCount] - 1;
    }

    template<>
    void ZGeom::addCSR<float,MKL_INT>(const SparseMatrixCSR<float, MKL_INT>& m1, const SparseMatrixCSR<float, MKL_INT>& m2, float beta, SparseMatrixCSR<float, MKL_INT>& m3)
    {
        assert(m1.rowCount() == m2.rowCount());    

        MKL_INT rowCount = m1.rowCount();
        MKL_INT nzmax = m1.nonzeroCount() + m2.nonzeroCount();       //maximum count of non-zero elements

        m3.mVal = new float[nzmax];
        m3.mColIdx = new MKL_INT[nzmax];
        m3.mRowPtr = new MKL_INT[rowCount+1];  
        m3.mRowCount = rowCount; 

        char trans = 'N';
        MKL_INT request = 0;
        MKL_INT sort = 3;
        MKL_INT m = rowCount, n = rowCount;
        MKL_INT info;

        mkl_scsradd(&trans, &request, &sort, &m, &n, m1.nzVal(), m1.colIdx(), m1.rowPtr(), &beta, 
            m2.nzVal(), m2.colIdx(), m2.rowPtr(), m3.mVal, m3.mColIdx, m3.mRowPtr, &nzmax, &info);    

        m3.mNonzeroCount = m3.mRowPtr[rowCount] - 1;
    }
}