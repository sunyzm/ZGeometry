#include "VecN.h"
#include "SparseMatrix.h"
#include <cassert>
#include <mkl.h>

namespace ZGeom
{
    template<>
    VecN<double> mulMatVec(const SparseMatrix<double>& mat, const VecN<double>& vec, bool matIsSym)
    {
        assert(mat.colCount() == vec.size() && mat.rowCount() == mat.colCount());
        MKL_INT order = mat.rowCount();
        VecN<double> rv(order);

        if (matIsSym) {
            std::vector<double> vals;
            std::vector<MKL_INT> rowPtr, colIdx;
            mat.convertToCSR(vals, colIdx, rowPtr, MAT_UPPER);
            char uplo = 'U';
            MKL_INT rowCount = mat.rowCount();
            mkl_dcsrsymv(&uplo, &order, &vals[0], &rowPtr[0], &colIdx[0], vec.c_ptr(), rv.c_ptr());
        }
        else {
            std::vector<double> vals;
            std::vector<MKL_INT> rows, cols;
            mat.convertToCOO(rows, cols, vals, MAT_FULL);
            char transa = 'N';
            MKL_INT nnz = vals.size();
            mkl_dcoogemv (&transa, &order, &vals[0], &rows[0], &cols[0], &nnz, vec.c_ptr(), rv.c_ptr());
        }

        return rv;
    }

    template<>
    VecN<float> mulMatVec(const SparseMatrix<float>& mat, const VecN<float>& vec, bool matIsSym)
    {
        assert(mat.colCount() == vec.size() && mat.rowCount() == mat.colCount());
        MKL_INT order = mat.rowCount();
        VecN<float> rv(order);

        if (matIsSym) {
            std::vector<float> vals;
            std::vector<MKL_INT> rowPtr, colIdx;
            mat.convertToCSR(vals, colIdx, rowPtr, MAT_UPPER);
            char uplo = 'U';
            MKL_INT rowCount = mat.rowCount();
            mkl_scsrsymv(&uplo, &order, &vals[0], &rowPtr[0], &colIdx[0], vec.c_ptr(), rv.c_ptr());
        }
        else {
            std::vector<float> vals;
            std::vector<MKL_INT> rows, cols;
            mat.convertToCOO(rows, cols, vals, MAT_FULL);
            char transa = 'N';
            MKL_INT nnz = vals.size();
            mkl_scoogemv (&transa, &order, &vals[0], &rows[0], &cols[0], &nnz, vec.c_ptr(), rv.c_ptr());
        }

        return rv;
    }
}