#include "MatVecArithmetic.h"
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
            MKL_INT nnz = (int)vals.size();
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
            MKL_INT nnz = (int)vals.size();
            mkl_scoogemv (&transa, &order, &vals[0], &rows[0], &cols[0], &nnz, vec.c_ptr(), rv.c_ptr());
        }

        return rv;
    }

    double innerProduct(const std::vector<double>& v1, const std::vector<double>& v2)
    {
        assert(v1.size() == v2.size());

        MKL_INT n = (int)v1.size();
        MKL_INT xinc = 1, yinc = 1;

        return ddot(&n, &v1[0], &xinc, &v2[0], &yinc);
    }

    double innerProductSym(const std::vector<double>& v1, SparseMatVecMultiplier* mulA, const std::vector<double>& v2)
    {
        assert(v1.size() == mulA->getOrder() && v1.size() == v2.size());
        MKL_INT m = (int)v1.size();
        double *pv1 = const_cast<double*>(&v1[0]);
        double *pv2 = const_cast<double*>(&v2[0]);

        std::vector<double> vy(m);
        (*mulA)(pv2, &vy[0]);

        MKL_INT xinc = 1, yinc = 1;    
        return ddot(&m, pv1, &xinc, &vy[0], &yinc);
    }

    double innerProductSym(const std::vector<double>& v1, const SparseMatrixCSR<double, MKL_INT>& A, const std::vector<double>& v2)
    {
        assert(v1.size() == A.rowCount() && A.rowCount() == v2.size());

        MKL_INT m = (int)v1.size();
        MKL_INT *ja = A.colIdx();
        MKL_INT *ia = A.rowPtr();
        double *a = A.nzVal();
        char *uplo = "U";
        double *pv1 = const_cast<double*>(&v1[0]);
        std::vector<double> vy(m);

        mkl_dcsrsymv(uplo, &m, a, ia, ja, pv1, &vy[0]);

        MKL_INT xinc = 1, yinc = 1;
        double *pv2 = const_cast<double*>(&v2[0]);

        return ddot(&m, &vy[0], &xinc, pv2, &yinc);
    }

    double innerProductSym( const VecNd& v1, const SparseMatrixCSR<double, int>& A, const VecNd& v2 )
    {
        assert(v1.size() == A.rowCount() && A.rowCount() == v2.size());

        MKL_INT m = v1.size();
        MKL_INT *ja = A.colIdx();
        MKL_INT *ia = A.rowPtr();
        double *a = A.nzVal();
        char *uplo = "U";
        double *pv1 = v1.c_ptr();
        std::vector<double> vy(m);

        mkl_dcsrsymv(uplo, &m, a, ia, ja, pv1, &vy[0]);

        MKL_INT xinc = 1, yinc = 1;
        double *pv2 = v2.c_ptr();

        return ddot(&m, &vy[0], &xinc, pv2, &yinc);
    }

}