#ifndef ZGEOM_MAT_VEC_ARITHMETIC_H
#define ZGEOM_MAT_VEC_ARITHMETIC_H

#include "VecN.h"
#include "SparseMatrix.h"
#include "SparseMatVecMultiplier.h"

namespace ZGeom
{
    template<typename T> VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

    double innerProduct(const std::vector<double>& v1, const std::vector<double>& v2);

    double innerProductSym(const std::vector<double>& v1, SparseMatVecMultiplier* mulA, const std::vector<double>& v2);

    double innerProductSym(const std::vector<double>& v1, const SparseMatrixCSR<double, int>& A, const std::vector<double>& v2);

    double innerProductSym(const VecNd& v1, const SparseMatrixCSR<double, int>& A, const VecNd& v2);
}

#endif