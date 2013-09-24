#ifndef ZGEOM_MAT_VEC_ARITHMETIC_H
#define ZGEOM_MAT_VEC_ARITHMETIC_H

#include "Vec3.h"
#include "VecN.h"
#include "SparseMatrix.h"
#include "SparseMatVecMultiplier.h"
#include "DenseMatVecMultiplier.h"

namespace ZGeom
{
    template<typename T> 
    inline T cotVec3(const Vec3<T>& v1, const Vec3<T>& v2)
    {
        return ZGeom::dot(v1, v2) / ZGeom::cross(v1, v2).length(); 
    }

    template<typename T> 
    inline T cosVec3(const Vec3<T>& v1, const Vec3<T>& v2)
    {
        return ZGeom::dot(v1, v2) / (v1.length() * v2.length());
    }

    template<typename T> VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

    double innerProduct(const std::vector<double>& v1, const std::vector<double>& v2);

    double innerProductSym(const std::vector<double>& v1, SparseMatVecMultiplier* mulA, const std::vector<double>& v2);

    double innerProductSym(const std::vector<double>& v1, const SparseMatrixCSR<double, int>& A, const std::vector<double>& v2);

    double innerProductSym(const VecNd& v1, const SparseMatrixCSR<double, int>& A, const VecNd& v2);
}

#endif