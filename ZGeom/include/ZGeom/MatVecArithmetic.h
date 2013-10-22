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
	VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

	void mulMatMat(const SparseMatrix<double>& mat1, const SparseMatrix<double>& mat2, SparseMatrix<double>& mat3);
	void addMatMat(const SparseMatrix<double>& mat1, const SparseMatrix<double>& mat2, double beta, SparseMatrix<double>& mat3);	// compute mat3 = mat1 + beta * mat2
	
	double innerProduct(const std::vector<double>& v1, const std::vector<double>& v2);

	double innerProductSym(const std::vector<double>& v1, SparseMatVecMultiplier* mulA, const std::vector<double>& v2);

	double innerProductSym(const std::vector<double>& v1, const SparseMatrixCSR<double, int>& A, const std::vector<double>& v2);

	double innerProductSym(const VecNd& v1, const SparseMatrixCSR<double, int>& A, const VecNd& v2);
}

#endif