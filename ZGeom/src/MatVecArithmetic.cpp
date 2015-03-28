#include "MatVecArithmetic.h"
#include <cassert>
#include <iostream>
#include <mkl.h>

namespace ZGeom {

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

double innerProductStd( const std::vector<double>& v1, const std::vector<double>& v2 )
{
	assert(v1.size() == v2.size());
	return cblas_ddot((int)v1.size(), &v1[0], 1, &v2[0], 1);
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

double RegularProductFunc(const VecN<double>& v1, const VecN<double>& v2)
{
	assert(v1.size() == v2.size());
	return cblas_ddot(v1.size(), v1.c_ptr(), 1, v2.c_ptr(), 1);
}

void mulMatMat( const SparseMatrix<double>& mat1, const SparseMatrix<double>& mat2, SparseMatrix<double>& mat3 )
{
	if (mat1.colCount() != mat2.rowCount()) 
		throw std::runtime_error("Matrix dimension not compatible for multiplication!");

	std::vector<int> ja, ia, jb, ib;
	std::vector<double> a, b;
	mat1.convertToCSR(a, ja, ia, MAT_FULL);
	mat2.convertToCSR(b, jb, ib, MAT_FULL);

	char trans = 'N';
	int sort = 0;
	int m = mat1.rowCount();
	int n = mat1.colCount();
	int k = mat2.colCount();
	int nzmax = 0;
	int info = 0;
	double *c = NULL;
	int *jc = NULL;
	int *ic = new int[m+1];
	int request = 1;

	mkl_dcsrmultcsr(&trans, &request, &sort, &m, &n, &k, &a[0], &ja[0], &ia[0],
					&b[0], &jb[0], &ib[0], c, jc, ic, &nzmax, &info);
	if (info != 0) {
		std::cout << "mkl_dcsrmultcsr error code: " << info;
		delete []ic;
		return;
	}

	nzmax = ic[m] - 1;
	c = new double[nzmax];
	jc = new int[nzmax];
	request = 2;
		
	mkl_dcsrmultcsr(&trans, &request, &sort, &m, &n, &k, &a[0], &ja[0], &ia[0],
					&b[0], &jb[0], &ib[0], c, jc, ic, &nzmax, &info);
	
    if (info != 0) std::cout << "mkl_dcsrmultcsr error code: " << info;	
    else mat3.convertFromCSR(m, k, c, jc, ic);		

	delete []c;
	delete []jc;
	delete []ic;
}

void addMatMat( const SparseMatrix<double>& mat1, const SparseMatrix<double>& mat2, double beta, SparseMatrix<double>& mat3 )
{
	if (mat1.rowCount() != mat2.rowCount() || mat1.colCount() != mat2.colCount()) 
		throw std::runtime_error("Matrix dimension not compatible for addition!");

	char trans = 'N';
	int request = 1;
	int sort = 0;
	int m = (int)mat1.rowCount();
	int n = (int)mat1.colCount();
	std::vector<double> a, b;
	std::vector<int> ja, ia, jb, ib;
	mat1.convertToCSR(a, ja, ia, MAT_FULL);
	mat2.convertToCSR(b, jb, ib, MAT_FULL);
	int nzmax = 0;
	double *c = NULL;
	int *jc = NULL; 
	int *ic = new int[m+1];
	int info = 0;

	mkl_dcsradd(&trans, &request, &sort, &m, &n, &a[0], &ja[0], &ia[0], &beta,
			    &b[0], &jb[0], &ib[0], c, jc, ic, &nzmax, &info);
	if (info != 0) {
		std::cout << "mkl_dcsrmultcsr error code: " << info;
		delete []ic;
		return;
	}

	nzmax = ic[m] - 1;
	c = new double[nzmax];
	jc = new int[nzmax];
	request = 2;

	mkl_dcsradd(&trans, &request, &sort, &m, &n, &a[0], &ja[0], &ia[0], &beta,
				&b[0], &jb[0], &ib[0], c, jc, ic, &nzmax, &info);
	if (info != 0) {
		std::cout << "mkl_dcsrmultcsr error code: " << info;
		delete []ic;
		delete []jc;
		delete []c;
		return;
	}

	mat3.convertFromCSR(m, n, c, jc, ic);		
	delete []c;
	delete []jc;
	delete []ic;
}

} // end of namespace