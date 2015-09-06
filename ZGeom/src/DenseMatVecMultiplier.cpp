#include "DenseMatVecMultiplier.h"
#include <stdexcept>
#include <mkl.h>

namespace ZGeom
{
	void DenseMatVecMultiplier::mul( const VecNd& vin, VecNd& vout )
	{
		if (vin.size() != mCol) throw std::runtime_error("Mat-Vec multiplier not compatible with input vector");
		vout.resize(mRow, 0.0); 
		mul(vin.c_ptr(), vout.c_ptr());
	}

	void DenseMatVecMultiplier::mul( double *in, double *out )
	{
		char trans = 'N';
		double alpha = 1.0, beta = 0.0;
		int incx = 1, incy = 1;
		cblas_dgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, mRow, mCol, alpha, mMat, mCol, in, incx, beta, out, incy);
	}

	DenseMatVecMultiplier::DenseMatVecMultiplier( const DenseMatrixd& denseMat )
	{
		mRow = denseMat.rowCount();
		mCol = denseMat.colCount();
		mMat = denseMat.raw_ptr();
	}

}
