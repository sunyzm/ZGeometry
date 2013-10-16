#ifndef ZGEOM_SPARSE_MAT_VEC_MULTIPLIER_H
#define ZGEOM_SPARSE_MAT_VEC_MULTIPLIER_H

#include <mkl.h>
#include "MatVecFunctor.h"
#include "SparseMatrix.h"
#include "VecN.h"

namespace ZGeom
{

	class SparseMatVecMultiplier : public MatVecFunctor
	{
	public:
		virtual void operator() (double *in, double* out) { mul(in, out); }   
		SparseMatVecMultiplier(const SparseMatrix<double>& mat, bool isSymmetric = false);
		SparseMatVecMultiplier(int order, int nonzeroCount, int* rows, int *cols, double *vals, bool isSymmetric = false);
		~SparseMatVecMultiplier();

		void mul(double *in, double *out);
		void mul(const VecNd& vin, VecNd& vout);
		VecNd multiply(const VecNd& right);		

	private:
		bool mIsSymmetric;
		MKL_INT* mRowInd;
		MKL_INT* mColInd;
		double* mVal;
		MKL_INT mNonzeroCount;
	};
}

#endif