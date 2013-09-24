#ifndef ZGEOM_DENSE_MAT_VEC_MULTIPLIER_H
#define ZGEOM_DENSE_MAT_VEC_MULTIPLIER_H

#include "MatVecFunctor.h"
#include "DenseMatrix.h"
#include "VecN.h"

namespace ZGeom
{

class DenseMatVecMultiplier : public MatVecFunctor
{
public:
	DenseMatVecMultiplier(const DenseMatrixd& denseMat);
	virtual void operator () (double *in, double *out) { mul(in, out); }
	void mul(double *in, double *out);
	void mul(const VecNd& vin, VecNd& vout);
	VecNd multiplier(const VecNd& vin);

private:
	int mRow;
	int mCol;
	double *mMat;
};

}

 

#endif