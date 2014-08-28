#ifndef ZGEOM_SPARSE_SYMMETRIC_MAT_VEC_SOLVER_H
#define ZGEOM_SPARSE_SYMMETRIC_MAT_VEC_SOLVER_H
#include <mkl.h>
#include "SparseSymSolver.h"
#include "SparseMatrix.h"
#include "SparseMatrixCSR.h"
#include "MatVecFunctor.h"

namespace ZGeom {

class SparseSymMatVecSolver : public SparseSymSolver<double>, public MatVecFunctor
{
public:
	SparseSymMatVecSolver() {}
	SparseSymMatVecSolver(const SparseMatrix<double>& laplacianMat, bool isPositiveDefinite, bool verbose = false)
		: SparseSymSolver<double>(laplacianMat, isPositiveDefinite, verbose) {}

	virtual void operator() (double *in, double* out);
};

inline void SparseSymMatVecSolver::operator () (double *in, double* out)
{
	solve(1, in, out);
}

}   // end of namespace

#endif