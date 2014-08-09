#ifndef ZGEOM_SPARSE_SOLVER_H
#define ZGEOM_SPARSE_SOLVER_H
#include "VecN.h"
#include "MatlabEngineWrapper.h"
#include "SparseMatrix.h"


namespace ZGeom
{
	void solveSparse(MatlabEngineWrapper& matlabEngine, const SparseMatrix<double>& matA, const VecNd& vec_b, VecNd& vec_x);
	void solveSparseMultiColumn(MatlabEngineWrapper& matlabEngine, const SparseMatrix<double>& matA, const std::vector<VecNd>& vecB, std::vector<VecNd>& vecX);
}



#endif