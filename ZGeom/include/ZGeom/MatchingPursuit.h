#ifndef ZGEOM_MATCHING_PURSUIT_H
#define ZGEOM_MATCHING_PURSUIT_H
#include <tuple>
#include <functional>
#include <vector>
#include "VecN.h"
#include "Approximation.h"
#include "DenseMatrix.h"
#include "MatlabEngineWrapper.h"

namespace ZGeom
{
	void MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, const InnerProdcutFunc& innerProdFunc, int nSelected, FunctionApproximation& vPursuit );
	
	//// compute OMP with MKL
	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit);
	
	//// compute OMP with Matlab
	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit, MatlabEngineWrapper& engine);
}

#endif