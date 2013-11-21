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
	
	//// compute OMP with MKL and regular inner product
	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit);

	//// compute Multi-channel OMP with MKL and regular inner product
	void SimultaneousOMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuits);

	//// compute OMP with MKL and regular inner product
	void OrthogonalMatchingPursuit_AMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit);

	//// compute OMP with MKL and customized InnerProdcutFunc
	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit);
	
	//// compute OMP with Matlab and customized InnerProductFunc
	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit, MatlabEngineWrapper& engine);
}

#endif