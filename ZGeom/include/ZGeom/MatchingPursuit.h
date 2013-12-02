#ifndef ZGEOM_MATCHING_PURSUIT_H
#define ZGEOM_MATCHING_PURSUIT_H
#include <tuple>
#include <functional>
#include <vector>
#include <mkl.h>
#include "VecN.h"
#include "Approximation.h"
#include "DenseMatrix.h"
#include "MatlabEngineWrapper.h"

namespace ZGeom
{
	const InnerProdcutFunc RegularProductFunc = 
		[](const VecN<double>& v1, const VecN<double>& v2) 
	{
		assert(v1.size() == v2.size());
		return cblas_ddot(v1.size(), v1.c_ptr(), 1, v2.c_ptr(), 1);
	};

	void GeneralizedSimultaneousFourierApprox(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuits, const InnerProdcutFunc& innerProdFunc = RegularProductFunc);

	void MatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit);	
	void GeneralizedMP( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit, const InnerProdcutFunc& innerProdFunc);

	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit);
	void OMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit);
	void GeneralizedOMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit, InnerProdcutFunc innerProdFunc);

	void SimultaneousMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuit, double p = 2.);
	void GeneralizedSimultaneousMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuit, const InnerProdcutFunc& innerProdFunc, double p = 2);
	
	//// compute Multi-channel OMP with MKL and regular inner product
	void SimultaneousOMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuits, double p = 2.);
	void GeneralizedSimultaneousOMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuits, const InnerProdcutFunc& innerProdFunc, double p = 2.);

	//// compute OMP with MKL and regular inner product
	void OrthogonalMatchingPursuit_AMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit);
		
	//// compute OMP with Matlab and customized InnerProductFunc
	void GeneralizedOMP_MATLAB(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit, MatlabEngineWrapper& engine);
}

#endif