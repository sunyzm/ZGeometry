#ifndef ZGEOM_MATCHING_PURSUIT_H
#define ZGEOM_MATCHING_PURSUIT_H
#include <tuple>
#include <functional>
#include "VecN.h"
#include "MatlabEngineWrapper.h"

namespace ZGeom
{
	typedef std::tuple<double, int, double> PursuitApproxItem;

	void MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, const InnerProdcutFunc& innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit );
	void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit, MatlabEngineWrapper& engine);
}

#endif