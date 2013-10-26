#ifndef ZGEOM_MATCHING_PURSUIT_H
#define ZGEOM_MATCHING_PURSUIT_H
#include <tuple>
#include <functional>
#include "VecN.h"

namespace ZGeom
{
	typedef std::tuple<double, int, double> PursuitApproxItem;
	typedef std::function<double(const VecNd&, const VecNd&)> InnerProdcutFunc;

	void MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelectedBasis, std::vector<PursuitApproxItem>& vPursuit );

}

#endif