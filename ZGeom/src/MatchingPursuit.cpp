#include <MatchingPursuit.h>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <map>
#include <algorithm>
#include <ppl.h>
#include <concurrent_vector.h>

void ZGeom::MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit )
{
	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
#if 0
	std::set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);
	ZGeom::VecNd vRf = vSignal;
	for (int k = 0; k < nSelected; ++k) {
		double maxCoeff = 0;
		int iSelected = -1;
		for (int iBasis : availableBasis) {
			double candidateCoeff = innerProdFunc(vBasis[iBasis], vRf);
			if (std::fabs(candidateCoeff) > std::fabs(maxCoeff)) {
				maxCoeff = candidateCoeff;
				iSelected = iBasis;
			}
		}
		vRf = vRf - maxCoeff*vBasis[iSelected];
		vPursuit.push_back(std::make_tuple(vRf.norm2(), iSelected, maxCoeff));

		availableBasis.erase(iSelected);
	}
#endif

	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);
	ZGeom::VecNd vRf = vSignal;
	for (int k = 0; k < nSelected; ++k) {
		Concurrency::concurrent_vector<std::pair<int,double> > vCoeff;
		Concurrency::parallel_for_each(availableBasis.begin(), availableBasis.end(), [&](int iBasis) {
			vCoeff.push_back(std::make_pair(iBasis, innerProdFunc(vBasis[iBasis], vRf)));
		});

		double maxCoeff = 0;
		int iSelected = -1;
		for (auto bp : vCoeff) {
			if (std::fabs( bp.second) > std::fabs(maxCoeff)) {
				maxCoeff =  bp.second;
				iSelected = bp.first;
			}
		}
		vRf = vRf - maxCoeff*vBasis[iSelected];
		vPursuit.push_back(std::make_tuple(vRf.norm2(), iSelected, maxCoeff));

		availableBasis.erase(iSelected);
	}
}
