#include <MatchingPursuit.h>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <map>
#include <algorithm>
#include <ppl.h>
#include <concurrent_vector.h>
#include "DenseMatrix.h"

void ZGeom::MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, const InnerProdcutFunc& innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit )
{
	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);
	ZGeom::VecNd vRf = vSignal;

#if 0
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
#else
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
#endif
}

void ZGeom::OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit, MatlabEngineWrapper& engine )
{
	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	const int signalSize = vSignal.size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);

	engine.addColVec(vSignal, "vSignal");

	VecNd vRf = vSignal;
	double *matBasis = new double[nSelected * signalSize];
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
		vPursuit.push_back(std::make_tuple(0, iSelected, maxCoeff));
		availableBasis.erase(iSelected);

		std::copy_n(vBasis[iSelected].c_ptr(), signalSize, matBasis + k * signalSize);
		engine.addArray(matBasis, signalSize, k + 1, false, "matBasis");

		engine.eval("xcoeff = matBasis\\vSignal");
		double *xcoeff = engine.getDblVariablePtr("xcoeff");

		VecNd vNewRf = vSignal;
		for (int b = 0; b <= k; ++b) {
			std::get<2>(vPursuit[b]) = xcoeff[b];
			vNewRf -= xcoeff[b] * vBasis[std::get<1>(vPursuit[b])];
		}
		std::get<0>(vPursuit.back()) = vNewRf.norm2();	
		vRf = vNewRf;
	}

	delete []matBasis;
}
