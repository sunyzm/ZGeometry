#include <MatchingPursuit.h>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <map>
#include <algorithm>
#include <ppl.h>
#include <concurrent_vector.h>
#include <amp.h>
#include "DenseMatrix.h"
#include <mkl.h>

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

void ZGeom::OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit)
{
	//std::cout << "For verification: " << innerProdFunc(vBasis[0], vBasis[0]) << ' ' << innerProdFunc(vBasis[1], vBasis[1]) << std::endl;

	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	const int signalSize = vSignal.size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);	

	double *matBasis = new double[nSelected * signalSize];
	double *work = new double[signalSize*17];
	double *a = new double[nSelected * signalSize];
	double *b = new double[signalSize];

	VecNd vRf = vSignal;
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

		double *xcoeff;
		// solve least square via MKL
		{
			int matrix_order = LAPACK_COL_MAJOR;
			char trans = 'N';
			int m = signalSize;
			int n = k + 1;
			int nrhs = 1;
			int lwork = m + 16*m;
			int lda = signalSize;
			int ldb = signalSize;
			std::memcpy(a, matBasis, sizeof(double)*m*n);
			std::memcpy(b, vSignal.c_ptr(), sizeof(double)*m);
			LAPACKE_dgels(matrix_order, trans, m, n, nrhs, a, lda, b, ldb);
			xcoeff = b;
		}

		VecNd vNewRf = vSignal;
		for (int b = 0; b <= k; ++b) {
			std::get<2>(vPursuit[b]) = xcoeff[b];
			vNewRf -= xcoeff[b] * vBasis[std::get<1>(vPursuit[b])];
		}
		std::get<0>(vPursuit.back()) = vNewRf.norm2();	
		vRf = vNewRf;
	}

	delete []work;
	delete []a;
	delete []b;
	delete []matBasis;
}


void ZGeom::OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, std::vector<PursuitApproxItem>& vPursuit, MatlabEngineWrapper& engine )
{
	//std::cout << "For verification: " << innerProdFunc(vBasis[0], vBasis[0]) << ' ' << innerProdFunc(vBasis[1], vBasis[1]) << std::endl;

	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	const int signalSize = vSignal.size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);	


	engine.addColVec(vSignal, "vSignal");
	double *matBasis = new double[nSelected * signalSize];
	VecNd vRf = vSignal;

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
		std::copy_n(vBasis[iSelected].c_ptr(), signalSize, matBasis + k * signalSize);	//append to new selected basis to matBasis matrix

		double *xcoeff;
		// solve least square via matlab
		engine.addArray(matBasis, signalSize, k + 1, false, "matBasis");
		engine.eval("xcoeff = matBasis\\vSignal");
		xcoeff = engine.getDblVariablePtr("xcoeff");

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
