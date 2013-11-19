#include "MatchingPursuit.h"
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <ppl.h>
#include <concurrent_vector.h>
#include <mkl.h>
#include "tbb/tbb.h"
#include "tbb/concurrent_vector.h"
#include "DenseMatrix.h"
#include <amp.h>

namespace ZGeom
{
	void MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, const InnerProdcutFunc& innerProdFunc, int nSelected, FunctionApproximation& vPursuit )
	{
		if (nSelected <= 0 ||  nSelected > vBasis.size())
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
			vPursuit.addItem(vRf.norm2(), iSelected, maxCoeff);

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
			vPursuit.addItem(vRf.norm2(), iSelected, maxCoeff);

			availableBasis.erase(iSelected);
		}
#endif
	}


	void OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit )
	{
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
			const int nAvaliableBasis = (int)availableBasis.size();
			double maxCoeff = 0;
			int iSelected = -1;

			tbb::concurrent_vector< std::pair<int,double> > vCoeff;
			vCoeff.reserve(availableBasis.size());
			for (int iBasis : availableBasis) vCoeff.push_back(std::make_pair(iBasis, 0));

			tbb::parallel_for_each(vCoeff.begin(), vCoeff.end(), [&](std::pair<int,double>& cp) {
				cp.second = cblas_ddot(signalSize, vBasis[cp.first].c_ptr(), 1, vRf.c_ptr(), 1);
			});		

			for (auto bp : vCoeff) {
				if (std::fabs( bp.second) > std::fabs(maxCoeff)) {
					maxCoeff =  bp.second;
					iSelected = bp.first;
				}
			}
			
			vPursuit.addItem(0, iSelected, maxCoeff);
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
				vPursuit[b].coeff() = xcoeff[b];
				vNewRf -= xcoeff[b] * vBasis[vPursuit[b].index()];
			}
			vPursuit.back().res() = vNewRf.norm2();	
			vRf = vNewRf;
		}

		delete []work;
		delete []a;
		delete []b;
		delete []matBasis;
	}

	void OrthogonalMatchingPursuit_AMP( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, FunctionApproximation& vPursuit )
	{
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

		using namespace Concurrency;
		{
			Concurrency::accelerator default_acc;
			std::wcout << default_acc.device_path << "\n";
			std::wcout << default_acc.dedicated_memory << "\n";
			std::wcout << (default_acc.supports_double_precision ? 
				"double precision: true" : "double precision: false") << "\n";
		}

		//// the following basis matrix should be created and assigned directly on accelerator memory 
		double *matAllBasis = new double[vBasis.size() * signalSize];
		for (size_t i = 0; i < vBasis.size(); ++i)
			std::memcpy((void*)(matAllBasis + i*signalSize), (void*)vBasis[i].c_ptr(), sizeof(double)*signalSize);
		array_view<double, 2> basisView((int)vBasis.size(), signalSize, matAllBasis);
		//////////////////////////////////////////////////////////////////////////

		for (int k = 0; k < nSelected; ++k) {
			const int nAvaliableBasis = (int)availableBasis.size();
			double maxCoeff = 0;
			int iSelected = -1;

			std::vector<int> vBasisIdx(availableBasis.begin(), availableBasis.end());
			array_view<int, 1> basisSelectView(nAvaliableBasis, &vBasisIdx[0]);
			array_view<double, 1> RfView(signalSize, vRf.c_ptr());
			std::vector<double> vBasisCoeff(nAvaliableBasis);
			array_view<double, 1> basisCoeffView(nAvaliableBasis, &vBasisCoeff[0]);

			basisCoeffView.discard_data();
			parallel_for_each(basisCoeffView.extent, [=](index<1> idx) restrict(amp) {
				int selectedBasisIdx = basisSelectView(idx[0]);
				basisCoeffView[idx] = 0;
				for (int j = 0; j < signalSize; ++j) {
					basisCoeffView[idx] += basisView(selectedBasisIdx, j) * RfView(j);
				}
			});
			basisCoeffView.synchronize();

			for (int i = 0; i < nAvaliableBasis; ++i) {
				if (std::fabs(vBasisCoeff[i]) > std::fabs(maxCoeff)) {
					maxCoeff = vBasisCoeff[i];
					iSelected = vBasisIdx[i];
				}
			}

			vPursuit.addItem(0, iSelected, maxCoeff);
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
				vPursuit[b].coeff() = xcoeff[b];
				vNewRf -= xcoeff[b] * vBasis[vPursuit[b].index()];
			}
			vPursuit.back().res() = vNewRf.norm2();	
			vRf = vNewRf;
		}

		delete []matAllBasis;
		delete []work;
		delete []a;
		delete []b;
		delete []matBasis;
	}

	
	void OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit)
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
			tbb::concurrent_vector< std::pair<int,double> > vCoeff;
			vCoeff.reserve(availableBasis.size());

			for (int iBasis : availableBasis) vCoeff.push_back(std::make_pair(iBasis, 0));
			tbb::parallel_for_each(vCoeff.begin(), vCoeff.end(), [&](std::pair<int,double>& cp) {
				cp.second = innerProdFunc(vBasis[cp.first], vRf);
			});			
			
			/*tbb::parallel_for_each(availableBasis.begin(), availableBasis.end(), [&](int iBasis) {
				vCoeff.push_back(std::make_pair(iBasis, innerProdFunc(vBasis[iBasis], vRf)));
			});*/


			double maxCoeff = 0;
			int iSelected = -1;
			for (auto bp : vCoeff) {
				if (std::fabs( bp.second) > std::fabs(maxCoeff)) {
					maxCoeff =  bp.second;
					iSelected = bp.first;
				}
			}
			vPursuit.addItem(0, iSelected, maxCoeff);
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
				vPursuit[b].coeff() = xcoeff[b];
				vNewRf -= xcoeff[b] * vBasis[vPursuit[b].index()];
			}
			vPursuit.back().res() = vNewRf.norm2();	
			vRf = vNewRf;
		}

		delete []work;
		delete []a;
		delete []b;
		delete []matBasis;
	}

	void OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, FunctionApproximation& vPursuit, MatlabEngineWrapper& engine )
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
			tbb::concurrent_vector<std::pair<int,double> > vCoeff;
			tbb::parallel_for_each(availableBasis.begin(), availableBasis.end(), [&](int iBasis) {
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
			vPursuit.addItem(0, iSelected, maxCoeff);
			availableBasis.erase(iSelected);
			std::copy_n(vBasis[iSelected].c_ptr(), signalSize, matBasis + k * signalSize);	//append to new selected basis to matBasis matrix

			double *xcoeff;
			// solve least square via matlab
			engine.addArray(matBasis, signalSize, k + 1, false, "matBasis");
			engine.eval("xcoeff = matBasis\\vSignal");
			xcoeff = engine.getDblVariablePtr("xcoeff");

			VecNd vNewRf = vSignal;
			for (int b = 0; b <= k; ++b) {
				vPursuit[b].coeff() = xcoeff[b];
				vNewRf -= xcoeff[b] * vBasis[vPursuit[b].index()];
			}
			vPursuit.back().res() = vNewRf.norm2();	
			vRf = vNewRf;
		}

		delete []matBasis;
	}

} // end of namespace
