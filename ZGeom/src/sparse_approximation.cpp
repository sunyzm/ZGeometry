#include "sparse_approximation.h"
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <ppl.h>
#include <amp.h>
#include <concurrent_vector.h>
#include <mkl.h>
#include "DenseMatrix.h"

namespace ZGeom {

void GeneralizedMultiChannelFourierApprox(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, InnerProdcutFunc innerProdFunc)
{
	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	assert(vSignals.size() == vPursuits.size());
	for (auto p : vPursuits) p->clear();

	const int nChannels = (int)vSignals.size();
	const int signalSize = vSignals[0].size();
	std::vector<VecNd> vRfs = vSignals;

	for (int c = 0; c < nChannels; ++c) {
		for (int k = 0; k < nSelected; ++k) {
			double coeff = innerProdFunc(vBasis[k], vSignals[c]);
			vRfs[c] -= coeff * vBasis[k];
			vPursuits[c]->addItem(k, coeff);
		}
	}
}


void MatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit )
{
	GeneralizedMP(vSignal, vBasis, nSelected, vPursuit, RegularProductFunc);
}


void GeneralizedMP( const VecNd& vSignal, const std::vector<VecNd>& vAtoms, int nSelected, SparseCoding& vPursuit, const InnerProdcutFunc& innerProdFunc)
{
	using Concurrency::concurrent_vector;
	using Concurrency::parallel_for;

	if (nSelected <= 0 ||  nSelected > vAtoms.size())
		throw std::logic_error("nSelected (basis) too small or too large!");
	vPursuit.clear();

	int dictSize = (int)vAtoms.size();
	ZGeom::VecNd vRf = vSignal;
	concurrent_vector<double> vCoeff(dictSize);

	for (int k = 0; k < nSelected; ++k) {
		parallel_for(0, dictSize, [&](int iBasis) {
			vCoeff[iBasis] = std::fabs(innerProdFunc(vAtoms[iBasis], vRf));
		});

		int iSelected = int(std::max_element(vCoeff.begin(), vCoeff.end()) - vCoeff.begin());

		double coeffSelected = innerProdFunc(vAtoms[iSelected], vRf);
		vRf -= coeffSelected * vAtoms[iSelected];
		vPursuit.addItem(iSelected, coeffSelected);
	}
}


void OrthogonalMatchingPursuit( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit )
{
	using namespace Concurrency;

	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	const int signalSize = vSignal.size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);	

	double *matBasis = new double[nSelected * signalSize];
	double *a = new double[nSelected * signalSize];
	double *b = new double[signalSize];

	VecNd vRf = vSignal;

	for (int k = 0; k < nSelected; ++k) {
		const int nAvaliableBasis = (int)availableBasis.size();
		double maxCoeff = 0;
		int iSelected = -1;

		concurrent_vector<std::pair<int,double> > vCoeff;
		vCoeff.reserve(availableBasis.size());
		for (int iBasis : availableBasis) vCoeff.push_back(std::make_pair(iBasis, 0));

		parallel_for_each(vCoeff.begin(), vCoeff.end(), [&](std::pair<int,double>& cp) {
			cp.second = cblas_ddot(signalSize, vBasis[cp.first].c_ptr(), 1, vRf.c_ptr(), 1);
		});		

		for (auto bp : vCoeff) {
			if (std::fabs( bp.second) > std::fabs(maxCoeff)) {
				maxCoeff =  bp.second;
				iSelected = bp.first;
			}
		}
			
		vPursuit.addItem(iSelected, maxCoeff);
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
			int lda = signalSize;
			int ldb = signalSize;
			std::memcpy(a, matBasis, sizeof(double)*m*n);
			std::memcpy(b, vSignal.c_ptr(), sizeof(double)*m);
			LAPACKE_dgels(matrix_order, trans, m, n, nrhs, a, lda, b, ldb);
			xcoeff = b;
		}

		VecNd vNewRf = vSignal;
		for (int j = 0; j <= k; ++j) {
			vPursuit[j].coeff() = xcoeff[j];
			vNewRf -= xcoeff[j] * vBasis[vPursuit[j].index()];
		}
		vRf = vNewRf;
	}

	delete []a;
	delete []b;
	delete []matBasis;
}


void OMP( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit )
{
	GeneralizedOMP(vSignal, vBasis, nSelected, vPursuit, RegularProductFunc);
}


void GeneralizedOMP( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit, InnerProdcutFunc innerProdFunc )
{
	using Concurrency::parallel_for_each;
	using Concurrency::parallel_for;
	using Concurrency::concurrent_vector;

	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	const int signalSize = vSignal.size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);	

	double *matBasis = new double[signalSize * nSelected];
	double *a = new double[signalSize * nSelected];
	double *b = new double[signalSize];
	VecNd vRf = vSignal;

	for (int k = 0; k < nSelected; ++k) 
	{
		int nAvailableBasis = (int)availableBasis.size();
		concurrent_vector<int> vAvailableBasis(availableBasis.begin(), availableBasis.end());
		concurrent_vector<double> vCoeff(nAvailableBasis);
		parallel_for(0, nAvailableBasis, [&](int iAB){
			vCoeff[iAB] = std::fabs(innerProdFunc(vBasis[vAvailableBasis[iAB]], vRf));
		});
		int iABSelected = int(std::max_element(vCoeff.begin(), vCoeff.end()) - vCoeff.begin());
		int iSelected = vAvailableBasis[iABSelected];
		double maxCoeff = innerProdFunc(vBasis[iSelected], vRf);

		vPursuit.addItem(iSelected, maxCoeff);
		availableBasis.erase(iSelected);
		std::copy_n(vBasis[iSelected].c_ptr(), signalSize, matBasis + k * signalSize);

		// solve least square via MKL
		{
			int matrix_order = LAPACK_COL_MAJOR;
			char trans = 'N';
			int m = signalSize;
			int n = k + 1;
			int nrhs = 1;
			int lda = signalSize;
			int ldb = signalSize;
			std::memcpy(a, matBasis, sizeof(double)*m*n);	//copy to 'a' necessary because 'a' will be overwritten as output of dgels
			std::memcpy(b, vSignal.c_ptr(), sizeof(double)*m);
			LAPACKE_dgels(matrix_order, trans, m, n, nrhs, a, lda, b, ldb);
		}

		VecNd vNewRf = vSignal;
		for (int j = 0; j <= k; ++j) {
			vPursuit[j].coeff() = b[j];
			vNewRf -= b[j] * vBasis[vPursuit[j].index()];
		}
		vRf = vNewRf;
	}

	delete []a;
	delete []b;
	delete []matBasis;
}


void SimultaneousMP( const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, double p /*= 2.*/ )
{
	GeneralizedSimultaneousMP(vSignals, vBasis, nSelected, vPursuits, RegularProductFunc, p);
}


void GeneralizedSimultaneousMP( const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vAtoms, 
		                        int nSelected, std::vector<SparseCoding*>& vPursuits, 
								const InnerProdcutFunc& innerProdFunc, double p /*= 1*/ )
{
	using Concurrency::concurrent_vector;
	using Concurrency::parallel_for;

	if (nSelected <= 0 || nSelected > vAtoms.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	assert(vSignals.size() == vPursuits.size());
				
	const int nChannels = (int)vSignals.size();
	const int signalSize = vSignals[0].size();
	const int dictSize = (int)vAtoms.size();

	for (auto pur : vPursuits) pur->clear();
	std::vector<VecNd> vRfs = vSignals;
	concurrent_vector<double> vCoeffNorm(dictSize);

	for (int k = 0; k < nSelected; ++k) {
		parallel_for(0, dictSize, [&](int atomIdx) {
			VecNd channelCoeff(nChannels);
			for (int c = 0; c < nChannels; ++c)
				channelCoeff[c] = innerProdFunc(vAtoms[atomIdx], vRfs[c]);
			vCoeffNorm[atomIdx] = channelCoeff.pNorm(p);
		});		

		int iSelected = int(std::max_element(vCoeffNorm.begin(), vCoeffNorm.end()) - vCoeffNorm.begin());
			
		for (int c = 0; c < nChannels; ++c) {
			double coeff = innerProdFunc(vAtoms[iSelected], vRfs[c]);
			vRfs[c] -= coeff * vAtoms[iSelected];
			vPursuits[c]->addItem(iSelected, coeff);
		}
	}
}


void SimultaneousOMP( const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, double p /*= 1*/ )
{
	GeneralizedSimultaneousOMP(vSignals, vBasis, nSelected, vPursuits, RegularProductFunc, p);
}


void GeneralizedSimultaneousOMP( const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, const InnerProdcutFunc& innerProdFunc, double p /*= 2.*/ )
{
	using namespace Concurrency;

	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	assert(vSignals.size() == vPursuits.size());
	int nChannels = (int)vSignals.size();

	for (auto ptr : vPursuits) ptr->clear();

	const int signalSize = vSignals[0].size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);	

	std::vector<VecNd> vRfs = vSignals;
	double *signalData = new double[signalSize * nChannels];
	for (int c = 0; c < nChannels; ++c)
		std::memcpy((void*)(signalData + signalSize * c), (void*)vSignals[c].c_ptr(), sizeof(double)*signalSize);

	double *matBasis = new double[nSelected * signalSize];
	double *a = new double[nSelected * signalSize];
	double *b = new double[signalSize * nChannels];		

	for (int k = 0; k < nSelected; ++k) {
		const int nAvaliableBasis = (int)availableBasis.size();
		double maxCoeff = 0;
		int iSelected = -1;

		concurrent_vector< std::pair<int,double> > vCoeff;
		vCoeff.reserve(availableBasis.size());
		for (int iBasis : availableBasis) vCoeff.push_back(std::make_pair(iBasis, 0));

		parallel_for_each(vCoeff.begin(), vCoeff.end(), [&](std::pair<int,double>& cp) {
			VecNd channelCoeff(nChannels);
			for (int c = 0; c < nChannels; ++c)
				channelCoeff[c] = innerProdFunc(vBasis[cp.first], vRfs[c]);

			cp.second = channelCoeff.pNorm(p);
		});		

		for (auto& bp : vCoeff) {
			if (bp.second > maxCoeff) {
				maxCoeff =  bp.second;
				iSelected = bp.first;
			}
		}

		for (auto v : vPursuits) v->addItem(iSelected, 0);
		availableBasis.erase(iSelected);
		std::copy_n(vBasis[iSelected].c_ptr(), signalSize, matBasis + k * signalSize);

		// solve least squares via MKL
		{
			int matrix_order = LAPACK_COL_MAJOR;
			char trans = 'N';
			int m = signalSize;
			int n = k + 1;
			int nrhs = nChannels;
			int lda = signalSize;
			int ldb = signalSize;
			std::memcpy(a, matBasis, sizeof(double)*m*n);
			std::memcpy(b, signalData, sizeof(double)*m*nChannels);
			LAPACKE_dgels(matrix_order, trans, m, n, nrhs, a, lda, b, ldb);

			for (int c = 0; c < nChannels; ++c) {
				VecNd vNewRf = vSignals[c];
				for (int j = 0; j <= k; ++j) {
					(*vPursuits[c])[j].coeff() = b[c*ldb + j];
					vNewRf -= (*vPursuits[c])[j].coeff() * vBasis[(*vPursuits[c])[j].index()];
				}

				vRfs[c] = vNewRf;
			}
		}
	}

	delete []a;
	delete []b;
	delete []matBasis;
	delete []signalData;
}
	

void OrthogonalMatchingPursuit_AMP( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit )
{
	if (nSelected <= 0 || nSelected > vBasis.size())
		throw std::logic_error("nSelectedBasis too small or too large!");
	vPursuit.clear();
	const int signalSize = vSignal.size();
	std::unordered_set<int> availableBasis;
	for (int i = 0; i < vBasis.size(); ++i) availableBasis.insert(i);	

	double *matBasis = new double[nSelected * signalSize];
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

		vPursuit.addItem(iSelected, maxCoeff);
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
		vRf = vNewRf;
	}

	delete []matAllBasis;
	delete []a;
	delete []b;
	delete []matBasis;
}
	
void GeneralizedOMP_MATLAB( const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, SparseCoding& vPursuit, MatlabEngineWrapper& engine )
{
	using namespace Concurrency;

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
		concurrent_vector<std::pair<int,double> > vCoeff;
		parallel_for_each(availableBasis.begin(), availableBasis.end(), [&](int iBasis) {
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
		vPursuit.addItem(iSelected, maxCoeff);
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
		vRf = vNewRf;
	}

	delete []matBasis;
}

void OMP_SPAMS(MatlabEngineWrapper& engine, const VecNd& vSignal, const std::vector<VecNd>& vAtoms, int supportSize, SparseCoding& fa )
{
	int signalSize = (int)vSignal.size();
	int atomCount = (int)vAtoms.size();

	engine.addColVec(vSignal, "X");
	double *pMatDict = new double[signalSize*atomCount];
	for (int i = 0; i < atomCount; ++i) {
		std::copy_n(vAtoms[i].c_ptr(), signalSize, pMatDict + signalSize*i);
	}
	engine.addArray(pMatDict, signalSize, atomCount, false, "D");
	engine.addDoubleScalar((double)supportSize, "L0");
	engine.eval("[idx,coeff]=spamsOMP(X,D,L0);");

	double *pIdx = engine.getDblVariablePtr("idx");
	double *pCoeff = engine.getDblVariablePtr("coeff");
	fa.clear();
	for (int i = 0; i < supportSize; ++i)
		fa.addItem((int)pIdx[i] - 1, pCoeff[i]);

	engine.removeVariable("X");
	engine.removeVariable("D");
	engine.removeVariable("L0");
	delete []pMatDict;
}

void LASSO_SPAMS(MatlabEngineWrapper& engine, const VecNd& vSignal, const std::vector<VecNd>& vAtoms, double lambda, SparseCoding& fa )
{
	int signalSize = (int)vSignal.size();
	int atomCount = (int)vAtoms.size();

	engine.addColVec(vSignal, "X");
	double *pMatDict = new double[signalSize*atomCount];
	for (int i = 0; i < atomCount; ++i) {
		std::copy_n(vAtoms[i].c_ptr(), signalSize, pMatDict + signalSize*i);
	}
	engine.addArray(pMatDict, signalSize, atomCount, false, "D");
	engine.addDoubleScalar(lambda, "L1");
	engine.eval("[idx,coeff,nnz]=spamsLasso(X,D,L1);");

	double *pNnz = engine.getDblVariablePtr("nnz");
	double *pIdx = engine.getDblVariablePtr("idx");
	double *pCoeff = engine.getDblVariablePtr("coeff");
	int supportSize = int(*pNnz);
	
	fa.clear();
	for (int i = 0; i < supportSize; ++i)
		fa.addItem((int)pIdx[i] - 1, pCoeff[i]);

	engine.removeVariable("X");
	engine.removeVariable("D");
	engine.removeVariable("L1");
	delete []pMatDict;
}

void singleChannelSparseApproximate(const VecNd& vSignal, const Dictionary& dict, SparseCoding& sc, SparseApproximationOptions opts)
{
    std::vector<SparseCoding> vsc;
    multiChannelSparseApproximate(std::vector< VecNd>{vSignal}, dict, vsc, opts);
    sc = vsc[0];
}

void multiChannelSparseApproximate(const std::vector<VecNd>& vSignals, const Dictionary& dict, std::vector<SparseCoding>& vCodings, SparseApproximationOptions opts)
{
	runtime_assert(vSignals[0].size() == dict.atomDim(), "signal and atom size not equal");
	runtime_assert(opts.mCodingSize > 0 && opts.mCodingSize <= dict.size(), "Illegal coding size");

	const int channelCount = (int)vSignals.size();
	const int signalDim = dict.atomDim();
	vCodings.resize(channelCount);

	if (opts.mApproxMethod == ZGeom::SA_Truncation)
	{
		for (int i = 0; i < opts.mCodingSize; ++i) {
			for (int c = 0; c < channelCount; ++c) {
				double innerProd = dict[i].dot(vSignals[c]);
				vCodings[c].addItem(i, innerProd);
			}
		}
	}
	else if (opts.mApproxMethod == ZGeom::SA_OMP)
	{
		for (int c = 0; c < channelCount; ++c) {
			ZGeom::OMP(vSignals[c], dict.getAtoms(), opts.mCodingSize, vCodings[c]);
		}
	}
    else if (opts.mApproxMethod == ZGeom::SA_OMP_SPAMS)
    {
        for (int c = 0; c < channelCount; ++c)
            ZGeom::OMP_SPAMS(*opts.mMatlabEngine, vSignals[c], dict.getAtoms(), opts.mCodingSize, vCodings[c]);
    }
	else if (opts.mApproxMethod == ZGeom::SA_SOMP)
	{
		std::vector<ZGeom::SparseCoding*> vCodingPtrs;
		for (int c = 0; c < channelCount; ++c) vCodingPtrs.push_back(&vCodings[c]);
		ZGeom::SimultaneousOMP(vSignals, dict.getAtoms(), opts.mCodingSize, vCodingPtrs);
	}
}

} // end of namespace
