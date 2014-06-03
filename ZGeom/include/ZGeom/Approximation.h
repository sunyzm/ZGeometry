#ifndef ZGEOM_APPROXIMATION_H
#define ZGEOM_APPROXIMATION_H
#include <vector>
#include <tuple>
#include <functional>
#include "util.h"
#include "VecN.h"
#include "DenseMatrix.h"
#include "MatlabEngineWrapper.h"
#include "Dictionary.h"

namespace ZGeom
{

struct ApproxItem
{
	ApproxItem() : mRes(0), mBasisIdx(-1), mCoeff(0) {}
	ApproxItem(int i, double c) : mRes(0), mBasisIdx(i), mCoeff(c) {}
	ApproxItem(double r, int i, double c) : mRes(r), mBasisIdx(i), mCoeff(c) {}

	double res() const { return mRes; }
	double& res() { return mRes; }
	int index() const { return mBasisIdx; }
	int& index() { return mBasisIdx; }
	double coeff() const { return mCoeff; }
	double& coeff() { return mCoeff; }

	double mRes;
	int mBasisIdx;
	double mCoeff;
};


class FunctionApproximation
{
public:
	FunctionApproximation() {}

	void addItem(double r, int i, double c)
	{
		mApproxItems.push_back(ApproxItem(r, i, c));
	}

	void clear() { mApproxItems.clear(); }
	const std::vector<ApproxItem>& getApproxItems() const { return mApproxItems; }
	ApproxItem& operator [] (int i) { return mApproxItems[i]; }
	const ApproxItem& operator [] (int i) const { return mApproxItems[i]; }
	ApproxItem& back() { return mApproxItems.back(); }
	void resize(int n) { mApproxItems.resize(n); }
	int size() const { return (int)mApproxItems.size(); }

	std::vector<int> getAllAtomIndex() const 
	{
		std::vector<int> vIdx;
		for (auto ai : mApproxItems) vIdx.push_back(ai.index());
		return vIdx;
	}

	friend std::ostream& operator << (std::ostream &os, const FunctionApproximation& fa) 
	{
		for (auto t : fa.mApproxItems) {
			os << t.index() << ", " << t.coeff() << '\n';
		}
		return os;
	}

private:
	std::vector<ApproxItem> mApproxItems;
};

double RegularProductFunc(const VecN<double>& v1, const VecN<double>& v2);

VecNd ReconstructApproximationSingleChannel(const Dictionary& dict, const FunctionApproximation& approx);
VecNd ReconstructApproximationSingleChannel(const std::vector<VecNd>& vAtoms, const FunctionApproximation& approx);
void ReconstructApproximationMultiChannel(const std::vector<VecNd>& vAtoms, const std::vector<FunctionApproximation>& vApprox, std::vector<VecNd>& vReconstructed);

void GeneralizedMultiChannelFourierApprox(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<FunctionApproximation*>& vPursuits, InnerProdcutFunc innerProdFunc);

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

void OMP_SPAMS(MatlabEngineWrapper& engine, const VecNd& vSignal, const std::vector<VecNd>& vAtoms, int supportSize, FunctionApproximation& fa);
void LASSO_SPAMS(MatlabEngineWrapper& engine, const VecNd& vSignal, const std::vector<VecNd>& vAtoms, double lambda, FunctionApproximation& fa);

} // end of namespace

#endif
