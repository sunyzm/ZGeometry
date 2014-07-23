#ifndef ZGEOM_APPROXIMATION_H
#define ZGEOM_APPROXIMATION_H
#include <vector>
#include <tuple>
#include <functional>
#include "util.h"
#include "VecN.h"
#include "DenseMatrix.h"
#include "MatVecArithmetic.h"
#include "MatlabEngineWrapper.h"
#include "SparseRepresentation.h"

namespace ZGeom {

void GeneralizedMultiChannelFourierApprox(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, InnerProdcutFunc innerProdFunc);

void MatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit);	
void GeneralizedMP( const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit, const InnerProdcutFunc& innerProdFunc);

void OrthogonalMatchingPursuit(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit);
void OMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit);
void GeneralizedOMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit, InnerProdcutFunc innerProdFunc);

void SimultaneousMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuit, double p = 2.);
void GeneralizedSimultaneousMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuit, const InnerProdcutFunc& innerProdFunc, double p = 2);

//// compute Multi-channel OMP with MKL and regular inner product
void SimultaneousOMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, double p = 2.);
void GeneralizedSimultaneousOMP(const std::vector<VecNd>& vSignals, const std::vector<VecNd>& vBasis, int nSelected, std::vector<SparseCoding*>& vPursuits, const InnerProdcutFunc& innerProdFunc, double p = 2.);

//// compute OMP with MKL and regular inner product
void OrthogonalMatchingPursuit_AMP(const VecNd& vSignal, const std::vector<VecNd>& vBasis, int nSelected, SparseCoding& vPursuit);

//// compute OMP with Matlab and customized InnerProductFunc
void GeneralizedOMP_MATLAB(const VecNd& vSignal, const std::vector<VecNd>& vBasis, InnerProdcutFunc innerProdFunc, int nSelected, SparseCoding& vPursuit, MatlabEngineWrapper& engine);

void OMP_SPAMS(MatlabEngineWrapper& engine, const VecNd& vSignal, const std::vector<VecNd>& vAtoms, int supportSize, SparseCoding& fa);
void LASSO_SPAMS(MatlabEngineWrapper& engine, const VecNd& vSignal, const std::vector<VecNd>& vAtoms, double lambda, SparseCoding& fa);

enum SparseApproxMethod { SA_Truncation, SA_SMP, SA_OMP, SA_SOMP, SA_LASSO };

struct SparseApproximationOptions
{
	int mCodingSize;
	bool mSimultaneousCoding;
	SparseApproxMethod mApproxMethod;
};

void multiChannelSparseApproximate(const std::vector<VecNd>& vSignals, const Dictionary& dict, std::vector<SparseCoding>& vCodings, SparseApproximationOptions opts);
void singleChannelSparseRecontruct(const Dictionary& dict, const SparseCoding& coding, VecNd& reconstructed);
void multiChannelSparseReconstruct(const Dictionary& dict, const std::vector<SparseCoding>& vCodings, std::vector<VecNd>& vReconstructed);

} // end of namespace

#endif
