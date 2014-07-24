#include "MCA.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ZGeom {

double FuncAbs(double x) { return std::fabs(x); }

class FuncTresholding
{
public:
	FuncTresholding(MCAoptions::ThresholdingMode m, double t) : mMode(m), mThresh(std::fabs(t)) {}
	double operator()(double x) {
		if (mMode == MCAoptions::HARD_THRESH) {
			return std::fabs(x) > mThresh ? x : 0;
		} else if (mMode == MCAoptions::SOFT_THRESH) {
			if (x > mThresh) return x - mThresh;
			else if (x < -mThresh) return x + mThresh;
			else return 0;
		}
		else throw std::logic_error("Undefined ThroshodlingMode");
	}

private:
	MCAoptions::ThresholdingMode mMode;
	double mThresh;
};

void singleChannelMCA(const VecNd& vSignal, const std::vector<const Dictionary*>& vDicts, std::vector<SparseCoding>& vCodings, MCAoptions* opts)
{
	/// initialization
	int dictCount = (int)vDicts.size();
	int signalSize = vSignal.size();
	vCodings.resize(dictCount);	
	std::vector<VecNd> vDenseCoeff(dictCount);
	std::vector<VecNd> vComponents(dictCount);
	for (int k = 0; k < dictCount; ++k) {
		vDenseCoeff[k].resize(vDicts[k]->size(), 0);
		vComponents[k].resize(signalSize, 0);
	}
	VecNd vRes(vSignal);
	double lambdaThresh;

	// initial threshold lambda
	std::vector<double> vMaxCoeff(dictCount);
	for (int k = 0; k < dictCount; ++k) {
		int dictSizeK = vDicts[k]->size();
		VecNd vCoeffK = dictVecInnerProducts(*vDicts[k], vRes);
		std::transform(vCoeffK.c_ptr(), vCoeffK.c_ptr_end(), vCoeffK.c_ptr(), FuncAbs);
		vMaxCoeff[k] = *std::max_element(vCoeffK.c_ptr(), vCoeffK.c_ptr_end());
	}

	switch (opts->threshStrategy)
	{
	case MCAoptions::MEAN_OF_MAX:
		lambdaThresh = std::accumulate(vMaxCoeff.begin(), vMaxCoeff.end(), 0.) / vMaxCoeff.size();
		break;
	case MCAoptions::MIN_OF_MAX:
		lambdaThresh = *std::min_element(vMaxCoeff.begin(), vMaxCoeff.end());
		break;
	case MCAoptions::SECOND_TO_MAX:
		vMaxCoeff.erase(std::max_element(vMaxCoeff.begin(), vMaxCoeff.end()));
		lambdaThresh = *std::max_element(vMaxCoeff.begin(), vMaxCoeff.end());
		break;
	case MCAoptions::MAX_OF_MAX:
		lambdaThresh = *std::max_element(vMaxCoeff.begin(), vMaxCoeff.end());
		break;
	}

	lambdaThresh /= 2;
	double deltaThresh = (lambdaThresh - opts->stopCriterion) / (opts->nIter - 1);
	double minResMag;

	/// main iteration
	for (int t = 0; t < opts->nIter; ++t) 
	{
		std::vector<VecNd> vTmpCoeff(dictCount);
		for (int k = 0; k < dictCount; ++k) {
			// compute marginal residual
			VecNd vResK = vRes + vComponents[k];
			// update k-th component coefficients by thresholding 
			vTmpCoeff[k] = dictVecInnerProducts(*vDicts[k], vResK);
			MCAoptions::ThresholdingMode threshMode = MCAoptions::SOFT_THRESH;
			//if (t > 1) threshMode = MCAoptions::HARD_THRESH;
			
			std::transform(vTmpCoeff[k].c_ptr(), vTmpCoeff[k].c_ptr_end(), vTmpCoeff[k].c_ptr(), FuncTresholding(threshMode, lambdaThresh));
			// update k-th component
			vComponents[k] = singleChannelDenseReconstruct(*vDicts[k], vTmpCoeff[k]);
		}

		// update residual
		vRes = vSignal;
		for (int k = 0; k < dictCount; ++k) {
			vRes -= vComponents[k];
		}
		double resMag = vRes.norm2();
		if (t == 0 || resMag < minResMag) {
			minResMag = resMag;
			vDenseCoeff = vTmpCoeff;
		}
		else break;

		// update threshold lambda
		lambdaThresh -= deltaThresh;
	}

	for (int k = 0; k < dictCount; ++k) {
		vCodings[k].fromDense(vDenseCoeff[k], 0.1);
	}
}


} // end of namespace