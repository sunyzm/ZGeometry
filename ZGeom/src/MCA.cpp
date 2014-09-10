#include "MCA.h"
#include "Approximation.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>

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
    using std::vector;

	/// initialization
	int dictCount = (int)vDicts.size();
	int signalSize = vSignal.size();
	vCodings.resize(dictCount);	
	vector<VecNd> vDenseCoeff(dictCount);
	vector<VecNd> vComponents(dictCount);
	for (int k = 0; k < dictCount; ++k) {
		vDenseCoeff[k].resize(vDicts[k]->size(), 0);
		vComponents[k].resize(signalSize, 0);
	}
	VecNd vRes(vSignal);
	double lambdaThresh;

	// initial threshold lambda
	VecNd vMaxCoeff(dictCount);
	for (int k = 0; k < dictCount; ++k) {
		int dictSizeK = vDicts[k]->size();
		VecNd vCoeffK = dictVecInnerProducts(*vDicts[k], vRes);
        vMaxCoeff[k] = vCoeffK.inftyNorm();
	}

	switch (opts->threshStrategy)
	{
	case MCAoptions::MEAN_OF_MAX:
        lambdaThresh = vMaxCoeff.mean();
		break;
	case MCAoptions::MIN_OF_MAX:
        lambdaThresh = vMaxCoeff.min_element();
		break;
    case MCAoptions::MAX_OF_MAX:
        lambdaThresh = vMaxCoeff.max_element();
        break;
	case MCAoptions::SECOND_TO_MAX:
        std::sort(vMaxCoeff.c_ptr(), vMaxCoeff.c_ptr_end(), std::greater<double>());
        lambdaThresh = vMaxCoeff[1];
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

double calThresh(const std::vector<const Dictionary*>& vDicts, const VecNd& vRes, MCAoptions::ThresholdingStrategy ts)
{
    int dictCount = (int)vDicts.size();
    VecNd vMaxCoeff(dictCount);
    for (int k = 0; k < dictCount; ++k) {
        int dictSizeK = vDicts[k]->size();
        VecNd vCoeffK = dictVecInnerProducts(*vDicts[k], vRes);
        vMaxCoeff[k] = vCoeffK.inftyNorm();
    }

    double lambdaThresh;
    switch (ts)
    {
    case MCAoptions::MEAN_OF_MAX:
        lambdaThresh = vMaxCoeff.mean();
        break;
    case MCAoptions::MIN_OF_MAX:
        lambdaThresh = vMaxCoeff.min_element();
        break;
    case MCAoptions::MAX_OF_MAX:
        lambdaThresh = vMaxCoeff.max_element();
        break;
    case MCAoptions::SECOND_TO_MAX:
        std::sort(vMaxCoeff.c_ptr(), vMaxCoeff.c_ptr_end(), std::greater<double>());
        lambdaThresh = vMaxCoeff[1];
        break;
    }
    return lambdaThresh;
}

void singleChannelMCA2(const VecNd& vSignal, const std::vector<const Dictionary*>& vDicts, std::vector<SparseCoding>& vCodings, MCAoptions* opts)
{
    using std::vector;

    /// initialization
    int dictCount = (int)vDicts.size();
    int signalSize = vSignal.size();
    vCodings.resize(dictCount);
    vector<VecNd> vDenseCoeff(dictCount);
    vector<VecNd> vComponents(dictCount);
    for (int k = 0; k < dictCount; ++k) {
        vDenseCoeff[k].resize(vDicts[k]->size(), 0);
        vComponents[k].resize(signalSize, 0);
    }
    VecNd vRes(vSignal);
    double lambdaThresh = calThresh(vDicts, vRes, opts->threshStrategy);
    double minResNorm(0);    // minimum residual norm

    /// main iteration
    for (int t = 0; t < opts->nIter; ++t)
    {
        std::vector<VecNd> vTmpCoeff(dictCount);
        for (int k = 0; k < dictCount; ++k) {
            // compute marginal residual
            VecNd vResK = vRes + vComponents[k];
            // compute k-th component coefficients
            vTmpCoeff[k] = dictVecInnerProducts(*vDicts[k], vResK);
            // update k-th component coefficients by thresholding 
            MCAoptions::ThresholdingMode threshMode = opts->threshMode;
            std::transform(vTmpCoeff[k].c_ptr(), vTmpCoeff[k].c_ptr_end(), vTmpCoeff[k].c_ptr(), FuncTresholding(threshMode, lambdaThresh));
            // update k-th component
            vComponents[k] = singleChannelDenseReconstruct(*vDicts[k], vTmpCoeff[k]);
        }

        // update residual
        vRes = vSignal;
        for (int k = 0; k < dictCount; ++k) vRes -= vComponents[k];        
        double resNorm = vRes.norm2();
        if (t == 0 || resNorm < minResNorm) {
            minResNorm = resNorm;
            vDenseCoeff = vTmpCoeff;
        }
        else {
            std::cout << "Total MCA iterations#: " << t << "\t";
            break;
        }

        // update threshold lambda
        lambdaThresh = calThresh(vDicts, vRes, opts->threshStrategy);
    }

    for (int k = 0; k < dictCount; ++k) {
        vCodings[k].fromDense(vDenseCoeff[k], 0.01);
        std::cout << "MCA" << k+1 << "#: " << vCodings[k].size() << '\t';
    }
    std::cout << "\n";
}

void singleChannelMCA3(const VecNd& vSignal, const std::vector<const Dictionary*>& vDicts, std::vector<SparseCoding>& vCodings, MCAoptions* opts)
{
    using std::vector;
    /// initialization
    int dictCount = (int)vDicts.size();
    int signalSize = vSignal.size();
    vCodings.resize(dictCount);
    vector<VecNd> vDenseCoeff(dictCount);
    vector<VecNd> vComponents(dictCount);
    for (int k = 0; k < dictCount; ++k) {
        vDenseCoeff[k].resize(vDicts[k]->size(), 0);
        vComponents[k].resize(signalSize, 0);
    }
    VecNd vRes(vSignal);
    double lambdaThresh = calThresh(vDicts, vRes, opts->threshStrategy);
    double minResNorm(0);    // minimum residual norm
    SparseApproximationOptions approxOpts;
    approxOpts.mCodingSize = 50;
    approxOpts.mApproxMethod = ZGeom::SA_OMP;
    MCAoptions::ThresholdingMode threshMode = opts->threshMode;

    /// main iteration
    for (int t = 0; t < opts->nIter; ++t)
    {
        std::vector<VecNd> vTmpCoeff(dictCount);
        for (int k = 0; k < dictCount; ++k) {
            // compute marginal residual
            //VecNd vResK = vRes + vComponents[k];
            VecNd vResK(vSignal);
            for (int j = 0; j < dictCount; ++j)
                if (j != k) vResK -= vComponents[j];
            // compute k-th component coefficients
            SparseCoding sc;
            ZGeom::OMP(vResK, vDicts[k]->getAtoms(), approxOpts.mCodingSize, sc);
            vTmpCoeff[k] = sc.toDense(vDicts[k]->size());
            //vTmpCoeff[k] = dictVecInnerProducts(*vDicts[k], vResK);
            // update k-th component coefficients by thresholding 
            //std::transform(vTmpCoeff[k].c_ptr(), vTmpCoeff[k].c_ptr_end(), vTmpCoeff[k].c_ptr(), FuncTresholding(threshMode, lambdaThresh));
            // update k-th component
            vComponents[k] = singleChannelDenseReconstruct(*vDicts[k], vTmpCoeff[k]);
        }

        // update residual
        vRes = vSignal;
        for (int k = 0; k < dictCount; ++k) vRes -= vComponents[k];
        double resNorm = vRes.norm2();
        if (t == 0 || fabs(resNorm - minResNorm) > 0.01) 
        {
            minResNorm = resNorm;
            vDenseCoeff = vTmpCoeff;
            std::cout << "(Iteration " << t << ") residual: " << resNorm << "\tlambda: " << lambdaThresh << '\n';
        }
        else {
            std::cout << "Total MCA iterations#: " << t << "\t";
            break;
        }

        // update threshold lambda
        lambdaThresh = calThresh(vDicts, vRes, opts->threshStrategy);
    }

    for (int k = 0; k < dictCount; ++k) {
        vCodings[k].fromDense(vDenseCoeff[k]);
        std::cout << "MCA" << k + 1 << "#: " << vCodings[k].size() << '\t';
    }
    std::cout << "\n";
}

} // end of namespace