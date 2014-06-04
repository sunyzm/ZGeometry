#include "SparseRepresentation.h"
#include "MatVecArithmetic.h"

namespace ZGeom {

ZGeom::VecNd DenseReconstructSingleChannel(const Dictionary& dict, const VecNd& vCoeff)
{
	assert(dict.atomCount() == vCoeff.size());
	int atomCount = dict.atomCount();
	int signalSize = dict.atomDim();
	VecNd vResult(signalSize, 0);
	for (int i = 0; i < atomCount; ++i) {
		vResult += dict[i] * vCoeff[i];
	}
	return vResult;
}

std::vector<VecNd> DenseReconstructMultiChannel(const Dictionary& dict, const std::vector<VecNd>& vCoeffs)
{
	int channelCount = (int)vCoeffs.size();
	std::vector<VecNd> vResults(channelCount);
	for (int c = 0; c < channelCount; ++c) {
		vResults[c] = DenseReconstructSingleChannel(dict, vCoeffs[c]);
	}
	return vResults;
}

ZGeom::VecNd SparseReconstructSingleChannel(const Dictionary& dict, const SparseCoding& approx)
{
	ZGeom::VecNd vApprox(dict.atomDim(), 0);
	for (auto item : approx.getApproxItems())
	{
		vApprox += dict[item.index()] * item.coeff();
	}

	return vApprox;
}

ZGeom::VecNd SparseReconstructSingleChannel(const std::vector<ZGeom::VecNd>& vAtoms, const SparseCoding& approx)
{
	ZGeom::VecNd vApprox(vAtoms[0].size());
	for (auto item : approx.getApproxItems()) {
		vApprox += vApprox[item.index()] * item.coeff();
	}
	return vApprox;
}

void ReconstructApproximationMultiChannel(const std::vector<VecNd>& vAtoms, const std::vector<SparseCoding>& vApprox, std::vector<VecNd>& vReconstructed)
{
	int nChannels = (int)vApprox.size();
	vReconstructed.resize(nChannels);
	for (int i = 0; i < nChannels; ++i)
		vReconstructed[i] = SparseReconstructSingleChannel(vAtoms, vApprox[i]);
}

ZGeom::VecNd DictVecInnerProducts(const Dictionary& dict, const VecNd& vec)
{
	int vecSize = vec.size();
	int dictSize = dict.size();
	VecNd vCoeff(dictSize);
	for (int i = 0; i < dictSize; ++i)
		vCoeff[i] = RegularProductFunc(dict[i], vec);
	return vCoeff;
}

} // end of namespace