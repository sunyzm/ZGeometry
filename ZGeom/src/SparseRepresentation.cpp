#include "SparseRepresentation.h"
#include "MatVecArithmetic.h"

namespace ZGeom {

ZGeom::VecNd singleChannelDenseReconstruct(const Dictionary& dict, const VecNd& vCoeff)
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
		vResults[c] = singleChannelDenseReconstruct(dict, vCoeffs[c]);
	}
	return vResults;
}

ZGeom::VecNd singleChannelSparseReconstruct(const Dictionary& dict, const SparseCoding& coding)
{
	return singleChannelSparseReconstruct(dict.getAtoms(), coding);
}

ZGeom::VecNd singleChannelSparseReconstruct(const std::vector<ZGeom::VecNd>& vDictAtoms, const SparseCoding& coding)
{
	ZGeom::VecNd vApprox(vDictAtoms[0].size());
	for (auto item : coding.getApproxItems()) {
		vApprox += vDictAtoms[item.index()] * item.coeff();
	}
	return vApprox;
}

void singleChannelSparseReconstruct(const Dictionary& dict, const SparseCoding& coding, VecNd& signalReconstructed)
{
	signalReconstructed = singleChannelSparseReconstruct(dict, coding);
}

void multiChannelSparseReconstruct(const Dictionary& dict, const std::vector<SparseCoding>& vCodings, std::vector<VecNd>& vReconstructed)
{
	int channelCount = (int)vCodings.size();
	vReconstructed.resize(channelCount);

	for (int c = 0; c < channelCount; ++c) {
		singleChannelSparseReconstruct(dict, vCodings[c], vReconstructed[c]);
	}
}

void multiChannelSparseReconstruct(const std::vector<VecNd>& vAtoms, const std::vector<SparseCoding>& vCodings, std::vector<VecNd>& vReconstructed)
{
	int nChannels = (int)vCodings.size();
	vReconstructed.resize(nChannels);

	for (int i = 0; i < nChannels; ++i) {
		vReconstructed[i] = singleChannelSparseReconstruct(vAtoms, vCodings[i]);
	}
}

ZGeom::VecNd dictVecInnerProducts(const Dictionary& dict, const VecNd& vec)
{
	int vecSize = vec.size();
	int dictSize = dict.size();
	VecNd vCoeff(dictSize);
	for (int i = 0; i < dictSize; ++i)
		vCoeff[i] = RegularProductFunc(dict[i], vec);
	return vCoeff;
}

void splitSparseCoding(const SparseCoding& coding1, SparseCoding& coding2, SparseCoding& coding3, std::function<bool(const SparseCodingItem&)> funcSplit)
{
	coding2.clear();
	coding3.clear();

	for (const SparseCodingItem& codingItem : coding1.getApproxItems()) {
		if (funcSplit(codingItem)) coding2.addItem(codingItem);
		else coding3.addItem(codingItem);
	}
}

void filterSparseCoding(const SparseCoding& coding1, SparseCoding& coding2, std::function<bool(const SparseCodingItem&)> funcFilter)
{
	coding2.clear();
	for (auto codingItem : coding1.getApproxItems()) {
		if (funcFilter(codingItem)) coding2.addItem(codingItem);
	}
}

void multiChannelSplitSparseCoding(const std::vector<SparseCoding>& vCodings1, std::vector<SparseCoding>& vCodings2, std::vector<SparseCoding>& vCodings3, std::function<bool(const SparseCodingItem&)> funcSplit)
{
	int nChannels = (int)vCodings1.size();
	vCodings2.resize(nChannels);
	vCodings3.resize(nChannels);

	for (int c = 0; c < nChannels; ++c)
		splitSparseCoding(vCodings1[c], vCodings2[c], vCodings3[c], funcSplit);
}

} // end of namespace