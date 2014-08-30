#include <algorithm>
#include <ppl.h>
#include <concurrent_vector.h>
#include <mkl.h>
#include "SparseRepresentation.h"
#include "MatVecArithmetic.h"
#include "DenseMatrix.h"



namespace ZGeom {

void Dictionary::mergeDict(const Dictionary& d2)
{
	for (const VecNd& atom : d2.mAtoms) mAtoms.push_back(atom);
	mDim = (int)mAtoms.size();
}

double Dictionary::calCoherence() const
{
    using namespace concurrency;
    int dictSize = (int)mAtoms.size();
    int signalSize = mDim;
    double **atomPtrs = new double*[dictSize];
    for (int i = 0; i < dictSize; ++i) atomPtrs[i] = mAtoms[i].c_ptr();
    concurrent_vector<double> vMaxCohere(dictSize - 1);
    parallel_for(0, dictSize - 1, [&](int i){
        double maxCohere = 0;
        for (int j = i + 1; j < dictSize; ++j) {
            double co = fabs(cblas_ddot(mDim, atomPtrs[i], 1, atomPtrs[j], 1));
            maxCohere = std::max(maxCohere, co);
        }
        vMaxCohere[i] = maxCohere;
    });

    delete[]atomPtrs;
    return *std::max_element(vMaxCohere.begin(), vMaxCohere.end());
}

double Dictionary::calCoherence2() const
{
    using namespace concurrency;
    int dictSize = (int)mAtoms.size();
    int signalSize = mDim;
    double **atomPtrs = new double*[dictSize];
    for (int i = 0; i < dictSize; ++i) atomPtrs[i] = mAtoms[i].c_ptr();

    concurrent_vector<double> vMaxCohere(dictSize - 1);
    parallel_for(0, dictSize - 1, [&](int i){
        double maxCohere = 0;
        for (int j = i + 1; j < dictSize; ++j) {
            double co = fabs(cblas_ddot(mDim, atomPtrs[i], 1, atomPtrs[j], 1));
            maxCohere = std::max(maxCohere, co);
        }
        vMaxCohere[i] = maxCohere;
    });

    delete[]atomPtrs;
    return *std::max_element(vMaxCohere.begin(), vMaxCohere.end());
}

double Dictionary::mutualCoherence(const Dictionary& dict2) const
{
    assert(mDim == dict2.mDim);
    using namespace concurrency;
    int dictSize1 = (int)mAtoms.size(), dictSize2 = dict2.atomCount();
    int signalSize = mDim;

    double **atomPtrs1 = new double*[dictSize1], **atomPtrs2 = new double*[dictSize2];
    for (int i = 0; i < dictSize1; ++i) atomPtrs1[i] = mAtoms[i].c_ptr();
    for (int j = 0; j < dictSize2; ++j) atomPtrs2[j] = dict2.mAtoms[j].c_ptr();

    concurrent_vector<double> vMaxCohere(dictSize1);
    parallel_for(0, dictSize1, [&](int i){
        double maxCohere = 0;
        for (int j = 0; j < dictSize2; ++j) {
            double co = fabs(cblas_ddot(mDim, atomPtrs1[i], 1, atomPtrs2[j], 1));
            maxCohere = std::max(maxCohere, co);
        }
        vMaxCohere[i] = maxCohere;
    });
    
    delete[]atomPtrs1; delete[]atomPtrs2;
    return *std::max_element(vMaxCohere.begin(), vMaxCohere.end());
}

void combineDictionary(const Dictionary& d1, const Dictionary& d2, Dictionary& d3)
{
	assert(d1.atomDim() == d2.atomDim());
	d3.mDim = d1.atomDim();
	d3.mAtoms.clear();
	d3.mAtoms.reserve(d1.atomCount() + d2.atomCount());
	for (const VecNd& atom : d1.mAtoms) d3.mAtoms.push_back(atom);
	for (const VecNd& atom : d2.mAtoms) d3.mAtoms.push_back(atom);
}

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

std::vector<VecNd> multiChannelDenseReconstruct(const Dictionary& dict, const std::vector<VecNd>& vCoeffs)
{
	int channelCount = (int)vCoeffs.size();
	std::vector<VecNd> vResults(channelCount);
	for (int c = 0; c < channelCount; ++c) {
		vResults[c] = singleChannelDenseReconstruct(dict, vCoeffs[c]);
	}
	return vResults;
}

ZGeom::VecNd singleChannelSparseReconstruct(const Dictionary& dict, const SparseCoding& coding, int nCode)
{
	return singleChannelSparseReconstruct(dict.getAtoms(), coding, nCode);
}

ZGeom::VecNd singleChannelSparseReconstruct(const std::vector<ZGeom::VecNd>& vDictAtoms, const SparseCoding& coding, int nCode)
{
	if (nCode < 0 || nCode > coding.size()) nCode = coding.size();
	ZGeom::VecNd vApprox(vDictAtoms[0].size());
	for (int i = 0; i < nCode; ++i) {
		auto item = coding[i];
		vApprox += vDictAtoms[item.index()] * item.coeff();
	}
	return vApprox;
}

void multiChannelSparseReconstruct(const Dictionary& dict, const std::vector<SparseCoding>& vCodings, std::vector<VecNd>& vReconstructed, int nCode)
{
	int channelCount = (int)vCodings.size();
	vReconstructed.resize(channelCount);

	for (int c = 0; c < channelCount; ++c) {
		vReconstructed[c] = singleChannelSparseReconstruct(dict, vCodings[c], nCode);
	}
}

void multiChannelSparseReconstruct(const std::vector<VecNd>& vAtoms, const std::vector<SparseCoding>& vCodings, std::vector<VecNd>& vReconstructed, int nCode)
{
	int nChannels = (int)vCodings.size();
	vReconstructed.resize(nChannels);

	for (int c = 0; c < nChannels; ++c) {
		vReconstructed[c] = singleChannelSparseReconstruct(vAtoms, vCodings[c], nCode);
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

}	// end of namespace

