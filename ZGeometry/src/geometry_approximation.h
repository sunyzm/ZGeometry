#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/sparse_approximation.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include "MeshHelper.h"

enum DictionaryType {
	DT_UNIT, DT_Fourier, DT_FourierSpikes, DT_SGW1, DT_SGW2, DT_SGW3, DT_SGW4, 
	DT_SGW5, DT_SGW3MHB, DT_SGW4MHB, DT_SGW5MHB, DT_HK
};

std::vector<ZGeom::SparseCoding> multiExtractFront(const std::vector<ZGeom::SparseCoding>& vCodings, int n);
std::vector<ZGeom::SparseCoding> multiExtractAfter(const std::vector<ZGeom::SparseCoding>& vCodings, int n);

void computeDictionary(DictionaryType dictType, const ZGeom::EigenSystem& es, ZGeom::Dictionary& dict);
void calSGWDict(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::Dictionary& dict);
void calHKDict(const ZGeom::EigenSystem& es, double timescale, ZGeom::Dictionary& dict);
void multiDictSparseDecomposition(const MeshCoordinates& coordInput, const std::vector<const ZGeom::Dictionary*>& vDicts, const std::vector<int>& vNNZ, std::vector<std::vector<ZGeom::SparseCoding> > &vFinalCodings);

ZGeom::VecNd singleChannelSparseInpaint(const ZGeom::VecNd& vSignal, const std::vector<int>& vMask, const ZGeom::Dictionary& dict, ZGeom::SparseCoding& sc);

struct SparseCodingOptions
{
	SparseCodingOptions() : mApproxMethod(ZGeom::SA_SOMP), mCodingAtomCount(-1) {}
	ZGeom::SparseApproxMethod mApproxMethod;
	int mCodingAtomCount;
	double mEpsilon;
	double lambda1, lambda2;
	double mResidual;
};

class SubMeshApprox
{
public:
	friend class ShapeApprox;

    SubMeshApprox(){}
    void init() { mMeshHelper.init(&mSubMesh); }
	int subMeshSize() const { return mSubMesh.vertCount(); }
	const std::vector<int>& mappedIdx() const { return mMappedIdx; }
	void prepareEigenSystem(LaplacianType laplacianType, int mEigenCount);
	void constructDict(DictionaryType dictType);
	void doSparseCoding(ZGeom::SparseApproxMethod approxMethod, int codingAtomCount);
	void sparseReconstruct(int reconstructAtomCount);
	void sparseReconstructStep(int step);
	int dictSize() const { return mDict.atomCount(); }
	int codingSize() const { return mCoding[0].size(); }
	const ZGeom::Dictionary& getDict() const { return mDict; }
	void computeSparseCoding(const std::vector<double>& vSignal, SparseCodingOptions& opts, ZGeom::SparseCoding& vApprox);
	void computeSparseCoding(const ZGeom::VecNd& vSignal, SparseCodingOptions& opts, ZGeom::SparseCoding& vApprox);
	const std::vector<ZGeom::SparseCodingItem>& getSparseCoding(int c) const { return mCoding[c]; }

private:
    SubMeshApprox(const SubMeshApprox&);

private:
	CMesh mSubMesh;
	MeshHelper mMeshHelper;
	std::vector<int> mMappedIdx;
	ZGeom::EigenSystem mEigenSystem;
	ZGeom::Dictionary mDict;
	std::vector<ZGeom::SparseCodingItem> mCoding[3];
	MeshCoordinates mReconstructedCoord;
};

class ShapeApprox
{
public:
	friend class ShapeEditor;	

	ShapeApprox() : mOriginalMesh(nullptr) {}
	void init(CMesh* mesh);
	void doSegmentation(int maxSize);
	void doEigenDecomposition(LaplacianType lapType, int eigenCount);
	void constructDictionaries(DictionaryType dictType);
	void findSparseRepresentationBySize(ZGeom::SparseApproxMethod codingMethod, int codingSize);
	void findSparseRepresentationByRatio(ZGeom::SparseApproxMethod codingMethod, double basisRatio, bool exploitSparsity);
	void findSparseRepresentationByBasisRatio(ZGeom::SparseApproxMethod codingMethod, double basisRatio);
	void findSparseRepresentationByCompressionRatio(ZGeom::SparseApproxMethod codingMethod, double compressionRatio);
	void doSparseReconstructionBySize(int reconstructSize, MeshCoordinates& approxCoord);
	void doSparseReconstructionByRatio(double basisRatio, MeshCoordinates& approxCoord, bool exploitSparsity);
	void doSparseReconstructionByBasisRatio(double basisRatio, MeshCoordinates& approxCoord);
	void doSparseReconstructionByCompressionRatio(double compressionRatio, MeshCoordinates& approxCoord);
	void doSparseReconstructionStepping(int totalSteps, std::vector<MeshCoordinates>& contCoords);
	void integrateSubmeshApproximation(MeshCoordinates& integratedApproxCoord);
	int  partitionCount() const { return mSubMeshApprox.size(); }
	ZGeom::Dictionary getSubDictionary(int idx) { return mSubMeshApprox[idx]->mDict; }
	const ZGeom::EigenSystem& getSubEigenSystem(int idx) { return mSubMeshApprox[idx]->mEigenSystem; }

private:
    void clear();
    ShapeApprox(const ShapeApprox&);

private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox*> mSubMeshApprox;
	std::vector<int> mPartIdx;
};