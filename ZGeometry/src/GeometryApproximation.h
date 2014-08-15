#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/Approximation.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"

using ZGeom::SparseApproxMethod;
enum DictionaryType {
	DT_UNIT, DT_Fourier, DT_FourierSpikes, DT_SGW1, DT_SGW2, DT_SGW3, DT_SGW4, 
	DT_SGW5, DT_SGW3MHB, DT_SGW4MHB, DT_SGW5MHB, DT_HK
};

void computeDictionary(DictionaryType dictType, const ZGeom::EigenSystem& es, ZGeom::Dictionary& dict);
void calSGWDict(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::Dictionary& dict);
void calHKDict(const ZGeom::EigenSystem& es, double timescale, ZGeom::Dictionary& dict);


struct SparseCodingOptions
{
	SparseCodingOptions() : mApproxMethod(ZGeom::SA_SOMP), mCodingAtomCount(-1) {}

	SparseApproxMethod mApproxMethod;
	int mCodingAtomCount;
	double mEpsilon;
	double lambda1, lambda2;
	double mResidual;
};

class SubMeshApprox
{
public:
	friend class ShapeApprox;

	void init() { mMeshProcessor.init(&mSubMesh); }
	int subMeshSize() const { return mSubMesh.vertCount(); }
	const std::vector<int>& mappedIdx() const { return mMappedIdx; }
	void prepareEigenSystem(LaplacianType laplacianType, int mEigenCount);
	void constructDict(DictionaryType dictType);
	void doSparseCoding(SparseApproxMethod approxMethod, int codingAtomCount);
	void sparseReconstruct(int reconstructAtomCount);
	void sparseReconstructStep(int step);
	int dictSize() const { return mDict.atomCount(); }
	int codingSize() const { return mCoding[0].size(); }
	const ZGeom::Dictionary& getDict() const { return mDict; }
	void computeSparseCoding(const std::vector<double>& vSignal, SparseCodingOptions& opts, ZGeom::SparseCoding& vApprox);
	void computeSparseCoding(const ZGeom::VecNd& vSignal, SparseCodingOptions& opts, ZGeom::SparseCoding& vApprox);
	const std::vector<ZGeom::SparseCodingItem>& getSparseCoding(int c) const { return mCoding[c]; }

private:
	CMesh mSubMesh;
	DifferentialMeshProcessor mMeshProcessor;
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
	ShapeApprox() : mOriginalMesh(NULL) {}
	void init(CMesh* mesh);
	void doSegmentation(int maxSize);
	void doEigenDecomposition(LaplacianType lapType, int eigenCount);
	void constructDictionaries(DictionaryType dictType);
	void findSparseRepresentationBySize(SparseApproxMethod codingMethod, int codingSize);
	void findSparseRepresentationByRatio(SparseApproxMethod codingMethod, double basisRatio, bool exploitSparsity);
	void findSparseRepresentationByBasisRatio(SparseApproxMethod codingMethod, double basisRatio);
	void findSparseRepresentationByCompressionRatio(SparseApproxMethod codingMethod, double compressionRatio);
	void doSparseReconstructionBySize(int reconstructSize, MeshCoordinates& approxCoord);
	void doSparseReconstructionByRatio(double basisRatio, MeshCoordinates& approxCoord, bool exploitSparsity);
	void doSparseReconstructionByBasisRatio(double basisRatio, MeshCoordinates& approxCoord);
	void doSparseReconstructionByCompressionRatio(double compressionRatio, MeshCoordinates& approxCoord);
	void doSparseReconstructionStepping(int totalSteps, std::vector<MeshCoordinates>& contCoords);
	void integrateSubmeshApproximation(MeshCoordinates& integratedApproxCoord);
	int  partitionCount() const { return mSubMeshApprox.size(); }
	ZGeom::Dictionary getSubDictionary(int idx) { return mSubMeshApprox[idx].mDict; }
	const ZGeom::EigenSystem& getSubEigenSystem(int idx) { return mSubMeshApprox[idx].mEigenSystem; }

private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox> mSubMeshApprox;
	std::vector<int> mPartIdx;
};