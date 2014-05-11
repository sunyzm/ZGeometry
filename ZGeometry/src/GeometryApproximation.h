#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"

enum DictionaryType {DT_Fourier, DT_FourierSpikes, DT_SGW3, DT_SGW4, DT_SGW5, DT_SGW3MHB, DT_SGW4MHB, DT_SGW5MHB};
enum SparseApproxMethod {SA_Truncation, SA_SMP, SA_SOMP, SA_LASSO};

struct SparseCodingOptions
{
	SparseCodingOptions() : mApproxMethod(SA_SOMP), mCodingAtomCount(-1) {}

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
	void computeSparseCoding(const std::vector<double>& vSignal, SparseCodingOptions& opts, ZGeom::FunctionApproximation& vApprox);
	void computeSparseCoding(const ZGeom::VecNd& vSignal, SparseCodingOptions& opts, ZGeom::FunctionApproximation& vApprox);

private:
	CMesh mSubMesh;
	DifferentialMeshProcessor mMeshProcessor;
	std::vector<int> mMappedIdx;
	ZGeom::EigenSystem mEigenSystem;
	ZGeom::Dictionary mDict;
	std::vector<ZGeom::ApproxItem> mCoding[3];
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
	
private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox> mSubMeshApprox;
	std::vector<int> mPartIdx;
};