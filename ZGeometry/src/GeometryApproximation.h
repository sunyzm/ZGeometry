#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"

enum DictionaryType {DT_Fourier, DT_FourierSpikes, DT_SGW3, DT_SGW4, DT_SGW5, DT_SGW3MHB, DT_SGW4MHB, DT_SGW5MHB};
enum SparseApproxMethod {SA_Truncation, SA_SMP, SA_SOMP};

class SubMeshApprox
{
public:
	friend class ShapeApprox;
	struct SparseCoeff 
	{
		SparseCoeff() : mIdx(-1), mCoeff(0) {}
		SparseCoeff(int i, double c) : mIdx(i), mCoeff(c) {}
		int mIdx;
		double mCoeff;
	};
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

private:
	CMesh mSubMesh;
	DifferentialMeshProcessor mMeshProcessor;
	std::vector<int> mMappedIdx;
	ZGeom::EigenSystem mEigenSystem;
	ZGeom::Dictionary mDict;
	std::vector<SparseCoeff> mCoding[3];
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
	//const MeshCoordinates& getApproxCoord() const { return mApproxCoord; }
	int partitionCount() const { return mSubMeshApprox.size(); }
	
private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox> mSubMeshApprox;
	std::vector<int> mPartIdx;
};