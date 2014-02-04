#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"

enum DictionaryType {DT_Fourier, DT_SGW1, DT_MixedFourierSGW};
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
	void doEigenDecomposition(int eigenCount);
	void findSparseRepresentation(DictionaryType dictType, SparseApproxMethod codingMethod, int codingSize);
	void sparseReconstruction(int reconstructSize);
	void sparseReconstructionStepping(int totalSteps, std::vector<MeshCoordinates>& contCoords);
	void integrateSubmeshApproximation(MeshCoordinates& integratedApproxCoord);
	const MeshCoordinates& getApproxCoord() const { return mApproxCoord; }
	int partitionCount() const { return mSubMeshApprox.size(); }

private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox> mSubMeshApprox;
	std::vector<int> mPartIdx;
	MeshCoordinates mApproxCoord;
};