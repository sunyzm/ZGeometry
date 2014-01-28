#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"

class Dictionary : public ZGeom::DenseMatrixd
{
public:
	int atomCount() const { return colCount(); }
};

class ShapeApprox
{
public:
	friend class ShapeEditor;	
	enum DictionaryType {Fourier, SGW, MixedFourierSGW};
	enum SparseApproxMethod {Truncation, SMP, SOMP};

	class SubMeshApprox
	{
	public:
		friend class ShapeApprox;
		struct SparseCoeff {
			int mIdx;
			double mCoeff;
		};
		
		void init() { mMeshProcessor.init(&mSubMesh, &gEngineWrapper); }
		int subMeshSize() const { return mSubMesh.vertCount(); }
		const std::vector<int>& mappedIdx() const { return mMappedIdx; }

		void prepareEigenSystem(LaplacianType laplacianType, int mEigenCount);
		void constructDict(DictionaryType dictType);
		void sparseApproximate(SparseApproxMethod approxMethod);

	private:
		CMesh mSubMesh;
		DifferentialMeshProcessor mMeshProcessor;
		std::vector<int> mMappedIdx;
		ZGeom::EigenSystem mEigenSystem;
		Dictionary mDict;
		std::vector<SparseCoeff> mCoding;
		MeshCoordinates mReconstructedCoord;
	};

	ShapeApprox() : mOriginalMesh(NULL) {}
	void init(CMesh* mesh) { mOriginalMesh = mesh; }
	void doSegmentation(int maxSize);
	void doSparseApproximation();
	void evaluate();

private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox> mSubMeshApprox;

};

class ShapeEditor
{
	friend class QZGeometryWindow;
public:
	ShapeEditor() : mMesh(nullptr), mProcessor(nullptr) {}
	void init(DifferentialMeshProcessor* processor);
	void revertCoordinates();
	void changeCoordinates();
	void continuousReconstruct(int selectedApprox, int atomCount);
	const MeshCoordinates& getOldMeshCoord() const { return mCoords[0]; }

	void addNoise(double phi);

	void manifoldHarmonicsReconstruct(int nEig);
	void meanCurvatureFlow(double tMultiplier, int nRepeat = 1);

	void deformSimple();
	void deformLaplacian();
	void deformLaplacian_v2();
	void deformBiLaplacian();
	void deformMixedLaplacian(double ks, double kb);
	void deformThinShell2(double ks, double kb);
	void deformSpectralWavelet();	

private:
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;
	void reconstructSpectralWavelet();
	void evalReconstruct(const MeshCoordinates& newCoord) const;
	
	void reconstructionTest1();
	void reconstructionTest2(bool doWavelet = true);
	void approximationTest1();
	void editTest2();

	void updateEditBasis(const std::vector<ZGeom::VecNd>& vAtoms, const std::vector<int>& vSelectedIdx);
	void computeApproximations(const std::vector<ZGeom::VecNd>& vAtoms, 
		                       ZGeom::FunctionApproximation* vApproxCoeff[3], 
							   int nReconstruct, 
							   std::vector<MeshCoordinates>& continuousCoords, 
							   MeshCoordinates& finalCoord);

	CMesh* mMesh;	
	DifferentialMeshProcessor* mProcessor;
	ShapeApprox mShapeApprox;
	ZGeom::MatlabEngineWrapper* mEngine;
	std::vector<ZGeom::VecNd> mEditBasis;	
	std::vector<ZGeom::VecNd> mAtoms;
	int mTotalScales;

	ZGeom::FunctionApproximation mApproxCoeff[3];
	std::vector<MeshCoordinates> mContReconstructCoords[3];
	std::vector<MeshCoordinates> mCoords;
	int mCoordSelect;
};