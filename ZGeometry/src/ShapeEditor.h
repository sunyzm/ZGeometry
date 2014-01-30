#pragma once
#include <QObject>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"
#include "Palette.h"
using ZGeom::VecNd;


enum DictionaryType {DT_Fourier, DT_SGW, DT_MixedFourierSGW};
enum SparseApproxMethod {SA_Truncation, SA_SMP, SA_SOMP};

class Dictionary
{
	friend class SubMeshApprox;
public:
	void resize(int N, int m) 
	{ 
		mAtoms.resize(N); 
		for (VecNd& v : mAtoms) v.resize(m);
	}
	void clear() { mAtoms.clear(); }
	void resize(int N) { mAtoms.resize(N); }
	VecNd& operator[] (int i) { return mAtoms[i]; }
	int atomCount() const { return mAtoms.size(); }
	const std::vector<VecNd>& getAtoms() const { return mAtoms; }

private:
	std::vector<VecNd> mAtoms;
};

class SubMeshApprox
{
	friend class ShapeApprox;
public:
	struct SparseCoeff {
		SparseCoeff() : mIdx(-1), mCoeff(0) {}
		SparseCoeff(int i, double c) : mIdx(i), mCoeff(c) {}
		int mIdx;
		double mCoeff;
	};
	void init() { mMeshProcessor.init(&mSubMesh, &g_engineWrapper); }
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
	Dictionary mDict;
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
	const Palette& getPalette() const { return mSegmentPalette; }
	const MeshCoordinates& getApproxCoord() const { return mApproxCoord; }

private:
	CMesh* mOriginalMesh;	
	std::vector<SubMeshApprox> mSubMeshApprox;
	Palette mSegmentPalette;
	MeshCoordinates mApproxCoord;
};

class ShapeEditor : public QObject
{
	Q_OBJECT

public:
	friend class QZGeometryWindow;
	ShapeEditor() : mMesh(nullptr), mProcessor(nullptr) {}
	void init(DifferentialMeshProcessor* processor);
	void revertCoordinates();
	void changeCoordinates(int coordID);
	void nextCoordinates();
	void continuousReconstruct(int selectedApprox, int atomCount);
	const MeshCoordinates& getOldMeshCoord() const { return mCoords[0]; }

	void addNoise(double phi);

	void fourierReconstruct(int nEig);
	void meanCurvatureFlow(double tMultiplier, int nRepeat = 1);

	void deformSimple();
	void deformLaplacian();
	void deformLaplacian_v2();
	void deformBiLaplacian();
	void deformMixedLaplacian(double ks, double kb);
	void deformThinShell2(double ks, double kb);
	void deformSpectralWavelet();	

signals:
	void approxStepsChanged(int index, int newSize);

private:
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;
	void reconstructSpectralWavelet();
	void evalReconstruct(const MeshCoordinates& newCoord) const;
	
	void monolithicApproximationTest1(bool doWavelet = true);	// use graph Laplacian
	void monolithicApproximationTest2(bool doWavelet = true);	// use CotFormula Laplacian
	void partitionedApproximationTest1();
	void editTest2();
	void evaluateApproximation(const MeshCoordinates& newCoord, const std::string leadText);

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

	std::vector<MeshCoordinates> mContReconstructCoords[3];
	int mCurCoordID;
	std::vector<MeshCoordinates> mCoords;
};