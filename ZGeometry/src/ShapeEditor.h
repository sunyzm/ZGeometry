#pragma once
#include <QObject>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"
#include "GeometryApproximation.h"
#include "Palette.h"
#include "global.h"

class ShapeEditor : public QObject
{
	Q_OBJECT

public:
	friend class QZGeometryWindow;
	ShapeEditor() : mMesh(nullptr), mProcessor(nullptr) {}
	void init(DifferentialMeshProcessor* processor);
	void runTests();
	void revertCoordinates();
	void changeCoordinates(int coordID);
	void nextCoordinates();
	void continuousReconstruct(int selectedApprox, int atomCount);
	const MeshCoordinates& getOldMeshCoord() const { return mOriginalCoord; }
	const Palette& getPalette() const { return mSegmentPalette; }

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
	void signatureComputed(QString sigName);
	void coordinateSelected(int selectedApprox, int coordIdx);

private:
	void addColorSignature(const std::string& colorSigName, const std::vector<ZGeom::Colorf>& vColors);

	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;
	void reconstructSpectralWavelet();
	
	void testSparseCompression();	// test compression performance using graph Laplacian basis
	void sparseDecompositionTest();  // test shape decomposition via sparse coding
	void sparseDecompositionTest2(); // test shape decomposition and separation
	void testArtificialShapeMCA();
	void testShapeMCA();
	void testDictionaryForDecomposition();

	void sparseFeatureFindingTest1();  // test feature finding and correspondence using cot formula basis
	void spectrumTest1();
	
	
	void evaluateApproximation(const MeshCoordinates& newCoord, const std::string leadText);
	void updateEditBasis(const std::vector<ZGeom::VecNd>& vAtoms, const std::vector<int>& vSelectedIdx);
	void computeApproximations(const std::vector<ZGeom::VecNd>& vAtoms, 
		                       ZGeom::SparseCoding* vApproxCoeff[3], 
							   int nReconstruct, 
							   std::vector<MeshCoordinates>& continuousCoords, 
							   MeshCoordinates& finalCoord);
	void setStoredCoordinates(const MeshCoordinates& newCoord, int storeIdx);
	MeshCoordinates& getStoredCoordinate(int idx);
	const MeshCoordinates& getApproximateCoordinate(int selctedApprox, int coordIdx) {
		return mContReconstructCoords[selctedApprox][coordIdx];
	}

// private fields
	CMesh* mMesh;	
	DifferentialMeshProcessor* mProcessor;
	ShapeApprox mShapeApprox;
	Palette mSegmentPalette;

	std::vector<MeshCoordinates> mContReconstructCoords[5];
	int mCurCoordID;
	MeshCoordinates mOriginalCoord;
	std::vector<MeshCoordinates> mStoredCoordinates;
	
	std::vector<ZGeom::VecNd> mEditBasis;	
	std::vector<ZGeom::VecNd> mAtoms;
	int mTotalScales;	
};
