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
	const MeshCoordinates& getOldMeshCoord() const { return mOriginalCoord; }
	const Palette& getPalette() const { return mSegmentPalette; }

	void continuousReconstruct(int selectedApprox, int atomCount);
    MeshCoordinates getNoisyCoord(double phi);
	void fourierReconstruct(int nEig);
	void meanCurvatureFlow(double tMultiplier, int nRepeat = 1);
	void deformSimple();
	void deformLaplacian();
	void deformLaplacian_v2();
	void deformBiLaplacian();
	void deformMixedLaplacian(double ks, double kb);
	void deformThinShell2(double ks, double kb);
    
    void fillHoles(bool skipExternalBoundary);
    void fillBoundedHole(const std::vector<int>& boundaryLoopEdges);
    void fillHole();
    void holeFairing();
    void visualizeBoundaries();

signals:
	void approxStepsChanged(int index, int newSize);
	void coordinateSelected(int selectedApprox, int coordIdx);
    void meshSignatureAdded();
    void meshPointFeatureChanged();    
    void meshLineFeatureChanged();
    void showSignature(QString sigName);

private:
	void addColorSignature(const std::string& colorSigName, const std::vector<ZGeom::Colorf>& vColors);
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;
	void testSparseCompression();	// test compression performance using graph Laplacian basis
	void testSparseDecomposition();  // test shape decomposition via sparse coding
	void testSparseDecomposition2(); // test shape decomposition and separation
	void testArtificialShapeMCA();	// MHB signal + SGW noise / MHB & SGW dictionary
	void testArtificailShapeMCA2();	// original signal + SGW noise / MHB & SGW dictionary
	void testArtificialShapeMCA3(); // original signal + (unrestricted) artificial noise
	void testDictionaryForDecomposition();
	void testSparseFeatureFinding();  // test feature finding and correspondence using cot formula basis
	void testSparseInpainting();
    void testDenoisingDLRS();
    void testDictionaryCoherence();
    void testWaveletAnalysis();
    void testWaveletComputation();
    void testSurfaceArea();

	void evaluateApproximation(const MeshCoordinates& newCoord, const std::string leadText);
	void updateEditBasis(const std::vector<ZGeom::VecNd>& vAtoms, const std::vector<int>& vSelectedIdx);
	void computeApproximations(const std::vector<ZGeom::VecNd>& vAtoms, 
		                       ZGeom::SparseCoding* vApproxCoeff[3], 
							   int nReconstruct, 
							   std::vector<MeshCoordinates>& continuousCoords, 
							   MeshCoordinates& finalCoord);
	void setStoredCoordinates(const MeshCoordinates& newCoord, int storeIdx);
	MeshCoordinates& getStoredCoordinate(int idx);
	const MeshCoordinates& getApproximateCoordinate(int selctedApprox, int coordIdx) { return mContReconstructCoords[selctedApprox][coordIdx]; }

    /* private fields */
	CMesh* mMesh;	
	DifferentialMeshProcessor* mProcessor;
	ShapeApprox mShapeApprox;
	Palette mSegmentPalette;

	std::vector<MeshCoordinates> mContReconstructCoords[5];
	int mCurCoordID;
	MeshCoordinates mOriginalCoord;
	std::vector<MeshCoordinates> mStoredCoordinates;
	
	std::vector<ZGeom::VecNd> mEditBasis;	
	int mTotalScales;	

    /* fields for boundaries */
    struct FilledHoleVerts {
        std::vector<int> vert_on_boundary;
        std::vector<int> vert_inside;
    };

    std::vector<FilledHoleVerts> filled_boundaries;
};
