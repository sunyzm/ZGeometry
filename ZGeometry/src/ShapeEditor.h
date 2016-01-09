#pragma once
#include <QObject>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "MeshHelper.h"
#include "geometry_approximation.h"
#include "Palette.h"
#include "global.h"


class ShapeEditor : public QObject
{
	Q_OBJECT

public:
	friend class QZGeometryWindow;
    static std::string strOriginalCoord;

	ShapeEditor() : mMesh(nullptr), mMeshHelper(nullptr) {}
	void init(MeshHelper& processor);
	void runTests();
    void resetStoredCoordinates();
	void revertCoordinates();
	void changeCoordinates(int coordID);
    void changeCoordinates(std::string coordName);
    void setStoredCoordinates(const MeshCoordinates& newCoord, int storeIdx, std::string coordName = "");
    void addCoordinate(const MeshCoordinates &newCoord, std::string coordName);
    MeshCoordinates getStoredCoordinate(int idx);
	const MeshCoordinates& getOldMeshCoord() const { return mStoredCoordinates[0].second; }
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
    
    void fillHoles(bool skipExternalBoundary);
    void fillBoundedHole(const std::vector<int>& boundaryLoopEdges);
    void visualizeBoundaries();
    void holeFairing();
    void holeFairingFourierOMP();
    void holeFairingFourierLARS();
    void holeFairingFourierLS();
    void holeFairingLeastSquare();
    void holeEstimateCurvature();
    void holeEstimateNormals();

    void generateNoise(const std::vector<int>& selectedVerts, double sigma = 0.02);
    void inpaintHolesLARS(const std::vector<int>& selectedVerts, double eps = 1e-4);
    void inpaintHolesL1LS(const std::vector<int>& selctedVerts, double lambda = 1e-3, double tol = 1e-3);

signals:
	void approxStepsChanged(int index, int newSize);
	void coordinateSelected(int selectedApprox, int coordIdx);
    void meshSignatureAdded();
    void meshPointFeatureChanged();    
    void meshLineFeatureChanged();
    void showSignature(QString sigName);

private:
	void addColorSignature(const std::string& colorSigName, const std::vector<ZGeom::Colorf>& vColors);
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<ZGeom::Vec3d>& anchorPos) const;
	void testSparseCompression();	// test compression performance using graph Laplacian basis
	void testSparseDecomposition();  // test shape decomposition via sparse coding
	void testSparseDecomposition2(); // test shape decomposition and separation
	void testArtificialShapeMCA();	// MHB signal + SGW noise / MHB & SGW dictionary
	void testDictionaryForDecomposition();
	void testSparseFeatureFinding();  // test feature finding and correspondence using cot formula basis
	void testSparseInpainting();
    void testDenoisingDLRS();
    void testDictionaryCoherence();
    void testWaveletAnalysis();
    void testWaveletComputation();
    void testSurfaceInpainting();

	void evaluateApproximation(const MeshCoordinates& newCoord, const std::string leadText);
	void updateEditBasis(const std::vector<ZGeom::VecNd>& vAtoms, const std::vector<int>& vSelectedIdx);
	void computeApproximations(const std::vector<ZGeom::VecNd>& vAtoms, 
		                       ZGeom::SparseCoding* vApproxCoeff[3], 
							   int nReconstruct, 
							   std::vector<MeshCoordinates>& continuousCoords, 
							   MeshCoordinates& finalCoord);
	const MeshCoordinates& getApproximateCoordinate(int selctedApprox, int coordIdx) { return mContReconstructCoords[selctedApprox][coordIdx]; }

    /* private fields */
	CMesh* mMesh;	
	MeshHelper* mMeshHelper;
	ShapeApprox mShapeApprox;
	Palette mSegmentPalette;

	std::vector<MeshCoordinates> mContReconstructCoords[5];
	int mCurCoordID;
    std::vector<std::pair<std::string, MeshCoordinates>> mStoredCoordinates;
	
	std::vector<ZGeom::VecNd> mEditBasis;	
	int mTotalScales;	

public:
    std::vector<ZGeom::MeshRegion> filled_boundaries;

};

