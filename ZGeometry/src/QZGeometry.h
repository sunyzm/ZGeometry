#pragma once
#include "ui_ZGeometry.h"
#include <vector>
#include <string>
#include <QMainWindow>
#include <QSignalMapper>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"
#include "ShapeMatcher.h"
#include "ShapeEditor.h"
#include "global.h"

class QZGeometryWindow : public QMainWindow
{
	Q_OBJECT

public:
	QZGeometryWindow(QWidget *parent = 0,  Qt::WindowFlags flags = 0);
	~QZGeometryWindow();
	bool initialize(const std::string& mesh_list_name);

signals:
	void refPoint1Changed(int n);

private slots:
	void addMesh();
	void saveSignature();
	void saveMatchingResult();
	void loadMatchingResult();

	void setTaskRegistration();
	void setTaskEditing();

	void setEditModeMove();
	void setEditModePick();
	void setEditModeDrag();

	void clearHandles();
	void captureGL();
	void captureGLAs();

	void computeSimilarityMap(int simType);
	void computeLaplacian(int laplacianType);
	void computeEigenfunction();
	void computeEditBasis();
	void computeDictAtom();
	void displayBasis(int idx);
	void computeHK();
	void computeHKS();	
	void computeHKSFeatures();
	void computeMHW();
	void computeMHWS();
	void computeMHWFeatures();
	void computeSGW();
	void computeSGWSFeatures();
	void computeBiharmonic();
	void computeCurvatureMean();
	void computeCurvatureGauss();
	void computeGeodesics();
	void computeHeatTransfer();

	void continuousApprox1(int level);
	void continuousApprox2(int level);
	void continuousApprox3(int level);

	void displaySignature(QString sigName );
	void setSignatureMode(const QString& sigModeName);
	void updateSignatureMin(int sMin);
	void updateSignatureMax(int sMax);
	void displayNeighborVertices();
	void displayDiffPosition();

	void setDisplayPointCloud();
	void setDisplayWireframe();
	void setDisplayMesh();
	
	void selectObject(int index);
	void setRefPoint1(int vn);
	void setRefPoint2(int vn);
	void setCommonParameter(int p);
	void setFeaturePointSize(int v);
	void setLaplacianType(const QString& laplacianTypeName);

	void toggleShowSignature(bool show = false);
	void toggleShowRefPoint(bool show = false);
	void toggleShowColorLegend(bool show = false);
	void toggleShowFeatures(bool show = false);
	void toggleDrawMatching(bool show = false);
	void toggleShowMatchingLines(bool show = false);
	void toggleDrawRegistration(bool show = false);
	
	void clone();
	void revert();
	void nextCoordinate();
	void nextBasisScale();
	void addNoise();
	void reconstructMHB();
	void reconstructSGW();
	void deformSimple();
	void deformLaplace();
	void deformLaplace2();
	void deformBiLaplace();
	void deformMixedLaplace();
	void deformSGW();
	void diffusionFlow();
	void runTests();

	void registerAutomatic();
	void buildHierarchy();
	void detectFeatures();
	void matchFeatures();
	void registerStep();
	void registerFull();
	void registerTest();
	void showFiner();		// lower level
	void showCoarser();		// higher level

	void updateDisplaySignatureMenu();
	void openOutputLocation();
	void resizeApproxSlider(int slider, int newSize);

private:	
	/* methods */
	void makeConnections();
	void loadInitialMeshes(const std::string& initial_mesh_list);
	void registerPreprocess();
	void keyPressEvent(QKeyEvent *event);
	void repeatOperation();	// repeat previous operation
	void updateReferenceMove(int obj);
	void constructLaplacians(LaplacianType laplacianType = CotFormula);
	void decomposeLaplacians(LaplacianType laplacianType = CotFormula);
	bool laplacianRequireDecompose(int obj, int nEigVec, LaplacianType laplacianType) const;
	void decomposeSingleLaplacian(int obj, int nEigVec, LaplacianType laplacianType = CotFormula);
	void allocateStorage(int newMeshCount);
	void updateSignature(SignatureMode smode);
	/* helper functions */
	void evalDistance();
	void computeFunctionMaps(int num);
	void verifyAreas() const;
	void addColorSignature(int obj, const std::vector<double>& vVals, const std::string& sigName); 
	void updateSegmentationColors(int obj, const std::vector<int>& vPartIdx, const Palette& partPalette);
	double parameterFromSlider(double sDefault, double sMin, double sMax, bool verbose = false);
	void signaturesToColors(const std::vector<double>& vOriSig, std::vector<ZGeom::Colorf>& vColors, SignatureMode smode = SignatureMode::SM_Normalized);

	/* fields */
	Ui::ZGeometryClass	  ui;
	QSignalMapper*		  m_laplacianSignalMapper;	
	QSignalMapper*		  m_simlaritySignalMapper;	
	QSignalMapper*		  m_signatureSignalMapper;
	std::vector<QAction*> m_actionComputeLaplacians;
	std::vector<QAction*> m_actionComputeSimilarities;
	std::vector<QAction*> m_actionDisplaySignatures;
	QLabel mStatusLabel1;

	std::vector<CMesh*>	                    mMeshes;
	std::vector<DifferentialMeshProcessor*>	mProcessors;
	std::vector<RenderSettings*>			mRenderManagers;
	ShapeMatcher mShapeMatcher;
	ShapeEditor	 mShapeEditor;

	DeformType			mDeformType;
	LaplacianType       mActiveLalacian;
	SignatureMode		mSignatureMode;
	int					mObjInFocus;	
	int					mMeshCount;
	int					mCurrentBasisScale;
	int					mCommonParameter;

	enum {Compute_HKS, Compute_HK, 
		  Compute_MHWS, Compute_MHW, 
		  Compute_SGWS, Compute_SGW,
		  Compute_Eig_Func, Compute_Biharmonic, 
		  Compute_Edit_Basis, Compute_Dict_Atom,
		  Compute_Geodesics, Compute_Heat,
		  None
	} mLastOperation;
};
