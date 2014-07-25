#pragma once
#include "ui_ZGeometry.h"
#include <vector>
#include <string>
#include <QMainWindow>
#include <QSignalMapper>
#include <ZGeom/Mesh.h>
#include <ZGeom/Color.h>
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
	void displayQtVersion();

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
	
	void listMeshAttributes();

	/* computation */
	void computeLaplacian(int laplacianType);
	void computeEigenfunction();
	void computeCurvatureMean();
	void computeCurvatureGauss();
	void computeVertNormals();
	void computeFaceNormals();
	void computeEditBasis();
	void computeDictAtom();
	void computeHK();
	void computeHKS();	
	void computeHKSFeatures();
	void computeMHW();
	void computeMHWS();
	void computeSGW();
	void computeBiharmonic();
	void computeGeodesics();
	void computeHeatTransfer();

	/* interact with shapeEditor */
	void continuousApprox1(int level);
	void continuousApprox2(int level);
	void continuousApprox3(int level);
	void continuousApprox4(int level);
	void visualizeCompression(int selctedApprox, int coordIdx);

	void displaySignature(QString sigName);
	void displayFeature(QString featureName);
	void updateDisplaySignatureMenu();
	void updateDisplayFeatureMenu();
	void setSignatureMode(const QString& sigModeName);
	void updateSignatureMin(int sMin);
	void updateSignatureMax(int sMax);
	void displayNeighborVertices();
	void displayDiffPosition();
	void displayBasis(int idx);
	
	void setDisplayPointCloud();
	void setDisplayWireframe();
	void setDisplayMesh();
	
	void selectObject(int index);
	void setMesh1RefPoint(int vn);
	void setMesh2RefPoint(int vn);
	void setCommonParameter(int p);
	void setFeaturePointSize(int v);
	void setLaplacianType(const QString& laplacianTypeName);

	void toggleShowSignature(bool show = false);
	void toggleShowRefPoint(bool show = false);
	void toggleShowColorLegend(bool show = false);
	void toggleShowFeatures(bool show = false);
	void toggleShowWireframeOverlay(bool show = false);
	void toggleShowVectors(bool show = false);
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

	void openOutputLocation();
	void openSreenshotLocation();
	void resizeApproxSlider(int slider, int newSize);

private:	
	/* methods */
	void makeConnections();
	void loadInitialMeshes(const std::string& initial_mesh_list);
	void registerPreprocess();
	void keyPressEvent(QKeyEvent *event);
	void repeatOperation();	// repeat previous operation
	void updateReferenceMove(int obj);
	bool laplacianRequireDecompose(int obj, int nEigVec, LaplacianType laplacianType) const;
	void decomposeSingleLaplacian(int obj, int nEigVec, LaplacianType laplacianType = CotFormula);
	void allocateStorage(int newMeshCount);
	void updateSignature(SignatureMode smode);

	/* helper functions */
	bool isMeshSelected(int obj);
	void evalDistance();
	void computeFunctionMaps(int num);
	void verifyAreas() const;
	double parameterFromSlider(double sDefault, double sMin, double sMax, bool verbose = false);
	void addColorSignature(int obj, const std::vector<double>& vVals, const std::string& sigName); 
	void signaturesToColors(const std::vector<double>& vOriSig, std::vector<ZGeom::Colorf>& vColors, SignatureMode smode = SignatureMode::SM_Normalized);

	/* fields */
	Ui::ZGeometryClass	  ui;
	QSignalMapper*		  m_laplacianSignalMapper;	
	QSignalMapper*		  m_signatureSignalMapper;
	QSignalMapper*        m_featureSignalMapper;
	std::vector<QAction*> m_actionComputeLaplacians;
	std::vector<QAction*> m_actionDisplaySignatures;
	std::vector<QAction*> m_actionDisplayFeatures;
	QLabel				  mStatusLabel1;

	std::vector<CMesh*>	                    mMeshes;
	std::vector<DifferentialMeshProcessor*>	mProcessors;
	std::vector<RenderSettings*>			mRenderManagers;

	ShapeMatcher mShapeMatcher;
	ShapeEditor	 mShapeEditor;

	DeformType			mDeformType;
	LaplacianType       mActiveLalacian;
	SignatureMode		mSignatureMode;
	ZGeom::ColorMapType mColorMapType;

	int					mMeshCount;
	int					mObjInFocus;	
	int					mCurrentBasisScale;
	int					mCommonParameter;
	int                 mSelectedApprox, mCoordIdx;
	double              mDiffMax;

	enum { 
		Compute_HKS, Compute_HK, 
		Compute_MHWS, Compute_MHW, 
		Compute_SGWS, Compute_SGW,
		Compute_Eig_Func, Compute_Biharmonic, 
		Compute_Edit_Basis, Compute_Dict_Atom,
		Compute_Geodesics, Compute_Heat,
		None
	} mLastOperation;
};
