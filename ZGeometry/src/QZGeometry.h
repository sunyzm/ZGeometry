#pragma once
#include "ui_ZGeometry.h"
#include <vector>
#include <string>
#include <QMainWindow>
#include <QSignalMapper>
#include <ZGeom/Mesh.h>
#include <ZGeom/Color.h>
#include "MeshHelper.h"
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
    void computeVertNormals();
    void computeFaceNormals();
    void computeCurvatures();
    void computeShapeIndex();
    void computeGraphLaplacian();
    void computeGeoLaplacian();
    void computeLaplacian(LaplacianType laplacianType);
	void computeEigenfunction();
    void computeBiharmonicDistField();
    void computeSGW();
	void computeHK();
	void computeHKS();	
	void computeHksFeatures();
	
	void computeGeodesics();
	void computeHeatTransfer();

    void regionByDistanceField();

    /* interact with shapeEditor */
	void continuousApprox1(int level);
	void continuousApprox2(int level);
	void continuousApprox3(int level);
	void continuousApprox4(int level);
	void visualizeCompression(int selctedApprox, int coordIdx);
	void displaySignature(QString sigName);
	void displayFeature(QString featureName);
    void displayLine(QString lineName);
	void updateMenuDisplaySignature();
	void updateMenuDisplayFeatures();
    void updateMenuDisplayLines();
    void updateUI();
	void updateSignatureMin(int sMin);
	void updateSignatureMax(int sMax);
	void displayNeighborVertices();
    void computeHoleNeighbors();
	void displayDiffPosition();
	
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
	void toggleShowBoundingBox(bool show = false);
	void toggleShowLines(bool show = false);
	void toggleDrawMatching(bool show = false);
	void toggleShowMatchingLines(bool show = false);
	void toggleDrawRegistration(bool show = false);
    void toggleShowHoles(bool show = false);
    void toggleShowHoleBoundary(bool show = false);
    void toggleShowHoleHollow(bool show = false);
    void toggleShowHoleErrors(bool show = false);
    void toggleShowSurrounding(bool show = false);

    void setColor();
    void setLegend();
	
	void clone();
	void revertCoord();
	void switchToNextCoordinate();
    void switchToPrevCoordinate();
    void deleteCurrentCoordinate();
	void reconstructMHB();
	void deformSimple();
	void deformLaplace();
	void deformLaplace2();
	void deformBiLaplace();
	void deformMixedLaplace();
	void diffusionFlow();
    void fillHoles();

    /* hole filling related */
    void generateHoles();
    void generateRingHoles();
    void generateBandHole();
    void autoGenerateHoles();
    void ignoreOuterBoundary();
    void degradeHoles();
    void cutHoles();
    void cutToSelected();
    void triangulateHoles();
    void refineHoles();
    void refineHoles2();
    void refineHolesByVertNum();
    void copyMeshWithHoles();
    void evaluateCurrentInpainting();    
    void evaluateInpainting2();
    void fairHoleLeastSquares();
    void fairHoleThinPlateEnergy();
    void fairHoleL1LS();
    void smoothingHoleDLRS();

    void curveSignature();

    void switchToNextMesh();
    void switchToPreviousMesh();
    void removeCurrentMesh();

	void registerAutomatic();
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

    /* experiments */
    void doExperiment1();

private:	
	/* methods */
	void makeConnections();
	void loadInitialMeshes(const std::string& initial_mesh_list);
    void loadMesh(std::string mesh_filename, int obj, bool scale_to_unit = false);
	void registerPreprocess();
	void keyPressEvent(QKeyEvent *event);
	void repeatOperation();	// repeat previous operation
	void updateReferenceMove(int obj);
	bool laplacianRequireDecompose(int obj, int nEigVec, LaplacianType laplacianType);
	void decomposeSingleLaplacian(int obj, int nEigVec, LaplacianType laplacianType = CotFormula);
	void allocateStorage(int newMeshCount);
    void updateSignature();

	/* helper functions */
    CMesh* getMesh(int i) { return mMeshHelper[i].getMesh(); }
	bool isMeshSelected(int obj);
	void computeFunctionMaps(int num);
	double parameterFromSlider(double sDefault, double sMin, double sMax, bool verbose = false);
	void addColorSignature(int obj, const std::vector<double>& vVals, const std::string& sigName); 
	/* fields */
	Ui::ZGeometryClass	  ui;
	QSignalMapper*		  m_signatureSignalMapper;
	QSignalMapper*        m_featureSignalMapper;
    QSignalMapper*        m_linesSignalMapper;
	std::vector<QAction*> m_actionDisplaySignatures;
	std::vector<QAction*> m_actionDisplayFeatures;
    std::vector<QAction*> m_actionDisplayLines;
	QLabel				  mStatusLabel;

	std::vector<MeshHelper>	            mMeshHelper;
	std::vector<RenderSettings>			mRenderManagers;
	ShapeMatcher                        mShapeMatcher;
	ShapeEditor	                        mShapeEditor;

	DeformType			    mDeformType;
	LaplacianType           active_lap_type;

	int					mMeshCount;
	int					mObjInFocus;	
	int					mCurrentBasisScale;
	int					mCommonParameter;
	int                 mSelectedApprox, mCoordIdx;
	double              mDiffMax;

	enum { 
		Compute_HKS, Compute_HK, Compute_HKS_Feature,
		Compute_MHWS, Compute_MHW, 
		Compute_SGWS, Compute_SGW,
		Compute_Eig_Func, Compute_Biharmonic, 
		Compute_Dict_Atom,
		Compute_Geodesics, Compute_Heat,
		None
	} mLastOperation;

};
