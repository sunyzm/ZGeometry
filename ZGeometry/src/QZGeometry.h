#pragma once
#include "ui_ZGeometry.h"
#include <vector>
#include <string>
#include <QtWidgets/QMainWindow>
#include <QSignalMapper>
#include <engine.h>
#include <ZMesh/ZMesh.h>
#include <ZUtil/SimpleConfigLoader.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include "OutputHelper.h"
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"
#include "ShapeMatcher.h"
#include "ShapeEditor.h"

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

	void computeSimilarityMap(int simType);
	void computeLaplacian(int laplacianType);
	void computeHK();
	void computeHKS();	
	void computeHKSFeatures();
	void computeMHW();
	void computeMHWS();
	void computeMHWFeatures();
	void computeSGW();
	void computeSGWS();
	void computeSGWSFeatures();
	void computeBiharmonic();
	void computeEigenfunction();
	void computeCurvatureMean();
	void computeCurvatureGauss();

	void displaySignature(int signatureID);
	void displayExperimental();
	void displayNeighborVertices();
	void displayDiffPosition();

	void setDisplayPointCloud();
	void setDisplayWireframe();
	void setDisplayMesh();
	
	void selectObject(int index);
	void setRefPoint1(int vn);
	void setRefPoint2(int vn);
	void setCommonParameter(int p);
	
	void toggleShowSignature(bool show = false);
	void toggleShowRefPoint(bool show = false);
	void toggleShowColorLegend(bool show = false);
	void toggleShowFeatures(bool show = false);
	void toggleDrawMatching(bool show = false);
	void toggleShowMatchingLines(bool show = false);
	void toggleDrawRegistration(bool show = false);
	
	void clone();
	void deformSimple();
	void deformLaplace();
	void deformSGW();
	void reconstructMHB();
	void reconstructSGW();
	void filterExperimental();

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

private:	
	/* methods */
	void makeConnections();
	void loadInitialMeshes(const std::string& initial_mesh_list);
	void registerPreprocess();
	void keyPressEvent(QKeyEvent *event);
	void repeatOperation();	// repeat previous operation
	void updateReferenceMove(int obj);
	void constructLaplacians(MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	void decomposeLaplacians(MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	bool laplacianRequireDecompose(int obj, int nEigVec, MeshLaplacian::LaplacianType laplacianType) const;
	void decomposeSingleLaplacian(int obj, int nEigVec, MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	void allocateStorage(int newMeshCount);

	/* helper functions */
	void evalDistance();
	void computeFunctionMaps(int num);
	void verifyAreas() const;

private:	
	/* fields */
	Ui::ZGeometryClass	ui;
	ZGeom::MatlabEngineWrapper mEngineWrapper;	
	DeformType			deformType;
	int					mObjInFocus;	
	int					mMeshCount;
	int					mCommonParameter;

	std::vector<CMesh*>	                    mMeshes;
	std::vector<DifferentialMeshProcessor*>	mProcessors;
	std::vector<RenderSettings*>			mRenderManagers;
	ShapeMatcher					mShapeMatcher;

	struct {int xMove, yMove, zMove; } refMove;
	enum {Compute_HKS, Compute_HK, 
		  Compute_MHWS, Compute_MHW, 
		  Compute_SGWS, Compute_SGW,
		  Compute_EIG_FUNC,
		  None} current_operation;

	QSignalMapper*		  laplacianSignalMapper;	
	QSignalMapper*		  simlaritySignalMapper;	
	QSignalMapper*		  signatureSignalMapper;
	std::vector<QAction*> m_actionComputeLaplacians;
	std::vector<QAction*> m_actionComputeSimilarities;
	std::vector<QAction*> m_actionDisplaySignatures;

	/*---- static members as constant parameters ----*/
	static int DEFAULT_EIGEN_SIZE;
	static int DEFAULT_DEFORM_RING; 
	static int LOAD_MHB_CACHE;
	static double MIN_HK_TIMESCALE;
	static double DEFUALT_HK_TIMESCALE;
	static double MAX_HK_TIMESCALE;
	static double PARAMETER_SLIDER_CENTER;
	static double DR_THRESH_INCREMENT;
};

