#pragma once
#include "ui_ZGeometry.h"
#include <QtGui/QMainWindow>
#include <ZMesh.h>
#include <engine.h>
#include <vector>
#include "OutputHelper.h"
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"
#include "DiffusionShapeMatcher.h"
#include "SimpleConfigLoader.h"


class QZGeometryWindow : public QMainWindow
{
	Q_OBJECT

public:	// methods
	friend OutputHelper;
	QZGeometryWindow(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QZGeometryWindow();
	bool initialize();
signals:
	void refPoint1Changed(int n);

private slots:
	void addMesh();
	void saveSignature();

	void setTaskRegistration();
	void setTaskEditing();

	void setEditModeMove();
	void setEditModePick();
	void setEditModeDrag();

	void computeLaplacian1(); 
	void computeLaplacian2();
	void computeLaplacian3();
	void computeLaplacian4();
	void computeLaplacian5();
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
	void computeSimilarityMap();
	void computeSimilarityMap2();
	void computeSimilarityMap3();

	void displayEigenfunction();
	void displaySignature(int signatureID);
	void displayHKS();
	void displayHK();
	void displayMHW();
	void displayMHWS();
	void displaySGWS();
	void displayBiharmonic();
	void displayCurvatureMean();
	void displayCurvatureGauss();
	void displayExperimental();
	void displayNeighborVertices();
	void displayDiffPosition();
	void displaySimilarityMap();

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
	void showFiner();		// lower level
	void showCoarser();		// higher level

private:	// methods
	void makeConnections();
	void keyPressEvent(QKeyEvent *event);
	void repeatOperation();	// repeat previous operation
	void updateReferenceMove(int obj);
	void calculateLaplacians(LaplacianType laplacianType = CotFormula);
	void calculateSingleLaplacian(int obj, LaplacianType laplacianType = CotFormula);
	// helper functions
	void evalDistance();
	
private:	// attributes
	Ui::ZGeometryClass ui;
	DeformType deformType;
	int objSelect;

	Engine *m_ep;
	
	int num_meshes;
	std::vector<CMesh> m_mesh;
	std::vector<bool> mesh_valid;
	std::vector<DifferentialMeshProcessor> vMP;
	std::vector<RenderSettings> vRS;

	DiffusionShapeMatcher shapeMatcher;

	struct {int xMove, yMove, zMove; } refMove;
	int m_commonParameter;
	enum {Compute_HKS, Compute_HK, 
		  Compute_MHWS, Compute_MHW, 
		  Compute_SGWS, Compute_SGW, None} current_operation;

	std::vector<QAction*> m_actionLaplacians;

	/*---- static members as constant parameters ----*/
	static int DEFAULT_EIGEN_SIZE;
	static int DEFAULT_DEFORM_RING; 
	static int LOAD_MHB_CACHE;
	static double MIN_HK_TIMESCALE;
	static double DEFUALT_HK_TIMESCALE;
	static double MAX_HK_TIMESCALE;
	static double PARAMETER_SLIDER_CENTER;
	static double DR_THRESH_INCREMENT;
	static double MATCHING_THRESHOLD;
};

