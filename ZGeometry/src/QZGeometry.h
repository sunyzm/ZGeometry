#pragma once
#include "ui_ZGeometry.h"
#include <QtGui/QMainWindow>
#include <ZMesh.h>
#include <engine.h>
#include <vector>
#include "OutputHelper.h"
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"

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
	void setEditModeMove();
	void setEditModePick();
	void setEditModeDrag();

	void computeLaplacian(); 
	void computeHKS();
	void computeHK();
	void computeSGW();
	void computeSGWSFeatures();

	void displayEigenfunction();
	void displayHKS();
	void displayHK();
	void displayMexicanHatWavelet1();
	void displayMexicanHatWavelet2();
	void displayCurvatureMean();
	void displayCurvatureGauss();
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
	
	void clone();
	void deformSimple();
	void deformLaplace();
	void deformSGW();
	void reconstructMHB();
	void reconstructSGW();
	void filterExperimental();

private:	// methods
	void makeConnections();
	void keyPressEvent(QKeyEvent *event);
	void repeatOperation();	// repeat previous operation
	void updateReferenceMove(int obj);

private:	// attributes
	Ui::ZGeometryClass ui;
	DeformType deformType;
//	QString statusBarMsg;
	int totalShapeNum;
	int objSelect;

	Engine *m_ep;
	CMesh mesh1, mesh2;
	bool selected[2];
	bool mesh_valid[2];
	DifferentialMeshProcessor vMP[2];
	RenderSettings vRS[2];
	
	struct {int xMove, yMove, zMove; } refMove;
	int m_commonParameter;
	enum {Compute_HKS, Compute_HK, None} current_operation;
};

