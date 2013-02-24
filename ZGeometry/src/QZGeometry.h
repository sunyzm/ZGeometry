#pragma once
#include "ui_ZGeometry.h"
#include <QtGui/QMainWindow>
#include <ZMesh.h>
#include <engine.h>
#include <vector>
#include "OutputHelper.h"
#include "DifferentialMeshProcessor.h"

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
	void computeLaplacian(int obj = 0); 
	void computeSGW();
	void computeSGWSFeatures();
	void displayEigenfunction();
	void displayMexicanHatWavelet1();
	void displayMexicanHatWavelet2();
	void displayCurvatureMean();
	void displayCurvatureGauss();
	void displayExperimental();
	void displayNeighborVertices();
	void setDisplayPointCloud();
	void setDisplayWireframe();
	void setDisplayMesh();
	void displayDiffPosition();
	void setRefPoint1(int vn);
	void setCommonParameter(int p);
	void toggleShowSignature(bool show = false);
	void toggleShowRefPoint(bool show = false);
	void toggleShowColorLegend(bool show = false);
	void toggleShowFeatures(bool show = false);
	void selectObject(int index);

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
	void updateReferenceMove();

private:	// attributes
	Ui::ZGeometryClass ui;
	DeformType deformType;
//	QString statusBarMsg;
	int totalShapeNum;
	int objSelect;

	Engine *m_ep;
	CMesh mesh1, mesh2;
	bool selected[2];
	DifferentialMeshProcessor vMP[2];

	struct {int xMove, yMove, zMove; } refMove;
	int m_commonParameter;
};

