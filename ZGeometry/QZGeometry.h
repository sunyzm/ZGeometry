#pragma once
#include "ui_ZGeometry.h"
#include <QtGui/QMainWindow>
#include <ZMesh.h>
#include <engine.h>
#include <vector>
#include "OutputHelper.h"
#include "WaveletMeshProcessor.h"

class QZGeometryWindow : public QMainWindow
{
	Q_OBJECT

public:	// methods
	friend OutputHelper;
	QZGeometryWindow(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QZGeometryWindow();
	bool initialize();

private slots:
	void computeLaplacian(int obj = 0);
	void computeSGWSFeatures();
	void displayEigenfunction();
	void displayMexicanHatWavelet1();
	void displayMexicanHatWavelet2();
	void displayCurvatureMean();
	void displayCurvatureGauss();
	void displayExperimental();
	void displayPointCloud();
	void displayWireframe();
	void displayMesh();
	void displayDiffPosition();
	void selectVertex1(int vn);
	void setShowSignature();
	void setShowRefPoint();
	void setShowColorLegend();
	void setShowFeatures();
	void selectObject(int index);
	void clone();
	void deformExperimental();
	void filterExperimental();
private:	// methods
	void makeConnections();
	void keyPressEvent(QKeyEvent *event);
	void updateReferenceMove();
private:	// attributes
	Ui::ZGeometryClass ui;
//	QString statusBarMsg;
	int totalShapeNum;
	int objSelect;

	Engine *m_ep;
	CMesh mesh1, mesh2;
	bool selected[2];
	WaveletMeshProcessor vMP[2];

	struct {int xMove, yMove, zMove; } refMove;
};
