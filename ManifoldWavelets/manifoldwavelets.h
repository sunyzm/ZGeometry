#ifndef MANIFOLDWAVELETS_H
#define MANIFOLDWAVELETS_H

#include "ui_manifoldwavelets.h"
#include <QtGui/QMainWindow>
#include <mesh/Mesh.h>
#include <util/OutputHelper.h>
#include <engine.h>
#include "WaveletMeshProcessor.h"
#include <vector>

class QManifoldWavelets : public QMainWindow
{
	Q_OBJECT

public:	// methods
	friend OutputHelper;
	QManifoldWavelets(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QManifoldWavelets();
	bool initialize();

private slots:
	void computeLaplacian(int obj = 0);
	void displayEigenfunction();
	void displayMexicanHatWavelet1();
	void displayMexicanHatWavelet2();
	void displayCurvatureMean();
	void displayCurvatureGauss();
	void displayExperimental();
	void displayPointCloud();
	void displayWireframe();
	void showSignature();
	void displayMesh();
	void selectVertex1(int vn);
	void setShowRefPoint();
	void selectObject(int index);
	void reconstruct();
	void deformExperimental();
private:	// methods
	void makeConnections();
	void keyPressEvent(QKeyEvent *event);
	void updateReferenceMove();
private:	// attributes
	Ui::ManifoldWaveletsClass ui;
//	QString statusBarMsg;
	int totalShapeNum;
	int objSelect;

	Engine *m_ep;
	CMesh mesh1, mesh2;
	bool selected[2];
	WaveletMeshProcessor vMP[2];

	struct {int xMove, yMove, zMove; } refMove;
};

#endif // MANIFOLDWAVELETS_H
