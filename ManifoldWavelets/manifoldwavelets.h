#ifndef MANIFOLDWAVELETS_H
#define MANIFOLDWAVELETS_H

#include "ui_manifoldwavelets.h"
#include <QtGui/QMainWindow>
#include <mesh/Mesh.h>
#include <util/OutputHelper.h>
#include <engine.h>
#include "MeshProcessor.h"
#include <vector>

class QManifoldWavelets : public QMainWindow
{
	Q_OBJECT

public:
	friend OutputHelper;
	QManifoldWavelets(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QManifoldWavelets();

	bool initialize();
private slots:
	void computeLaplace();
	void displayEigenfunction();
	void makeConnections();
private:
	Ui::ManifoldWaveletsClass ui;
//	QString statusBarMsg;
	int totalShapeNum;
	int objSelect;

	Engine *m_ep;
	CMesh mesh1;
	MeshProcessor vMP[2];
};

#endif // MANIFOLDWAVELETS_H
