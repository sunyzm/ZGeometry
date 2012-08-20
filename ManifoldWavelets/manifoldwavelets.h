#ifndef MANIFOLDWAVELETS_H
#define MANIFOLDWAVELETS_H

#include "ui_manifoldwavelets.h"
#include <QtGui/QMainWindow>
#include <mesh/Mesh.h>
#include <util/OutputHelper.h>

class QManifoldWavelets : public QMainWindow
{
	Q_OBJECT

public:
	friend OutputHelper;
	QManifoldWavelets(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QManifoldWavelets();

	bool initialize();

private:
	Ui::ManifoldWaveletsClass ui;
//	QString statusBarMsg;

	CMesh mesh1;
	void drawMesh();
};

#endif // MANIFOLDWAVELETS_H
