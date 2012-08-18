#include "manifoldwavelets.h"
#include <QtDebug>

QManifoldWavelets::QManifoldWavelets(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	statusBarMsg = "Computation and visualization for manifold harmonics";
	ui.statusBar->showMessage(statusBarMsg);
	ui.consoleOutput->insertPlainText("Test output\nSecond line");

	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));

}

QManifoldWavelets::~QManifoldWavelets()
{
}
