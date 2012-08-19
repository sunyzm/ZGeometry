#include "manifoldwavelets.h"

OutputHelper qout;

QManifoldWavelets::QManifoldWavelets(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	ui.centralWidget->setLayout(ui.mainLayout);
	statusBarMsg = "For computation and visualization of manifold wavelet";
	ui.statusBar->showMessage(statusBarMsg);

	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	
	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);
	qout.output("Welcome!", OUT2CONSOLE);
	qout.output("For computation and visualization of manifold wavelet", OUT2STATUS);
}

QManifoldWavelets::~QManifoldWavelets()
{
}
