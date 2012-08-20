#include <fstream>
#include <string>
#include <vector>
#include <QtGui/QMessageBox>
#include "manifoldwavelets.h"

using namespace std;

OutputHelper qout;
const char* g_meshListName = "meshfiles.cfg";

QManifoldWavelets::QManifoldWavelets(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	ui.centralWidget->setLayout(ui.mainLayout);

	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	
	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);
	qout.output("Welcome!", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);
	qout.output("For computation and visualization of manifold wavelet", OUT_STATUS);
}

QManifoldWavelets::~QManifoldWavelets()
{
}

bool QManifoldWavelets::initialize()
{
	//// load meshes ////
	ifstream meshfiles(g_meshListName);
	if (!meshfiles)
	{
		qout.output("Cannot open " + QString(g_meshListName), OUT_MSGBOX);
		return false;
	}
	vector<string> meshFileList;
	while (!meshfiles.eof())
	{
		string meshFileName;
		getline(meshfiles, meshFileName);
		if (meshFileName == "") continue;
		if (meshFileName[0] == '#') continue;
		else meshFileList.push_back(meshFileName);
	}
	meshfiles.close();
	
	mesh1.Load()

	qout.output(QString(meshFileList.front().c_str()));

	return true;
}
