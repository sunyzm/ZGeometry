#include <fstream>
#include <string>
#include <vector>
#include <QtGui/QMessageBox>
#include "manifoldwavelets.h"

using namespace std;

extern OutputHelper qout;
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
	vector<string> vMeshFiles;
	while (!meshfiles.eof())
	{
		string meshFileName;
		getline(meshfiles, meshFileName);
		if (meshFileName == "") continue;
		if (meshFileName[0] == '#') continue;
		else vMeshFiles.push_back(meshFileName);
	}
	meshfiles.close();
	
	if (!mesh1.Load(vMeshFiles.front()))
	{
		qout.output("Cannot open " + vMeshFiles.front(), OUT_MSGBOX);
		return false;
	}	
	qout.output("Load mesh: " + QString(mesh1.m_meshName.c_str()) + "    Size: " + QString::number(mesh1.getVerticesNum()));

	mesh1.scaleEdgeLenToUnit();
	mesh1.gatherStatistics();

	ui.glMeshWidget->setMesh(&mesh1, 0);
	ui.glMeshWidget->fieldView(mesh1.m_Center, mesh1.m_bBox);

	return true;
}
