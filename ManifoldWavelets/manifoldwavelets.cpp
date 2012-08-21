#include <fstream>
#include <string>
#include <vector>
#include <QtGui/QMessageBox>
#include <QTime>
#include "manifoldwavelets.h"
#define DEFAULT_EIGEN_SIZE 300

using namespace std;

extern OutputHelper qout;
extern QString qformat;
const char* g_meshListName = "meshfiles.cfg";
int g_objSelect = 0;	//-1 means all object; -2 means no object

QManifoldWavelets::QManifoldWavelets(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	m_ep = NULL;
	totalShapeNum = 1;

	ui.setupUi(this);
	ui.centralWidget->setLayout(ui.mainLayout);

	this->makeConnections();

	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);
	qout.output("******** Welcome! ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);
	qout.output("For computation and visualization of manifold wavelet", OUT_STATUS);

}

QManifoldWavelets::~QManifoldWavelets()
{
	engClose(m_ep);
}

void QManifoldWavelets::makeConnections()
{	
	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
////////	Compute	////////
	QObject::connect(ui.actionCompute_Laplacian, SIGNAL(triggered()), this, SLOT(computeLaplace()));
////////	Display	////////
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(displayEigenfunction()));
	
}

bool QManifoldWavelets::initialize()
{
	if (!(m_ep = engOpen("\0")))
	{
		qout.output("Can't start MATLAB engine", OUT_MSGBOX);
		return false;
	}
	else qout.output("Matlab engine opened!", OUT_CONSOLE);

////////    load meshes    ////////
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

////////    load meshes    ////////
	if (!mesh1.Load(vMeshFiles.front()))
	{
		qout.output("Cannot open " + vMeshFiles.front(), OUT_MSGBOX);
		return false;
	}	
	mesh1.scaleEdgeLenToUnit();
	mesh1.gatherStatistics();
	qout.output("Load mesh: " + QString(mesh1.m_meshName.c_str()) + "    Size: " + QString::number(mesh1.getVerticesNum()));

	vMP[0].mesh = &mesh1;
	vMP[0].ep = m_ep;
	vMP[0].m_size = mesh1.getVerticesNum();
	ui.glMeshWidget->addMesh(&vMP[0]);
	
	ui.glMeshWidget->fieldView(mesh1.m_Center, mesh1.m_bBox);
	return true;
}

void QManifoldWavelets::computeLaplace()
{
	QTime timer;
	timer.start();

	string pathMHB = "output/" + vMP[0].mesh->m_meshName + ".mhb";
	ifstream ifs(pathMHB.c_str());
	if (ifs)
	{
		ifs.close();
		vMP[0].readMHB(pathMHB);
	}
	else 
	{
		ifs.close();
		vMP[0].decomposeLaplacian(DEFAULT_EIGEN_SIZE);
		vMP[0].writeMHB(pathMHB);
		qout.output("MHB saved to " + pathMHB);

		string pathEVL = "output/" + vMP[0].mesh->m_meshName + ".evl";	//eigen values
		ofstream ofs(pathEVL.c_str(), ios::trunc);
		for (vector<ManifoldBasis>::iterator iter = vMP[0].mhb.m_func.begin(); iter != vMP[0].mhb.m_func.end(); ++iter)
		{
			ofs << iter->m_val << endl;
		}
		ofs.close();
	}
	
	qout.output("Laplacian decomposed in " + QString::number(timer.elapsed()/1000.0) + "sec");
	qout.output(qformat.sprintf("--Min Eig Val: %f, Max Eig Val: %f", vMP[0].mhb.m_func.front().m_val, vMP[0].mhb.m_func.back().m_val));
}

void QManifoldWavelets::displayEigenfunction()
{
	MeshProcessor& mp = vMP[g_objSelect];
	
	mp.normalizeFrom(mp.mhb.m_func[1].m_vec);
	
	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show eigenfunction 1");
}
