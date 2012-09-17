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
	QObject::connect(ui.actionComputeLaplacian, SIGNAL(triggered()), this, SLOT(computeLaplace()));

////////	Display	////////
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(displayEigenfunction()));
	QObject::connect(ui.actionMexicanHatWavelet1, SIGNAL(triggered()), this, SLOT(displayMexicanHatWavelet1()));
	QObject::connect(ui.actionMexicanHatWavelet2, SIGNAL(triggered()), this, SLOT(displayMexicanHatWavelet2()));
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), ui.horizontalSlider1, SLOT(setValue(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), ui.spinBox1, SLOT(setValue(int)));
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), this, SLOT(selectVertex1(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), this, SLOT(selectVertex1(int)));
	QObject::connect(ui.actionExperimental, SIGNAL(triggered()), this, SLOT(experimental()));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered()), this, SLOT(setShowRefPoint()));
	QObject::connect(ui.actionMeanCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureMean()));
	QObject::connect(ui.actionGaussCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureGauss()));
}

bool QManifoldWavelets::initialize()
{
	if (!(m_ep = engOpen("\0")))
	{
		qout.output("Can't start MATLAB engine", OUT_MSGBOX);
		return false;
	}
	else qout.output("Matlab engine opened!", OUT_CONSOLE);
//	engEvalString(m_ep, "addpath('./MATLAB')");   //need to figure out

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
	vMP[0].init(&mesh1, m_ep);
	ui.glMeshWidget->addMesh(&vMP[0]);
	
	ui.glMeshWidget->fieldView(mesh1.m_Center, mesh1.m_bBox);

	ui.spinBox1->setMinimum(0);
	ui.spinBox1->setMaximum(mesh1.getVerticesNum()-1);
	ui.horizontalSlider1->setMinimum(0);
	ui.horizontalSlider1->setMaximum(mesh1.getVerticesNum()-1);
	ui.spinBox1->setValue(0);

	this->computeLaplace();

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

	ui.actionComputeLaplacian->setChecked(true);
	
	qout.output("Laplacian decomposed in " + QString::number(timer.elapsed()/1000.0) + "sec");
	qout.output(qformat.sprintf("--Min Eig Val: %f, Max Eig Val: %f", vMP[0].mhb.m_func.front().m_val, vMP[0].mhb.m_func.back().m_val));
}

void QManifoldWavelets::selectVertex1( int vn )
{
	vMP[0].pRef = vn;
	ui.glMeshWidget->updateGL();
}

void QManifoldWavelets::displayEigenfunction()
{
	MeshProcessor& mp = vMP[g_objSelect];
	
	mp.normalizeFrom(mp.mhb.m_func[1].m_vec);
	
	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show eigenfunction 1");
}

void QManifoldWavelets::displayMexicanHatWavelet1()
{
	MeshProcessor& mp = vMP[g_objSelect];

	vector<double> vMHW;
	mp.computeMexicanHatWavelet(vMHW, 30, 1);

	mp.normalizeFrom(vMHW);

	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show MHW from point " + QString::number(mp.pRef));

}

void QManifoldWavelets::displayMexicanHatWavelet2()
{
	MeshProcessor& mp = vMP[g_objSelect];

	vector<double> vMHW;
	mp.computeMexicanHatWavelet(vMHW, 30, 2);

	mp.normalizeFrom(vMHW);

	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show MHW from point " + QString::number(mp.pRef));

}

void QManifoldWavelets::experimental()
{
	MeshProcessor& mp = vMP[g_objSelect];

	vector<double> vExp;
	mp.computeExperimentalWavelet(vExp, 30);

	mp.normalizeFrom(vExp);

	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show MHW from point " + QString::number(mp.pRef));

}

void QManifoldWavelets::setShowRefPoint(/*bool checked*/)
{
	bool bChecked = ui.actionShowRefPoint->isChecked();
	ui.glMeshWidget->bShowRefPoint = bChecked;
	ui.actionShowRefPoint->setChecked(bChecked);
	ui.glMeshWidget->updateGL();
}

void QManifoldWavelets::displayCurvatureMean()
{
	MeshProcessor& mp = vMP[g_objSelect];

	vector<double> vCurvature;
	mp.computeCurvature(vCurvature, 0);

	auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
	qout.output("Min curvature: " + QString::number(*mm.first) + "  Max curvature: " + QString::number(*mm.second));

	mp.bandCurveFrom(vCurvature, 0, 1);

	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show Mean Curvature");
}

void QManifoldWavelets::displayCurvatureGauss()
{
	MeshProcessor& mp = vMP[g_objSelect];

	vector<double> vCurvature;
	mp.computeCurvature(vCurvature, 1);

	auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
	qout.output("Min curvature: " + QString::number(*mm.first) + "  Max curvature: " + QString::number(*mm.second));

	mp.bandCurveFrom(vCurvature, -1, 1);

	ui.glMeshWidget->vSettings[g_objSelect].displayType = DisplaySettings::Signature;
	ui.glMeshWidget->updateGL();
	qout.output("Show Mean Curvature");
}
