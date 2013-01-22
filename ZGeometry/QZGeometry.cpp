#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <QtGui/QMessageBox>
#include <QTime>
#include <GL/glew.h>
#include "QZGeometry.h"

#define DEFAULT_EIGEN_SIZE 300
#define LOAD_MHB_CACHE 1

using namespace std;

extern OutputHelper qout;
extern QString qformat;
const char* g_meshListName = "meshfiles.cfg";
//int g_objSelect = 0;	//-1 means all object; -2 means no object

QZGeometryWindow::QZGeometryWindow(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	m_ep = NULL;
	totalShapeNum = 1;

	ui.setupUi(this);
	ui.centralWidget->setLayout(ui.mainLayout);

	this->makeConnections();

	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);
	
	// initialize GLEW
//	if(glewInit() != GLEW_OK)
	//	throw std::runtime_error("glewInit failed");
	// print out some info about the graphics drivers

// 	ostringstream oss;
// 	oss << "OpenGL version: " << glGetString(GL_VERSION) << std::endl;
// 	oss << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
// 	oss << "Vendor: " << glGetString(GL_VENDOR) << std::endl;
// 	oss << "Renderer: " << glGetString(GL_RENDERER) << std::endl;
//	qout.output(oss.str(), OUT_CONSOLE);
	// make sure OpenGL version 3.2 API is available
// 	if(!GLEW_VERSION_3_2)
// 		throw std::runtime_error("OpenGL 3.2 API is not available.");
	
	qout.output("******** Welcome! ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);
	qout.output("For computation and visualization of manifold wavelet", OUT_STATUS);

	refMove.xMove = refMove.yMove = refMove.zMove = 0;
}

QZGeometryWindow::~QZGeometryWindow()
{
	engClose(m_ep);
}

void QZGeometryWindow::makeConnections()
{	
	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
////////	Compute	////////
	QObject::connect(ui.actionComputeLaplacian, SIGNAL(triggered()), this, SLOT(computeLaplacian()));

////////	Display	////////
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(displayEigenfunction()));
	QObject::connect(ui.actionMexicanHatWavelet1, SIGNAL(triggered()), this, SLOT(displayMexicanHatWavelet1()));
	QObject::connect(ui.actionMexicanHatWavelet2, SIGNAL(triggered()), this, SLOT(displayMexicanHatWavelet2()));
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), ui.horizontalSlider1, SLOT(setValue(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), ui.spinBox1, SLOT(setValue(int)));
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), this, SLOT(selectVertex1(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), this, SLOT(selectVertex1(int)));
	QObject::connect(ui.actionExperimental, SIGNAL(triggered()), this, SLOT(displayExperimental()));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered()), this, SLOT(setShowRefPoint()));
	QObject::connect(ui.actionMeanCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureMean()));
	QObject::connect(ui.actionGaussCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureGauss()));
	QObject::connect(ui.objSelectBox, SIGNAL(activated(int)), this, SLOT(selectObject(int)));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
	QObject::connect(ui.actionDeformExperimental, SIGNAL(triggered()), this, SLOT(deformExperimental()));
	QObject::connect(ui.actionDisplayPointCloud, SIGNAL(triggered()), this, SLOT(displayPointCloud()));
	QObject::connect(ui.actionDisplayWireframe, SIGNAL(triggered()), this, SLOT(displayWireframe()));
	QObject::connect(ui.actionDisplayMesh, SIGNAL(triggered()), this, SLOT(displayMesh()));
	QObject::connect(ui.actionShowSignature, SIGNAL(triggered()), this, SLOT(setShowSignature()));
	QObject::connect(ui.actionSGWSFeatures, SIGNAL(triggered()), this, SLOT(computeSGWSFeatures()));
	QObject::connect(ui.actionFilterExperimental, SIGNAL(triggered()), this, SLOT(filterExperimental()));
	QObject::connect(ui.actionDiffPosition, SIGNAL(triggered()), this, SLOT(displayDiffPosition()));
	QObject::connect(ui.actionColorLegend, SIGNAL(triggered()), this, SLOT(setShowColorLegend()));
}

bool QZGeometryWindow::initialize()
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

	if (!mesh1.Load(vMeshFiles.front()))
	{
		qout.output("Cannot open " + vMeshFiles.front(), OUT_MSGBOX);
		return false;
	}
	mesh1.scaleEdgeLenToUnit();
	mesh1.gatherStatistics();
	qout.output("Load mesh: " + QString(mesh1.m_meshName.c_str()) + "    Size: " + QString::number(mesh1.getVerticesNum()));

/*	
	ifstream ifcoord("output/coordtrans4.dat");
	for (int i = 0; i < mesh1.getVerticesNum(); ++i)
	{
		double x(0),y(0),z(0);
		ifcoord >> x >> y >> z;
		if (i < 5) 
			qout.output("vertex: " + QString::number(x) + "," + QString::number(y) + "," + QString::number(z));
		mesh1.m_pVertex[i].m_vPosition = Vector3D(x,y,z);
	}
//*/	
	
	vMP[0].init(&mesh1, m_ep);
	ui.glMeshWidget->addMesh(&vMP[0]);
	ui.glMeshWidget->fieldView(mesh1.m_Center, mesh1.m_bBox);

	ui.spinBox1->setMinimum(0);
	ui.spinBox1->setMaximum(mesh1.getVerticesNum()-1);
	ui.horizontalSlider1->setMinimum(0);
	ui.horizontalSlider1->setMaximum(mesh1.getVerticesNum()-1);
	ui.spinBox1->setValue(0);	

	selected[0] = ui.glMeshWidget->vSettings[0].selected = true;
	selected[1] = ui.glMeshWidget->vSettings[1].selected = false; 

	this->computeLaplacian(0);
	
	return true;
}

void QZGeometryWindow::computeLaplacian( int obj /*= 0*/ )
{
	QTime timer;
	timer.start();

	string pathMHB = "output/" + vMP[obj].mesh->m_meshName + ".mhb";
	
	ifstream ifs(pathMHB.c_str());
	if (ifs && LOAD_MHB_CACHE)
	{
		vMP[obj].readMHB(pathMHB);
	}
	else 
	{
		vMP[obj].decomposeLaplacian(DEFAULT_EIGEN_SIZE);
//		vMP[obj].decomposeLaplacian(vMP[obj].m_size);
		vMP[obj].writeMHB(pathMHB);
		qout.output("MHB saved to " + pathMHB);

		string pathEVL = "output/" + vMP[obj].mesh->m_meshName + ".evl";	//eigenvalues
		ofstream ofs(pathEVL.c_str(), ios::trunc);
		for (vector<ManifoldBasis>::iterator iter = vMP[obj].mhb.m_func.begin(); iter != vMP[obj].mhb.m_func.end(); ++iter)
		{
			ofs << iter->m_val << endl;
		}
		ofs.close();
	}
	ifs.close();
	ui.actionComputeLaplacian->setChecked(true);
	
	qout.output("Laplacian decomposed in " + QString::number(timer.elapsed()/1000.0) + " (s)");
	qout.output(qformat.sprintf("--Min Eig Val: %f, Max Eig Val: %f", vMP[obj].mhb.m_func.front().m_val, vMP[obj].mhb.m_func.back().m_val));
}

void QZGeometryWindow::selectVertex1( int vn )
{
	vMP[0].pRef = vn;
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	updateReferenceMove();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayEigenfunction()
{
	if (selected[0])
	{
		MeshProcessor& mp = vMP[0];
		mp.normalizeSignatureFrom(mp.mhb.m_func[1].m_vec);
		ui.glMeshWidget->vSettings[0].showColorSignature = true;
		ui.actionShowSignature->setChecked(true);
	}

// 	if (selected[1] && vMP[1].mesh)
// 	{
// 		MeshProcessor& mp = vMP[1];		
// 		mp.normalizeFrom(mp.mhb.m_func[1].m_vec);
// 		ui.glMeshWidget->vSettings[1].showColorSignature = true;
// 	}
	
	ui.glMeshWidget->update();
	qout.output("Show eigenfunction 1");
}

void QZGeometryWindow::displayMexicanHatWavelet1()
{
	if (selected[0])
	{
		WaveletMeshProcessor& mp = vMP[0];

		vector<double> vMHW;
		mp.computeMexicanHatWavelet(vMHW, 30, 1);
		mp.normalizeSignatureFrom(vMHW);

		ui.glMeshWidget->vSettings[0].showColorSignature = true;
		ui.actionShowSignature->setChecked(true);
	}

// 	if (selected[1] && vMP[1].mesh)
// 	{
// 		WaveletMeshProcessor& mp = vMP[1];		
// 
// 		vector<double> vMHW;
// 		mp.computeMexicanHatWavelet(vMHW, 30, 1);
// 		mp.normalizeFrom(vMHW);
// 
// 		ui.glMeshWidget->vSettings[1].showColorSignature = true;
// 	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayMexicanHatWavelet2()
{
	if (selected[0])
	{
		WaveletMeshProcessor& mp = vMP[0];

		vector<double> vMHW;
		mp.computeMexicanHatWavelet(vMHW, 30, 2);
		mp.normalizeSignatureFrom(vMHW);

		ui.glMeshWidget->vSettings[0].showColorSignature = true;
		ui.actionShowSignature->setChecked(true);
	}

// 	if (selected[1] && vMP[1].mesh)
// 	{
// 		WaveletMeshProcessor& mp = vMP[1];		
// 
// 		vector<double> vMHW;
// 		mp.computeMexicanHatWavelet(vMHW, 30, 2);
// 		mp.normalizeFrom(vMHW);
// 
// 		ui.glMeshWidget->vSettings[1].showColorSignature = true;
// 	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayExperimental()
{
	WaveletMeshProcessor& mp = vMP[0];

 	vector<double> vExp;
 	mp.computeExperimentalWavelet(vExp, 30); 

 	mp.normalizeSignatureFrom(vExp);

 	ui.glMeshWidget->vSettings[0].showColorSignature = true;
	ui.actionShowSignature->setChecked(true);

 	ui.glMeshWidget->update();
 	qout.output("Show MHW from point " + QString::number(mp.pRef));
	
// 	qout.output("Start calculating wavelet of geometry...");
// 	QTime timer;
// 	timer.start();
// 	mp.calGeometryDWT();
// 	qout.output("Finished! Time cost: " + QString::number(timer.elapsed()/1000.0) + " (s)");
}

void QZGeometryWindow::setShowRefPoint()
{
	bool bChecked = ui.actionShowRefPoint->isChecked();
	ui.glMeshWidget->vSettings[0].showRefPoint = bChecked;
	ui.actionShowRefPoint->setChecked(bChecked);
	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayCurvatureMean()
{
	MeshProcessor& mp = vMP[0];

	vector<double> vCurvature;
	mp.computeCurvature(vCurvature, 0);

	auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
	qout.output("Min curvature: " + QString::number(*mm.first) + "  Max curvature: " + QString::number(*mm.second));

	mp.bandCurveSignatureFrom(vCurvature, 0, 1);

	ui.glMeshWidget->vSettings[0].showColorSignature = true;
	ui.actionShowSignature->setChecked(true);

	ui.glMeshWidget->update();
	qout.output("Show Mean Curvature");
}

void QZGeometryWindow::displayCurvatureGauss()
{
	MeshProcessor& mp = vMP[0];

	vector<double> vCurvature;
	mp.computeCurvature(vCurvature, 1);

	auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
	qout.output("Min curvature: " + QString::number(*mm.first) + "  Max curvature: " + QString::number(*mm.second));

	mp.bandCurveSignatureFrom(vCurvature, -1, 1);

	ui.glMeshWidget->vSettings[0].showColorSignature = true;
	ui.actionShowSignature->setChecked(true);

	ui.glMeshWidget->update();
	qout.output("Show Mean Curvature");
}

void QZGeometryWindow::selectObject( int index )
{
	QString text = ui.objSelectBox->itemText(index);
	qout.output("Selected object(s): " + text);
	if (text == "1") 
	{
		selected[0] = ui.glMeshWidget->vSettings[0].selected = true;
		selected[1] = ui.glMeshWidget->vSettings[1].selected = false; 
	}
	else if (text == "2")
	{
		selected[0] = ui.glMeshWidget->vSettings[0].selected = false;
		selected[1] = ui.glMeshWidget->vSettings[1].selected = true; 
	}
	else if (text == "All")
	{
		selected[0] = ui.glMeshWidget->vSettings[0].selected = true;
		selected[1] = ui.glMeshWidget->vSettings[1].selected = true; 
	}
	else if (text == "None")
	{
		selected[0] = ui.glMeshWidget->vSettings[0].selected = false;
		selected[1] = ui.glMeshWidget->vSettings[1].selected = false; 
	}
}

void QZGeometryWindow::clone()
{
	mesh2.cloneFrom(mesh1);
	mesh2.gatherStatistics();

	vMP[1].init(&mesh2, m_ep);
	ui.glMeshWidget->addMesh(&vMP[1]);
	qout.output("Mesh " + QString(mesh2.m_meshName.c_str()) + " constructed! Size: " + QString::number(mesh2.getVerticesNum()));
	
	ui.spinBox2->setMinimum(0);
	ui.spinBox2->setMaximum(mesh2.getVerticesNum()-1);
	ui.horizontalSlider2->setMinimum(0);
	ui.horizontalSlider2->setMaximum(mesh2.getVerticesNum()-1);
	ui.spinBox2->setValue(0);

//	vector<double> vx, vy, vz;
//	vMP[0].reconstructByMHB(300, vx, vy, vz);
//	vMP[0].reconstructByDifferential(vx, vy, vz);
//	vMP[0].reconstructExperimental1(vx, vy, vz);
//	mesh2.setVertexCoordinates(vx, vy, vz);

/*  to prove the effect of scalar product   

	ofstream ofs1("output/dotproduct.dat"), ofs2("output/scalarproduct.dat");
	for (int i = 0; i < vMP[0].mhb.m_nEigFunc; ++i)
	{
		for (int j = 0; j < vMP[0].mhb.m_nEigFunc; ++j)
		{
			const vector<double> &ef1 = vMP[0].mhb.m_func[i].m_vec, 
				                 &ef2 = vMP[0].mhb.m_func[j].m_vec;
			const vector<double> &sc = vMP[0].mLaplacian.getVerticesWeight();

			ofs1 << VectorDotProduct(ef1, ef2) << ' ';
			ofs2 << VectorScalarProduct(ef1, ef2, sc) << ' ';
		}
		ofs1 << endl;
		ofs2 << endl;
		if (i % 10 == 0)
			qout.output("Row " + QString::number(i) + " finished!"); 
	}
	ofs1.close();
	ofs2.close();
*/

	ui.glMeshWidget->update();
}

void QZGeometryWindow::keyPressEvent( QKeyEvent *event )
{
	switch (event->key())
	{
	case Qt::Key_1:
		if (event->modifiers() & Qt::ControlModifier)
		{
			ui.objSelectBox->setCurrentIndex(ui.objSelectBox->findText("1"));
			selectObject(ui.objSelectBox->findText("1"));
		}
		break;

	case Qt::Key_2:
		if (event->modifiers() & Qt::ControlModifier)
		{
			ui.objSelectBox->setCurrentIndex(ui.objSelectBox->findText("2"));
			selectObject(ui.objSelectBox->findText("2"));
		}
		break;

	case Qt::Key_0:
		if (event->modifiers() & Qt::ControlModifier)
		{
			ui.objSelectBox->setCurrentIndex(ui.objSelectBox->findText("All"));
			selectObject(ui.objSelectBox->findText("All"));
		}
		break;
	
	case Qt::Key_X:
		if (event->modifiers() & Qt::ShiftModifier)
			refMove.xMove--;
		else refMove.xMove++;
		updateReferenceMove();
		break;
	
	case Qt::Key_Y:
		if (event->modifiers() & Qt::ShiftModifier)
			refMove.yMove--;
		else refMove.yMove++;
		updateReferenceMove();
		break;
	
	case Qt::Key_Z:
		if (event->modifiers() & Qt::ShiftModifier)
			refMove.zMove--;
		else refMove.zMove++;
		updateReferenceMove();
		break;

	default: QWidget::keyPressEvent(event);
	}
}

void QZGeometryWindow::updateReferenceMove()
{
	double unitMove = (mesh1.getBoundingBox().x + mesh1.getBoundingBox().y + mesh1.getBoundingBox().z)/300.0;
	Vector3D originalPos = mesh1.getVertex(vMP[0].pRef)->getPos();
	vMP[0].posRef.x = originalPos.x + unitMove * refMove.xMove;
	vMP[0].posRef.y = originalPos.y + unitMove * refMove.yMove;
	vMP[0].posRef.z = originalPos.z + unitMove * refMove.zMove;

	ui.glMeshWidget->update();
}

void QZGeometryWindow::deformExperimental()
{
	vector<double> vx, vy, vz;
//	vMP[0].reconstructByMHB(300, vx, vy, vz);
//	vMP[0].reconstructByDifferential(vx, vy, vz, true);
	vMP[0].reconstructBySGW(vx, vy, vz, true);
//	vMP[0].reconstructExperimental1(vx, vy, vz, true);
	mesh2.setVertexCoordinates(vx, vy, vz);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayPointCloud()
{
	ui.actionDisplayPointCloud->setChecked(true);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(false);

	ui.glMeshWidget->vSettings[0].displayType = ui.glMeshWidget->vSettings[1].displayType = DisplaySettings::PointCloud;
	ui.glMeshWidget->vSettings[0].glPolygonMode = ui.glMeshWidget->vSettings[1].glPolygonMode = GL_POINT;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayWireframe()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(true);
	ui.actionDisplayMesh->setChecked(false);

	ui.glMeshWidget->vSettings[0].displayType = ui.glMeshWidget->vSettings[1].displayType = DisplaySettings::Wireframe;
	ui.glMeshWidget->vSettings[0].glPolygonMode = ui.glMeshWidget->vSettings[1].glPolygonMode = GL_LINE;
	ui.glMeshWidget->update();

}

void QZGeometryWindow::displayMesh()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(true);

	ui.glMeshWidget->vSettings[0].displayType = ui.glMeshWidget->vSettings[1].displayType = DisplaySettings::Mesh;
	ui.glMeshWidget->vSettings[0].glPolygonMode = ui.glMeshWidget->vSettings[1].glPolygonMode = GL_FILL;
	ui.glMeshWidget->update();

}

void QZGeometryWindow::setShowSignature()
{
	bool bToShow = ui.actionShowSignature->isChecked();
	ui.actionShowSignature->setChecked(bToShow);
	ui.glMeshWidget->vSettings[0].showColorSignature = bToShow;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeSGWSFeatures()
{
	vMP[0].vFeatures.clear();
	double timescales[4] = {5, 10, 20, 40};
	for (int s = 0; s < 4; ++s)
	{
		vector<double> vSig;
		vector<int> vFeatures;
		vMP[0].getSGWSignature(timescales[s], vSig);
		mesh1.extractExtrema(vSig, 2, 1e-5, vFeatures);
		for (vector<int>::iterator iter = vFeatures.begin(); iter != vFeatures.end(); ++iter)
		{
			vMP[0].vFeatures.push_back(ManifoldFeature(*iter, s));
		}
	}

	ui.glMeshWidget->vSettings[0].showFeatures = true;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::filterExperimental()
{
	vector<double> vx, vy, vz;
	vMP[0].filterBySGW(vx, vy, vz);
	mesh2.setVertexCoordinates(vx, vy, vz);

	double errorSum(0);
	for (int i = 0; i < mesh1.getVerticesNum(); ++i)
	{
		errorSum += (mesh1.getVertex_const(i)->getPos() - mesh2.getVertex_const(i)->getPos()).length();
	}
	errorSum /= mesh1.getVerticesNum() * mesh1.getAvgEdgeLength();
	qout.output("Average position error: " + QString::number(errorSum));

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayDiffPosition()
{
	assert(mesh1.getVerticesNum() == mesh2.getVerticesNum());
	int size = mesh1.getVerticesNum();
	vector<double> vDiff;
	vDiff.resize(size);
	
	for (int i = 0; i < mesh1.getVerticesNum(); ++i)
	{
		vDiff[i] = (mesh1.getVertex_const(i)->getPos() - mesh2.getVertex_const(i)->getPos()).length() / mesh1.getAvgEdgeLength();
	}

	vMP[0].normalizeSignatureFrom(vDiff);
	ui.glMeshWidget->vSettings[0].showColorSignature = true;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setShowColorLegend()
{
	bool bChecked = ui.actionColorLegend->isChecked();
	ui.glMeshWidget->m_bShowLegend = bChecked;
	ui.actionColorLegend->setChecked(bChecked);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setShowFeatures()
{
	// TODO
}
