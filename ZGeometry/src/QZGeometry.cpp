#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <QtGui/QMessageBox>
#include <QTime>
#include "QZGeometry.h"
#include <ZUtil.h>
#include <set>

#define DEFAULT_EIGEN_SIZE 500
#define DEFAULT_DEFORM_RING 5 
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
	
	m_commonParameter = 50;
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
}

QZGeometryWindow::~QZGeometryWindow()
{
	engClose(m_ep);
}

void QZGeometryWindow::makeConnections()
{	
	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	
	////////	compute	////////
	QObject::connect(ui.actionComputeLaplacian, SIGNAL(triggered()), this, SLOT(computeLaplacian()));
	QObject::connect(ui.actionComputeSGW, SIGNAL(triggered()), this, SLOT(computeSGW()));
	QObject::connect(ui.actionSGWSFeatures, SIGNAL(triggered()), this, SLOT(computeSGWSFeatures()));

	////////    Control	////////
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), ui.horizontalSlider1, SLOT(setValue(int)));
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), this, SLOT(setRefPoint1(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), ui.spinBox1, SLOT(setValue(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), this, SLOT(setRefPoint1(int)));
	QObject::connect(ui.spinBoxParameter, SIGNAL(valueChanged(int)), ui.horizontalSliderParamter, SLOT(setValue(int)));
	QObject::connect(ui.spinBoxParameter, SIGNAL(valueChanged(int)), this, SLOT(setCommonParameter(int)));
	QObject::connect(ui.horizontalSliderParamter, SIGNAL(valueChanged(int)), ui.spinBoxParameter, SLOT(setValue(int)));
	QObject::connect(ui.horizontalSliderParamter, SIGNAL(valueChanged(int)), this, SLOT(setCommonParameter(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked(int)), this, SLOT(setRefPoint1(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked(int)), ui.horizontalSlider1, SLOT(setValue(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked(int)), ui.spinBox1, SLOT(setValue(int)));
	QObject::connect(ui.boxObjSelect, SIGNAL(activated(int)), this, SLOT(selectObject(int)));
	QObject::connect(ui.actionEditMove, SIGNAL(triggered()), this, SLOT(setEditModeMove()));
	QObject::connect(ui.actionEditPick, SIGNAL(triggered()), this, SLOT(setEditModePick()));
	QObject::connect(ui.actionEditDrag, SIGNAL(triggered()), this, SLOT(setEditModeDrag()));

	////////	Edit	////////
	QObject::connect(ui.actionDeformSimple, SIGNAL(triggered()), this, SLOT(deformSimple()));
	QObject::connect(ui.actionDeformSGW, SIGNAL(triggered()), this, SLOT(deformSGW()));
	QObject::connect(ui.actionDeformLaplace, SIGNAL(triggered()), this, SLOT(deformLaplace()));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
	QObject::connect(ui.actionReconstructMHB, SIGNAL(triggered()), this, SLOT(reconstructByMHB()));
	QObject::connect(ui.actionReconstructSGW, SIGNAL(triggered()), this, SLOT(reconstructBySGW()));
	QObject::connect(ui.actionFilter_1, SIGNAL(triggered()), this, SLOT(filterExperimental()));
	 
	////////	Display	////////
	QObject::connect(ui.actionDisplayMesh, SIGNAL(triggered()), this, SLOT(setDisplayMesh()));
	QObject::connect(ui.actionDisplayWireframe, SIGNAL(triggered()), this, SLOT(setDisplayWireframe()));
	QObject::connect(ui.actionDisplayPointCloud, SIGNAL(triggered()), this, SLOT(setDisplayPointCloud()));
	QObject::connect(ui.actionDisplayNeighbors, SIGNAL(triggered()), this, SLOT(displayNeighborVertices()));
	QObject::connect(ui.actionShowFeatures, SIGNAL(triggered(bool)), this, SLOT(toggleShowFeatures(bool)));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered(bool)), this, SLOT(toggleShowRefPoint(bool)));
	QObject::connect(ui.actionShowSignature, SIGNAL(triggered(bool)), this, SLOT(toggleShowSignature(bool)));
	QObject::connect(ui.actionColorLegend, SIGNAL(triggered(bool)), this, SLOT(toggleShowColorLegend(bool)));
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(displayEigenfunction()));
	QObject::connect(ui.actionMexicanHatWavelet1, SIGNAL(triggered()), this, SLOT(displayMexicanHatWavelet1()));
	QObject::connect(ui.actionMexicanHatWavelet2, SIGNAL(triggered()), this, SLOT(displayMexicanHatWavelet2()));
	QObject::connect(ui.actionExperimental, SIGNAL(triggered()), this, SLOT(displayExperimental()));
	QObject::connect(ui.actionMeanCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureMean()));
	QObject::connect(ui.actionGaussCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureGauss()));
	QObject::connect(ui.actionDiffPosition, SIGNAL(triggered()), this, SLOT(displayDiffPosition()));
}

bool QZGeometryWindow::initialize()
{
	qout.output("******** Welcome ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);

	if (!(m_ep = engOpen("\0")))
	{
		qout.output("Can't start MATLAB engine!", OUT_MSGBOX);
		return false;
	}
	else qout.output("Matlab engine initialized!", OUT_CONSOLE);
	qout.output('*', 24);
	
	setDisplayMesh();
	setEditModeMove();

	/// load meshes    
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
	
	Vector3D center1 = mesh1.m_Center, bbox1 = mesh1.m_bBox;
	qout.output(qformat.sprintf("Center: (%f,%f,%f)\nDimension: (%f,%f,%f)", center1.x, center1.y, center1.z, bbox1.x, bbox1.y, bbox1.z));

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
	
	/// update ui
	vMP[0].init(&mesh1, m_ep);
	ui.glMeshWidget->addMesh(&vMP[0]);
	ui.glMeshWidget->fieldView(center1, bbox1);

	ui.spinBox1->setMinimum(0);
	ui.spinBox1->setMaximum(mesh1.getVerticesNum()-1);
	ui.horizontalSlider1->setMinimum(0);
	ui.horizontalSlider1->setMaximum(mesh1.getVerticesNum()-1);
	ui.spinBox1->setValue(0);	
	ui.spinBoxParameter->setMinimum(0);
	ui.spinBoxParameter->setMaximum(100);
	ui.horizontalSliderParamter->setMinimum(0);
	ui.horizontalSliderParamter->setMaximum(100);
	ui.spinBoxParameter->setValue(50);

	selected[0] = ui.glMeshWidget->vSettings[0].selected = true;
	selected[1] = ui.glMeshWidget->vSettings[1].selected = false; 

	this->computeLaplacian(0);	
	qout.output("Non-zeros of Laplacian: " + Int2String(vMP[0].mLaplacian.getNonzeroNum()));
	vMP[0].mLaplacian.dumpLaplacian("output/sparse_laplacian.csv");

	return true;
}

void QZGeometryWindow::keyPressEvent( QKeyEvent *event )
{
	switch (event->key())
	{
	case Qt::Key_1:
		if (event->modifiers() & Qt::ControlModifier)
		{
			ui.boxObjSelect->setCurrentIndex(ui.boxObjSelect->findText("1"));
			selectObject(ui.boxObjSelect->findText("1"));
		}
		break;

	case Qt::Key_2:
		if (event->modifiers() & Qt::ControlModifier)
		{
			ui.boxObjSelect->setCurrentIndex(ui.boxObjSelect->findText("2"));
			selectObject(ui.boxObjSelect->findText("2"));
		}
		break;

	case Qt::Key_0:
		if (event->modifiers() & Qt::ControlModifier)
		{
			ui.boxObjSelect->setCurrentIndex(ui.boxObjSelect->findText("All"));
			selectObject(ui.boxObjSelect->findText("All"));
		}
		break;

	case Qt::Key_C:
		clone();
		break;

	case Qt::Key_F:
		toggleShowFeatures();
		break;

	case Qt::Key_R:
		toggleShowRefPoint();
		break;

	case Qt::Key_S:
		toggleShowSignature();
		break;

	case Qt::Key_W:	// switch between display mode
	{
		if (ui.glMeshWidget->vSettings[0].displayType == DisplaySettings::Mesh)
			setDisplayWireframe();
		else if (ui.glMeshWidget->vSettings[0].displayType == DisplaySettings::Wireframe)
			setDisplayPointCloud();
		else setDisplayMesh();
		break;
	}	
		
	case Qt::Key_M:
		setEditModeMove();
		break;

	case Qt::Key_P:
		setEditModePick();
		break;

	case Qt::Key_D:
		setEditModeDrag();
		break;

	case Qt::Key_E:
		deformSimple();
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

void QZGeometryWindow::computeLaplacian( int obj /*= 0*/ )
{
	QTime timer;
	timer.start();

	string pathMHB = "output/" + vMP[obj].getMesh()->m_meshName + ".mhb";
	
	ifstream ifs(pathMHB.c_str());
	if (ifs && LOAD_MHB_CACHE)	// MHB cache available for the current mesh
	{
		vMP[obj].readMHB(pathMHB);
	}
	else 
	{
		int nEig = DEFAULT_EIGEN_SIZE;
		if (nEig > mesh1.getVerticesNum()-1)
			nEig = mesh1.getVerticesNum() - 1;
		vMP[obj].decomposeLaplacian(nEig);
//		vMP[obj].decomposeLaplacian(vMP[obj].m_size);
		vMP[obj].writeMHB(pathMHB);
		qout.output("MHB saved to " + pathMHB);

		std::string pathEVL = "output/" + vMP[obj].getMesh()->m_meshName + ".evl";	//eigenvalues
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
	qout.output(qformat.sprintf("Min Eig Val: %f, Max Eig Val: %f", vMP[obj].mhb.m_func.front().m_val, vMP[obj].mhb.m_func.back().m_val));
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
			vMP[0].vFeatures.push_back(MeshFeature(*iter, s));
		}
	}

	if (!ui.actionShowFeatures->isChecked())
		toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeSGW()
{
	CStopWatch timer;
	timer.startTimer();

	double timescales[4] = {5, 10, 20, 40};
	vector<double> vTimes;
	std::copy(timescales, timescales + 4, vTimes.begin());

	vMP[0].computeSGW(vTimes);

	timer.stopTimer();
	qout.output(QString("Time for compute SGW: ") + QString::number(timer.getElapsedTime()));
}

void QZGeometryWindow::selectObject( int index )
{
	QString text = ui.boxObjSelect->itemText(index);
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

void QZGeometryWindow::setRefPoint1( int vn )
{
	vMP[0].pRef = vn;
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	updateReferenceMove();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setCommonParameter( int p )
{
	m_commonParameter = p;
}

void QZGeometryWindow::setEditModeMove()
{
	ui.glMeshWidget->editMode = GLMeshWidget::QZ_MOVE;
	ui.actionEditMove->setChecked(true);
	ui.actionEditDrag->setChecked(false);
	ui.actionEditPick->setCheckable(false);

	qout.output("Edit Mode: Move", OUT_STATUS);
}

void QZGeometryWindow::setEditModePick()
{
	ui.glMeshWidget->editMode = GLMeshWidget::QZ_PICK;
	ui.actionEditMove->setChecked(false);
	ui.actionEditDrag->setChecked(false);
	ui.actionEditPick->setCheckable(true);

	qout.output("Edit Mode: Pick", OUT_STATUS);
}

void QZGeometryWindow::setEditModeDrag()
{
	ui.glMeshWidget->editMode = GLMeshWidget::QZ_DRAG;
	ui.actionEditMove->setChecked(false);
	ui.actionEditDrag->setChecked(true);
	ui.actionEditPick->setCheckable(false);

	qout.output("Edit Mode: Drag", OUT_STATUS);
}

void QZGeometryWindow::toggleShowRefPoint(bool x)
{
	bool bChecked = !ui.glMeshWidget->m_bShowRefPoint;
	ui.glMeshWidget->m_bShowRefPoint = bChecked;
	ui.actionShowRefPoint->setChecked(bChecked);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowColorLegend( bool show /*= false*/ )
{
	bool bChecked = !ui.glMeshWidget->m_bShowLegend;
	ui.glMeshWidget->m_bShowLegend = bChecked;
	ui.actionColorLegend->setChecked(bChecked);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowFeatures( bool show /*= false*/ )
{
	bool bChecked = !ui.glMeshWidget->m_bShowFeatures;
	ui.glMeshWidget->m_bShowFeatures = bChecked;
	ui.actionShowFeatures->setChecked(bChecked);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowSignature( bool show /*= false*/ )
{
	bool bToShow = !ui.glMeshWidget->m_bShowSignature;
	ui.glMeshWidget->m_bShowSignature = bToShow;
	ui.actionShowSignature->setChecked(bToShow);
	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayPointCloud()
{
	ui.actionDisplayPointCloud->setChecked(true);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(false);

	ui.glMeshWidget->vSettings[0].displayType = ui.glMeshWidget->vSettings[1].displayType = DisplaySettings::PointCloud;
	ui.glMeshWidget->vSettings[0].glPolygonMode = ui.glMeshWidget->vSettings[1].glPolygonMode = GL_POINT;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayWireframe()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(true);
	ui.actionDisplayMesh->setChecked(false);

	ui.glMeshWidget->vSettings[0].displayType = ui.glMeshWidget->vSettings[1].displayType = DisplaySettings::Wireframe;
	ui.glMeshWidget->vSettings[0].glPolygonMode = ui.glMeshWidget->vSettings[1].glPolygonMode = GL_LINE;
	ui.glMeshWidget->update();

}

void QZGeometryWindow::setDisplayMesh()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(true);

	ui.glMeshWidget->vSettings[0].displayType = ui.glMeshWidget->vSettings[1].displayType = DisplaySettings::Mesh;
	ui.glMeshWidget->vSettings[0].glPolygonMode = ui.glMeshWidget->vSettings[1].glPolygonMode = GL_FILL;
	ui.glMeshWidget->update();

}

void QZGeometryWindow::displayEigenfunction()
{
	int select_eig = (m_commonParameter - 49 > 1) ? (m_commonParameter - 49) : 1;

	if (selected[0])
	{
		DifferentialMeshProcessor& mp = vMP[0];
		mp.normalizeSignatureFrom(mp.mhb.m_func[select_eig].m_vec);
	}

// 	if (selected[1] && vMP[1].mesh)
// 	{
// 		MeshProcessor& mp = vMP[1];		
// 		mp.normalizeFrom(mp.mhb.m_func[1].m_vec);
// 		ui.glMeshWidget->vSettings[1].showColorSignature = true;
// 	}

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output("Show eigenfunction" + Int2String(select_eig));
}

void QZGeometryWindow::displayMexicanHatWavelet1()
{
	if (selected[0])
	{
		DifferentialMeshProcessor& mp = vMP[0];

		vector<double> vMHW;
		mp.computeMexicanHatWavelet(vMHW, 30, 1);
		mp.normalizeSignatureFrom(vMHW);

	}

// 	if (selected[1] && vMP[1].mesh)
// 	{
// 		DifferentialMeshProcessor& mp = vMP[1];		
// 
// 		vector<double> vMHW;
// 		mp.computeMexicanHatWavelet(vMHW, 30, 1);
// 		mp.normalizeFrom(vMHW);
// 
// 		ui.glMeshWidget->vSettings[1].showColorSignature = true;
// 	}

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayMexicanHatWavelet2()
{
	if (selected[0])
	{
		DifferentialMeshProcessor& mp = vMP[0];

		vector<double> vMHW;
		mp.computeMexicanHatWavelet(vMHW, 30, 2);
		mp.normalizeSignatureFrom(vMHW);
	}

// 	if (selected[1] && vMP[1].mesh)
// 	{
// 		DifferentialMeshProcessor& mp = vMP[1];		
// 
// 		vector<double> vMHW;
// 		mp.computeMexicanHatWavelet(vMHW, 30, 2);
// 		mp.normalizeFrom(vMHW);
// 
// 		ui.glMeshWidget->vSettings[1].showColorSignature = true;
// 	}

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayExperimental()
{
	DifferentialMeshProcessor& mp = vMP[0];

 	vector<double> vExp;
 	mp.computeExperimentalWavelet(vExp, 30); 

 	mp.normalizeSignatureFrom(vExp);

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output("Show MHW from vertex #" + QString::number(mp.pRef));
	
// 	qout.output("Start calculating wavelet of geometry...");
// 	QTime timer;
// 	timer.start();
// 	mp.calGeometryDWT();
// 	qout.output("Finished! Time cost: " + QString::number(timer.elapsed()/1000.0) + " (s)");
}

void QZGeometryWindow::displayCurvatureMean()
{
	DifferentialMeshProcessor& mp = vMP[0];

	vector<double> vCurvature;
	mp.computeCurvature(vCurvature, 0);

	auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
	qout.output("Min curvature: " + QString::number(*mm.first) + "  Max curvature: " + QString::number(*mm.second));

	mp.bandCurveSignatureFrom(vCurvature, 0, 1);

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output("Show Mean Curvature");
}

void QZGeometryWindow::displayCurvatureGauss()
{
	DifferentialMeshProcessor& mp = vMP[0];

	vector<double> vCurvature;
	mp.computeCurvature(vCurvature, 1);

	auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
	qout.output("Min curvature: " + QString::number(*mm.first) + "  Max curvature: " + QString::number(*mm.second));

	mp.bandCurveSignatureFrom(vCurvature, -1, 1);

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output("Show Mean Curvature");
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
	
	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();
	
	ui.glMeshWidget->update();
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

void QZGeometryWindow::clone()
{
	mesh2.cloneFrom(mesh1);
	mesh2.gatherStatistics();

	vMP[1].init(&mesh2, m_ep);
	ui.glMeshWidget->addMesh(&vMP[1]);
	qout.output("Mesh " + QString(mesh2.m_meshName.c_str()) + " constructed! Size: " + QString::number(mesh2.getVerticesNum()));

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

void QZGeometryWindow::reconstructByMHB()
{
//	assert(vMP[1].mesh);
	
	double ratio = min((double)m_commonParameter/50.0, 1.0);
	int nEig = vMP[0].mhb.m_nEigFunc * ratio;

	vector<double> vx, vy, vz;
	vMP[0].reconstructByMHB(nEig, vx, vy, vz);
	mesh2.setVertexCoordinates(vx, vy, vz);
	ui.glMeshWidget->update();

	double errorSum(0);
	for (int i = 0; i < mesh1.getVerticesNum(); ++i)
	{
		errorSum += (mesh1.getVertex_const(i)->getPos() - mesh2.getVertex_const(i)->getPos()).length();
	}
	errorSum /= mesh1.getVerticesNum() * mesh1.getAvgEdgeLength();
	qout.output("MHB reconstruction with " + Int2String(nEig) + " MHBs.");
	qout.output("Average position error: " + QString::number(errorSum));
}

void QZGeometryWindow::reconstructBySGW()
{
//	assert(vMP[1].mesh);

	vector<double> vx, vy, vz;
	vMP[0].reconstructBySGW(vx, vy, vz, true);
	mesh2.setVertexCoordinates(vx, vy, vz);

	double errorSum(0);
	for (int i = 0; i < mesh1.getVerticesNum(); ++i)
	{
		errorSum += (mesh1.getVertex_const(i)->getPos() - mesh2.getVertex_const(i)->getPos()).length();
	}
	errorSum /= mesh1.getVerticesNum() * mesh1.getAvgEdgeLength();
	qout.output("SGW reconstruction.");
	qout.output("Average position error: " + QString::number(errorSum));
}

void QZGeometryWindow::deformSimple()
{
	//	vector<double> vx, vy, vz;
	//	vMP[0].reconstructByMHB(300, vx, vy, vz);
	//	vMP[0].reconstructByDifferential(vx, vy, vz, true);
	//	vMP[0].reconstructBySGW(vx, vy, vz, true);
	//	vMP[0].reconstructExperimental1(vx, vy, vz, true);
	//	mesh2.setVertexCoordinates(vx, vy, vz);

	int activeHandle = vMP[0].active_handle; 
	vector<int> vHandle;
	vHandle.push_back(activeHandle);
	vector<Vector3D> vHandlePos;
	vHandlePos.push_back(vMP[0].mHandles[activeHandle]);
	vector<int> vFree = vMP[0].getMesh()->getNeighboringVertex(activeHandle, DEFAULT_DEFORM_RING);
	vector<Vector3D> vNewPos;

	vMP[0].deform(vHandle, vHandlePos, vFree, vNewPos, Simple);
	mesh2.setVertexCoordinates(vFree, vNewPos);
	mesh2.setVertexCoordinates(vHandle, vHandlePos);

	setEditModeMove();
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

void QZGeometryWindow::displayNeighborVertices()
{
	int ring = (m_commonParameter > 45) ? (m_commonParameter-45) : 1;

	int ref = vMP[0].pRef;
	std::vector<int> vn = vMP[0].getMesh()->getNeighboringVertex(ref, ring);
//	std::vector<int> vn = vMP[0].getMesh()->getRingVertex(ref, ring);
	MeshFeatureList *mfl = new MeshFeatureList;
	for (auto iter = vn.begin(); iter != vn.end(); ++iter)
	{
		mfl->m_vFeatures.push_back(MeshFeature(*iter));
		mfl->setIDandName(1, "Neighbors");
	}
	vMP[0].addProperty(mfl);

	vMP[0].vFeatures = mfl->m_vFeatures;
	
	if (!ui.actionShowFeatures->isChecked())
		toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::deformSGW()
{
	int activeHandle = vMP[0].active_handle; 
	vector<int> vHandle;
	vHandle.push_back(activeHandle);
	vector<Vector3D> vHandlePos;
	vHandlePos.push_back(vMP[0].mHandles[activeHandle]);
	vector<int> vFree = vMP[0].getMesh()->getNeighboringVertex(activeHandle, DEFAULT_DEFORM_RING);
	vector<Vector3D> vNewPos;

	vMP[0].deform(vHandle, vHandlePos, vFree, vNewPos, SGW);
	mesh2.setVertexCoordinates(vFree, vNewPos);
	mesh2.setVertexCoordinates(vHandle, vHandlePos);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::deformLaplace()
{
	vector<int> vHandle;
	vector<Vector3D> vHandlePos;
	vector<int> vFree;
// 	int activeHandle = vMP[0].active_handle; 	
// 	vHandle.push_back(activeHandle);	
// 	vHandlePos.push_back(vMP[0].mHandles[activeHandle]);
//	vFree = vMP[0].getMesh()->getNeighboringVertex(activeHandle, DEFAULT_DEFORM_RING);
	
	std::set<int> sFreeIdx;
	for (auto iter = vMP[0].mHandles.begin(); iter!= vMP[0].mHandles.end(); ++iter)
	{
		vHandle.push_back(iter->first);
		vHandlePos.push_back(iter->second);

		vector<int> vNeighbor = vMP[0].getMesh()->getNeighboringVertex(iter->first, DEFAULT_DEFORM_RING);
		sFreeIdx.insert(vNeighbor.begin(), vNeighbor.end());
	}
	
	vFree.insert(vFree.begin(), sFreeIdx.begin(), sFreeIdx.end());
	vector<Vector3D> vNewPos;

	try
	{
		vMP[0].deform(vHandle, vHandlePos, vFree, vNewPos, Laplace);
		mesh2.setVertexCoordinates(vFree, vNewPos);
		mesh2.setVertexCoordinates(vHandle, vHandlePos);
	}
	catch (runtime_error* e)
	{
		qout.output(e->what(), OUT_MSGBOX);
	}

	ui.glMeshWidget->update();
}
