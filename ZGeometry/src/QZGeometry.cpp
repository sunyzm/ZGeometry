#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <set>
#include <exception>
#include <stdexcept>
#include <QtGui/QMessageBox>
#include <QTime>
#include <QProcess>
#include <ZUtil.h>
#include <SimpleConfigLoader.h>
#include <ppl.h>
#include "QZGeometry.h"
#include "common.h"

using namespace std;

extern OutputHelper qout;
extern SimpleConfigLoader g_configMgr;
extern GeometryTask g_task;

int QZGeometryWindow::DEFAULT_EIGEN_SIZE = 300;
int QZGeometryWindow::DEFAULT_DEFORM_RING = 5 ;
int QZGeometryWindow::LOAD_MHB_CACHE = 1;
double QZGeometryWindow::MIN_HK_TIMESCALE = 1e-2;
double QZGeometryWindow::DEFUALT_HK_TIMESCALE = 40.0;
double QZGeometryWindow::MAX_HK_TIMESCALE = 2000.0;
double QZGeometryWindow::PARAMETER_SLIDER_CENTER = 50;
double QZGeometryWindow::DR_THRESH_INCREMENT  = 0.00001;
double QZGeometryWindow::MATCHING_THRESHOLD = 0.002;

QZGeometryWindow::QZGeometryWindow(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	m_ep = NULL;
	num_meshes = 1;

	//* read in configuration parameters from g_configMgr *//
	LOAD_MHB_CACHE = g_configMgr.getConfigValueInt("LOAD_MHB_CACHE");
	PARAMETER_SLIDER_CENTER = g_configMgr.getConfigValueInt("PARAMETER_SLIDER_CENTER");
	DEFUALT_HK_TIMESCALE = g_configMgr.getConfigValueDouble("DEFUALT_HK_TIMESCALE");

	mesh_valid[0] = mesh_valid[1] = false;
	m_commonParameter = PARAMETER_SLIDER_CENTER;
	current_operation = None;
	deformType = Simple;
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	
	ui.setupUi(this);
	ui.centralWidget->setLayout(ui.mainLayout);
	ui.horizontalSliderParamter->setMinimum(0);
	ui.horizontalSliderParamter->setMaximum(2*PARAMETER_SLIDER_CENTER-1);
	ui.horizontalSliderParamter->setSliderPosition(PARAMETER_SLIDER_CENTER);

	m_actionLaplacians.resize(LaplacianEnd);
	m_actionLaplacians[Umbrella] = ui.actionComputeLaplacian1;
	m_actionLaplacians[CotFormula] = ui.actionComputeLaplacian2;
	m_actionLaplacians[Anisotropic1] = ui.actionComputeLaplacian3;
	m_actionLaplacians[Anisotropic2] = ui.actionComputeLaplacian4;
	m_actionLaplacians[IsoApproximate] = ui.actionComputeLaplacian5;


	this->makeConnections();
	
	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);
}

QZGeometryWindow::~QZGeometryWindow()
{
	engClose(m_ep);
}

void QZGeometryWindow::makeConnections()
{	
	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	
	////////	compute	////////
	QObject::connect(ui.actionComputeLaplacian1, SIGNAL(triggered()), this, SLOT(computeLaplacian1()));
	QObject::connect(ui.actionComputeLaplacian2, SIGNAL(triggered()), this, SLOT(computeLaplacian2()));
	QObject::connect(ui.actionComputeLaplacian3, SIGNAL(triggered()), this, SLOT(computeLaplacian3()));
	QObject::connect(ui.actionComputeLaplacian4, SIGNAL(triggered()), this, SLOT(computeLaplacian4()));
	QObject::connect(ui.actionComputeHK, SIGNAL(triggered()), this, SLOT(computeHK()));
	QObject::connect(ui.actionComputeHKS, SIGNAL(triggered()), this, SLOT(computeHKS()));
	QObject::connect(ui.actionComputeHKSFeatures, SIGNAL(triggered()), this, SLOT(computeHKSFeatures()));
	QObject::connect(ui.actionComputeMHW, SIGNAL(triggered()), this, SLOT(computeMHW()));
	QObject::connect(ui.actionComputeMHWS, SIGNAL(triggered()), this, SLOT(computeMHWS()));
	QObject::connect(ui.actionComputeMHWSFeatures, SIGNAL(triggered()), this, SLOT(computeMHWFeatures()));
	QObject::connect(ui.actionComputeSGWS, SIGNAL(triggered()), this, SLOT(computeSGWS()));
	QObject::connect(ui.actionComputeSGW, SIGNAL(triggered()), this, SLOT(computeSGW()));
	QObject::connect(ui.actionComputeSGWSFeatures, SIGNAL(triggered()), this, SLOT(computeSGWSFeatures()));
	QObject::connect(ui.actionComputeBiharmonic, SIGNAL(triggered()), this, SLOT(computeBiharmonic()));
	QObject::connect(ui.actionComputeSimilarityMap, SIGNAL(triggered()), this, SLOT(computeSimilarityMap()));
	QObject::connect(ui.actionComputeSimilarityMap2, SIGNAL(triggered()), this, SLOT(computeSimilarityMap2()));

	////////    Control	////////
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked1(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked1(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked1(int)), this, SLOT(setRefPoint1(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked1(int)), ui.horizontalSlider1, SLOT(setValue(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked1(int)), ui.spinBox1, SLOT(setValue(int)));

	QObject::connect(ui.spinBox2, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked2(int)));
	QObject::connect(ui.horizontalSlider2, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked2(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked2(int)), this, SLOT(setRefPoint2(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked2(int)), ui.horizontalSlider2, SLOT(setValue(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked2(int)), ui.spinBox2, SLOT(setValue(int)));

	QObject::connect(ui.spinBoxParameter, SIGNAL(valueChanged(int)), ui.horizontalSliderParamter, SLOT(setValue(int)));
	QObject::connect(ui.spinBoxParameter, SIGNAL(valueChanged(int)), this, SLOT(setCommonParameter(int)));
	QObject::connect(ui.horizontalSliderParamter, SIGNAL(valueChanged(int)), ui.spinBoxParameter, SLOT(setValue(int)));
	QObject::connect(ui.horizontalSliderParamter, SIGNAL(valueChanged(int)), this, SLOT(setCommonParameter(int)));
	QObject::connect(ui.boxObjSelect, SIGNAL(activated(int)), this, SLOT(selectObject(int)));
	QObject::connect(ui.actionEditMove, SIGNAL(triggered()), this, SLOT(setEditModeMove()));
	QObject::connect(ui.actionEditPick, SIGNAL(triggered()), this, SLOT(setEditModePick()));
	QObject::connect(ui.actionEditDrag, SIGNAL(triggered()), this, SLOT(setEditModeDrag()));

	////////	Edit	////////
	QObject::connect(ui.actionDeformSimple, SIGNAL(triggered()), this, SLOT(deformSimple()));
	QObject::connect(ui.actionDeformSGW, SIGNAL(triggered()), this, SLOT(deformSGW()));
	QObject::connect(ui.actionDeformLaplace, SIGNAL(triggered()), this, SLOT(deformLaplace()));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
	QObject::connect(ui.actionReconstructMHB, SIGNAL(triggered()), this, SLOT(reconstructMHB()));
	QObject::connect(ui.actionReconstructSGW, SIGNAL(triggered()), this, SLOT(reconstructSGW()));
	QObject::connect(ui.actionFilter_1, SIGNAL(triggered()), this, SLOT(filterExperimental()));
	 
	////////	Display	////////
	QObject::connect(ui.actionDisplayMesh, SIGNAL(triggered()), this, SLOT(setDisplayMesh()));
	QObject::connect(ui.actionDisplayWireframe, SIGNAL(triggered()), this, SLOT(setDisplayWireframe()));
	QObject::connect(ui.actionDisplayPointCloud, SIGNAL(triggered()), this, SLOT(setDisplayPointCloud()));
	QObject::connect(ui.actionDisplayNeighbors, SIGNAL(triggered()), this, SLOT(displayNeighborVertices()));
	QObject::connect(ui.actionShowFeatures, SIGNAL(triggered(bool)), this, SLOT(toggleShowFeatures(bool)));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered(bool)), this, SLOT(toggleShowRefPoint(bool)));
	QObject::connect(ui.actionShowSignature, SIGNAL(triggered(bool)), this, SLOT(toggleShowSignature(bool)));
	QObject::connect(ui.actionShowColorLegend, SIGNAL(triggered(bool)), this, SLOT(toggleShowColorLegend(bool)));
	QObject::connect(ui.actionDrawMatching, SIGNAL(triggered(bool)), this, SLOT(toggleDrawMatching(bool)));
	QObject::connect(ui.actionShowMatchingLines, SIGNAL(triggered(bool)), this, SLOT(toggleShowMatchingLines(bool)));
	QObject::connect(ui.actionDrawRegistration, SIGNAL(triggered(bool)), this, SLOT(toggleDrawRegistration(bool)));
	QObject::connect(ui.actionDisplayEigenfunction, SIGNAL(triggered()), this, SLOT(displayEigenfunction()));
	QObject::connect(ui.actionDisplayHKS, SIGNAL(triggered()), this, SLOT(displayHKS()));
	QObject::connect(ui.actionDisplayHK, SIGNAL(triggered()), this, SLOT(displayHK()));
	QObject::connect(ui.actionDisplayMHW, SIGNAL(triggered()), this, SLOT(displayMHW()));
	QObject::connect(ui.actionDisplayMHWS, SIGNAL(triggered()), this, SLOT(displayMHWS()));
	QObject::connect(ui.actionDisplayBiharmonic, SIGNAL(triggered()), this, SLOT(displayBiharmonic()));
	QObject::connect(ui.actionDisplaySimilarityMap, SIGNAL(triggered()), this, SLOT(displaySimilarityMap()));

	QObject::connect(ui.actionExperimental, SIGNAL(triggered()), this, SLOT(displayExperimental()));
	QObject::connect(ui.actionMeanCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureMean()));
	QObject::connect(ui.actionGaussCurvature, SIGNAL(triggered()), this, SLOT(displayCurvatureGauss()));
	QObject::connect(ui.actionDiffPosition, SIGNAL(triggered()), this, SLOT(displayDiffPosition()));

	////////	Task	////////
	QObject::connect(ui.actionTaskRegistration, SIGNAL(triggered()), this, SLOT(setTaskRegistration()));
	QObject::connect(ui.actionTaskEditing, SIGNAL(triggered()), this, SLOT(setTaskEditing()));

	////////	Register	////////
	QObject::connect(ui.actionRegisterAutomatic, SIGNAL(triggered()), this, SLOT(registerAutomatic()));
	QObject::connect(ui.actionBuildHierarchy, SIGNAL(triggered()), this, SLOT(buildHierarchy()));
	QObject::connect(ui.actionDetectFeatures, SIGNAL(triggered()), this, SLOT(detectFeatures()));
	QObject::connect(ui.actionMatchFeatures, SIGNAL(triggered()), this, SLOT(matchFeatures()));
	QObject::connect(ui.actionRegisterStep, SIGNAL(triggered()), this, SLOT(registerStep()));
	QObject::connect(ui.actionRegisterFull, SIGNAL(triggered()), this, SLOT(registerFull()));

}

bool QZGeometryWindow::initialize()
{
	qout.output("******** Welcome ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);

#ifdef NDEBUG
	std::string mesh_list_name = g_configMgr.getConfigValue("MESH_LIST_NAME");
#else NDEBUG
	std::string mesh_list_name = g_configMgr.getConfigValue("MESH_LIST_NAME_DEBUG");
#endif

	int eng_open_time = time_call( [&](){
		m_ep = engOpen("\0");
	});

	if (!m_ep)
	{
		qout.output("Can't start MATLAB engine!", OUT_MSGBOX);
		return false;
	}
	else qout.output(QString().sprintf("Matlab engine initialized! (%fs)", eng_open_time/1000.), OUT_TERMINAL);
	qout.output('*', 24);
	
	setDisplayMesh();
	setEditModeMove();

	//// ---- load meshes ---- ////
	num_meshes = g_configMgr.getConfigValueInt("NUM_PRELOAD_MESHES");
	
	ifstream meshfiles(mesh_list_name);
	if (!meshfiles)
	{
		qout.output("Cannot open " + mesh_list_name, OUT_MSGBOX);
		return false;
	}
	deque<string> vMeshFiles;
	while (!meshfiles.eof())
	{
		string meshFileName;
		getline(meshfiles, meshFileName);
		if (meshFileName == "") continue;
		if (meshFileName[0] == '#') continue;
		else vMeshFiles.push_back(meshFileName);
	}
	meshfiles.close();
	if (vMeshFiles.size() < num_meshes)
	{
		qout.output("Not enough meshes found!", OUT_MSGBOX);
		return false;
	}

	CStopWatch timer;
	timer.startTimer();
	Concurrency::parallel_for(0, num_meshes, [&](int obj)
	{
		CMesh& mesh = m_mesh[obj];
		mesh.Load(vMeshFiles[obj]);
		mesh.scaleEdgeLenToUnit();
		mesh.gatherStatistics();
	});
	timer.stopTimer();
	cout << "Time to load meshes: " << timer.getElapsedTime() << "s" << endl;

	for (int obj = 0; obj < num_meshes; ++obj)
	{
		CMesh& mesh = m_mesh[obj];
		Vector3D center = mesh.getCenter(), bbox = mesh.getBoundingBox();
		qout.output(QString().sprintf("Load mesh: %s; Size: %d", mesh.getMeshName().c_str(), mesh.getVerticesNum()), OUT_CONSOLE);
		qout.output(QString().sprintf("Center: (%f,%f,%f)\nDimension: (%f,%f,%f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_CONSOLE);
		vMP[obj].init(&mesh, m_ep);
		vRS[obj].mesh_color = preset_colors[obj];
		ui.glMeshWidget->addMesh(&vMP[obj], &vRS[obj]);
		mesh_valid[obj] = true;
	}

	int init_laplacian_type = g_configMgr.getConfigValueInt("INIT_LAPLACIAN_TYPE");
	if (init_laplacian_type >= 0 && init_laplacian_type < LaplacianEnd)
	{
		timer.startTimer();
		switch((LaplacianType)init_laplacian_type)
		{
		case Umbrella:
			computeLaplacian1();
			break;
		case CotFormula:
			computeLaplacian2();
			break;
		case Anisotropic1:
			computeLaplacian3();
			break;
		case Anisotropic2:
			computeLaplacian4();
			break;
		case IsoApproximate:
			computeLaplacian5();
			break;
		}
		timer.stopTimer();
		cout << "Time to decompose initial Laplacian: " << timer.getElapsedTime() << "(s)" << endl;
//		qout.output("Non-zeros of Laplacian: " + Int2String(vMP[0].vMeshLaplacian[init_laplacian_type].getNonzeroNum()));
//		qout.output("Non-zeros of Laplacian: " + Int2String(vMP[1].mLaplacian.getNonzeroNum()));
	}
	
	// ---- update ui ---- //
	if (mesh_valid[0])
	{
		ui.glMeshWidget->fieldView(m_mesh[0].getCenter(), m_mesh[0].getBoundingBox());
		ui.spinBox1->setMinimum(0);
		ui.spinBox1->setMaximum(m_mesh[0].getVerticesNum()-1);
		ui.horizontalSlider1->setMinimum(0);
		ui.horizontalSlider1->setMaximum(m_mesh[0].getVerticesNum()-1);
	}
	if (mesh_valid[1])
	{
		ui.spinBox1->setValue(0);	
		ui.spinBox2->setMinimum(0);
		ui.spinBox2->setMaximum(m_mesh[1].getVerticesNum()-1);
		ui.horizontalSlider2->setMinimum(0);
		ui.horizontalSlider2->setMaximum(m_mesh[1].getVerticesNum()-1);
		ui.spinBox2->setValue(0);
	}

	ui.spinBoxParameter->setMinimum(0);
	ui.spinBoxParameter->setMaximum(2 * PARAMETER_SLIDER_CENTER);
	ui.horizontalSliderParamter->setMinimum(0);
	ui.horizontalSliderParamter->setMaximum(2 * PARAMETER_SLIDER_CENTER);
	ui.spinBoxParameter->setValue(PARAMETER_SLIDER_CENTER);

	vRS[0].selected = true;
	vRS[1].selected = false; 

	if (g_task == TASK_REGISTRATION)
	{
		shapeMatcher.initialize(&vMP[0], &vMP[1], m_ep);
		string rand_data_file = g_configMgr.getConfigValue("RAND_DATA_FILE");
		shapeMatcher.readInRandPair(rand_data_file);
		ui.glMeshWidget->setShapeMatcher(&shapeMatcher);
	}

//	evalDistance();
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

	case Qt::Key_BracketLeft:
		showFiner();
		break;

	case Qt::Key_BracketRight:
		showCoarser();
		break;

	case Qt::Key_C:
		clone();
		break;

	case Qt::Key_F:
		if (event->modifiers() & Qt::AltModifier)
			toggleShowFeatures();
		break;

	case Qt::Key_L:
		if (event->modifiers() & Qt::AltModifier)
			toggleShowMatchingLines();
		break;

	case Qt::Key_R:
		if (event->modifiers() & Qt::AltModifier)
		{
			toggleShowRefPoint();
		}
		else repeatOperation();
		break;

	case Qt::Key_S:
		if (event->modifiers() & Qt::AltModifier)
			toggleShowSignature();
		break;

	case Qt::Key_W:	// switch between display mode
	{
		if (vRS[0].displayType == RenderSettings::Mesh)
			setDisplayWireframe();
		else if (vRS[0].displayType == RenderSettings::Wireframe)
			setDisplayPointCloud();
		else setDisplayMesh();
		break;
	}	
		
	case Qt::Key_M:
		setEditModeMove();
		break;

	case Qt::Key_P:
		setEditModePick();
		if (!ui.glMeshWidget->m_bShowRefPoint)
			toggleShowRefPoint();
		break;

	case Qt::Key_D:
		setEditModeDrag();
		break;

	case Qt::Key_E:
		if (deformType == Simple)
			deformSimple();
		else if (deformType == SGW)
			deformSGW();
		else if (deformType == Laplace)
			deformLaplace();
		break;

	case Qt::Key_G:
		this->reconstructSGW();
		break;

	case Qt::Key_Minus:
		this->vMP[0].constrain_weight /= 2;
		qout.output("Constrain weight: " + QString::number(vMP[0].constrain_weight));
		break;

	case Qt::Key_Equal:
		this->vMP[0].constrain_weight *= 2;
		qout.output("Constrain weight: " + QString::number(vMP[0].constrain_weight));
		break;

	case Qt::Key_X:
		if (event->modifiers() & Qt::ShiftModifier)
			refMove.xMove--;
		else refMove.xMove++;
//		updateReferenceMove();
		break;

	case Qt::Key_Y:
		if (event->modifiers() & Qt::ShiftModifier)
			refMove.yMove--;
		else refMove.yMove++;
//		updateReferenceMove();
		break;

	case Qt::Key_Z:
		if (event->modifiers() & Qt::ShiftModifier)
			refMove.zMove--;
		else refMove.zMove++;
//		updateReferenceMove();
		break;

	default: QWidget::keyPressEvent(event);
	}
}

void QZGeometryWindow::computeSGWSFeatures()
{
	vector<double> vTimes;
	vTimes.push_back(10);
	vTimes.push_back(30);
	vTimes.push_back(90);
	vTimes.push_back(270);

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i])
			vMP[i].computeKernelSignatureFeatures(vTimes, SGW_KERNEL);
	}

	if (!ui.glMeshWidget->m_bShowFeatures)
		toggleShowFeatures();
}

void QZGeometryWindow::computeSGW()
{
	CStopWatch timer;
	timer.startTimer();

	double timescales[] = {20}; //{5, 10, 20, 40};
	int nScales = sizeof (timescales) / sizeof(double);
	vector<double> vTimes;
	vTimes.resize(nScales);
	std::copy(timescales, timescales + nScales, vTimes.begin());

	vMP[0].computeSGW(vTimes, &transferFunc1, true, &transferScalingFunc1);

	timer.stopTimer();
	qout.output(QString("Time for compute SGW: ") + QString::number(timer.getElapsedTime()));
}

void QZGeometryWindow::deformSimple()
{
	int activeHandle = vMP[0].active_handle; 
	vector<int> vHandle;
	vHandle.push_back(activeHandle);
	vector<Vector3D> vHandlePos;
	vHandlePos.push_back(vMP[0].mHandles[activeHandle]);
	vector<int> vFree = vMP[0].getMesh_const()->getNeighborVertexIndex(activeHandle, DEFAULT_DEFORM_RING);
	vector<Vector3D> vNewPos;

	vMP[0].deform(vHandle, vHandlePos, vFree, vNewPos, Simple);
	m_mesh[1].setVertexCoordinates(vFree, vNewPos);
	m_mesh[1].setVertexCoordinates(vHandle, vHandlePos);

	deformType = Simple;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::deformLaplace()
{
	vector<int> vHandle;
	vector<Vector3D> vHandlePos;
	vector<int> vFree;
	vector<Vector3D> vNewPos;

	// 	int activeHandle = vMP[0].active_handle; 	
	// 	vHandle.push_back(activeHandle);	
	// 	vHandlePos.push_back(vMP[0].mHandles[activeHandle]);
	//	vFree = vMP[0].getMesh()->getNeighboringVertex(activeHandle, DEFAULT_DEFORM_RING);

	std::set<int> sFreeIdx;
	for (auto iter = vMP[0].mHandles.begin(); iter!= vMP[0].mHandles.end(); ++iter)
	{
		vHandle.push_back(iter->first);
		vHandlePos.push_back(iter->second);

		vector<int> vNeighbor = vMP[0].getMesh_const()->getNeighborVertexIndex(iter->first, DEFAULT_DEFORM_RING);
		sFreeIdx.insert(vNeighbor.begin(), vNeighbor.end());
	}
	vFree.insert(vFree.begin(), sFreeIdx.begin(), sFreeIdx.end());

	try
	{
		vMP[0].deform(vHandle, vHandlePos, vFree, vNewPos, Laplace);
		m_mesh[1].setVertexCoordinates(vFree, vNewPos);
		m_mesh[1].setVertexCoordinates(vHandle, vHandlePos);
	}
	catch (runtime_error* e)
	{
		qout.output(e->what(), OUT_MSGBOX);
	}


	deformType = Laplace;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::deformSGW()
{
	if (!vMP[0].isSGWComputed())
		this->computeSGW();

	vector<Vector3D> vHandlePos;
	vector<int> vHandle;
	vector<Vector3D> vNewPos;
	vector<int> vFree;

	std::set<int> sFreeIdx;
	for (auto iter = vMP[0].mHandles.begin(); iter!= vMP[0].mHandles.end(); ++iter)
	{
		vHandle.push_back(iter->first);
		vHandlePos.push_back(iter->second);

		vector<int> vNeighbor = vMP[0].getMesh_const()->getNeighborVertexIndex(iter->first, DEFAULT_DEFORM_RING);
		sFreeIdx.insert(vNeighbor.begin(), vNeighbor.end());
	}
	vFree.insert(vFree.begin(), sFreeIdx.begin(), sFreeIdx.end());	

	try
	{
		vMP[0].deform(vHandle, vHandlePos, vFree, vNewPos, SGW);
		m_mesh[1].setVertexCoordinates(vFree, vNewPos);
		m_mesh[1].setVertexCoordinates(vHandle, vHandlePos);
	}
	catch (runtime_error* e)
	{
		qout.output(e->what(), OUT_MSGBOX);
	}

	deformType = SGW;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::selectObject( int index )
{
	QString text = ui.boxObjSelect->itemText(index);
	qout.output("Selected object(s): " + text);
	if (text == "1") 
	{
		vRS[0].selected = true;
		vRS[1].selected = false; 
	}
	else if (text == "2")
	{
		vRS[0].selected = false;
		vRS[1].selected = true; 
	}
	else if (text == "All")
	{
		vRS[0].selected = true;
		vRS[1].selected = true; 
	}
	else if (text == "None")
	{
		vRS[0].selected = false;
		vRS[1].selected = false; 
	}
}

void QZGeometryWindow::setRefPoint1( int vn )
{
	vMP[0].setRefPointIndex(vn);
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	updateReferenceMove(0);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setRefPoint2( int vn )
{
	vMP[1].setRefPointIndex(vn);
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	updateReferenceMove(1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setCommonParameter( int p )
{
	m_commonParameter = p;

	if (current_operation == Compute_HKS || current_operation == Compute_HK 
		|| current_operation == Compute_MHWS || current_operation == Compute_MHW
		|| current_operation == Compute_SGWS || current_operation == Compute_SGW)
	{
		double time_scale;
		if (m_commonParameter <= PARAMETER_SLIDER_CENTER) 
			time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)m_commonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
		else 
			time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(m_commonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 
		qout.output(QString().sprintf("HKS timescale %f", time_scale), OUT_STATUS);
	}
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
	ui.actionShowColorLegend->setChecked(bChecked);

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

void QZGeometryWindow::toggleDrawMatching(bool show)
{
	bool bToShow = !ui.glMeshWidget->m_bDrawMatching;
	ui.glMeshWidget->m_bDrawMatching = bToShow;
	ui.actionDrawMatching->setChecked(bToShow);

	ui.glMeshWidget->update();

}

void QZGeometryWindow::toggleShowMatchingLines(bool show)
{
	bool bToShow = !ui.glMeshWidget->m_bShowCorrespondenceLine;
	ui.glMeshWidget->m_bShowCorrespondenceLine = bToShow;
	ui.actionShowMatchingLines->setChecked(bToShow);

	ui.glMeshWidget->update();

}

void QZGeometryWindow::toggleDrawRegistration( bool show /*= false*/ )
{
	bool bToShow = !ui.glMeshWidget->m_bDrawRegistration;
	ui.glMeshWidget->m_bDrawRegistration = bToShow;
	ui.actionDrawRegistration->setChecked(bToShow);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayPointCloud()
{
	ui.actionDisplayPointCloud->setChecked(true);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(false);

	vRS[0].displayType = vRS[1].displayType = RenderSettings::PointCloud;
	vRS[0].glPolygonMode = vRS[1].glPolygonMode = GL_POINT;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayWireframe()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(true);
	ui.actionDisplayMesh->setChecked(false);

	vRS[0].displayType = vRS[1].displayType = RenderSettings::Wireframe;
	vRS[0].glPolygonMode = vRS[1].glPolygonMode = GL_LINE;
	ui.glMeshWidget->update();

}

void QZGeometryWindow::setDisplayMesh()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(true);

	vRS[0].displayType = vRS[1].displayType = RenderSettings::Mesh;
	vRS[0].glPolygonMode = vRS[1].glPolygonMode = GL_FILL;
	ui.glMeshWidget->update();

}

void QZGeometryWindow::displayEigenfunction()
{
	int select_eig = (m_commonParameter - 49 >= 1) ? (m_commonParameter - 49) : 1;

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed()) 
		{
			DifferentialMeshProcessor& mp = vMP[i];
			vRS[i].normalizeSignatureFrom(mp.getMHB().m_func[select_eig].m_vec);
			vRS[i].showColorSignature = true;
		}
	}

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output("Show eigenfunction" + Int2String(select_eig));
}

void QZGeometryWindow::displayExperimental()
{
	DifferentialMeshProcessor& mp = vMP[0];

 	vector<double> vExp;
 	mp.computeExperimentalWavelet(vExp, 30); 

 	vRS[0].normalizeSignatureFrom(vExp);

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output(QString().sprintf("Show MHW from vertex #%d", mp.getRefPointIndex()));
	
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

	vRS[0].bandCurveSignatureFrom(vCurvature, 0, 1);

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

	vRS[0].bandCurveSignatureFrom(vCurvature, -1, 1);

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
	qout.output("Show Mean Curvature");
}

void QZGeometryWindow::displayDiffPosition()
{
	assert(m_mesh[0].getVerticesNum() == m_mesh[1].getVerticesNum());
	int size = m_mesh[0].getVerticesNum();
	vector<double> vDiff;
	vDiff.resize(size);

	for (int i = 0; i < m_mesh[0].getVerticesNum(); ++i)
	{
		vDiff[i] = (m_mesh[0].getVertex_const(i)->getPosition() - m_mesh[1].getVertex_const(i)->getPosition()).length() / m_mesh[0].getAvgEdgeLength();
	}

	vRS[0].normalizeSignatureFrom(vDiff);
	
	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();
	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::updateReferenceMove( int obj )
{
	DifferentialMeshProcessor& mp = vMP[obj]; 

	double unitMove = (mp.getMesh_const()->getBoundingBox().x + mp.getMesh_const()->getBoundingBox().y + mp.getMesh_const()->getBoundingBox().z)/300.0;
	Vector3D originalPos = mp.getMesh_const()->getVertex_const(mp.getRefPointIndex())->getPosition();
	
	mp.setRefPointPosition(originalPos.x + unitMove * refMove.xMove,
						   originalPos.y + unitMove * refMove.yMove,
						   originalPos.z + unitMove * refMove.zMove);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::clone()
{
	m_mesh[1].cloneFrom(m_mesh[0]);
	m_mesh[1].gatherStatistics();

	vMP[1].init(&m_mesh[1], m_ep);
	vRS[1].mesh_color = preset_colors[1];
	ui.glMeshWidget->addMesh(&vMP[1], &vRS[1]);
	qout.output(QString().sprintf("Mesh %s constructed! Size: %d", m_mesh[1].getMeshName().c_str(), m_mesh[1].getVerticesNum()));
	
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

void QZGeometryWindow::reconstructMHB()
{
	double ratio = min((double)m_commonParameter/PARAMETER_SLIDER_CENTER, 1.0);
	int nEig = vMP[0].getMHB().m_nEigFunc * ratio;

	vector<double> vx, vy, vz;
	vMP[0].reconstructByMHB(nEig, vx, vy, vz);
	m_mesh[1].setVertexCoordinates(vx, vy, vz);
	ui.glMeshWidget->update();

	double errorSum(0);
	for (int i = 0; i < m_mesh[0].getVerticesNum(); ++i)
	{
		errorSum += (m_mesh[0].getVertex_const(i)->getPosition() - m_mesh[1].getVertex_const(i)->getPosition()).length();
	}
	errorSum /= m_mesh[0].getVerticesNum() * m_mesh[0].getAvgEdgeLength();
	qout.output("MHB reconstruction with " + Int2String(nEig) + " MHBs.");
	qout.output("Average position error: " + QString::number(errorSum));

	ui.glMeshWidget->update();
}

void QZGeometryWindow::reconstructSGW()
{
	if (!vMP[0].isSGWComputed())
		this->computeSGW();
	
	vector<double> vx, vy, vz;
	CStopWatch timer;
	timer.startTimer();
	try
	{
		vMP[0].reconstructBySGW(vx, vy, vz, true);
	}
	catch (logic_error* e)
	{
		qout.output(e->what(), OUT_MSGBOX);
		return;
	}
	timer.stopTimer();
	qout.output(QString("SGW reconstruct time: ") + QString::number(timer.getElapsedTime()));

	m_mesh[1].setVertexCoordinates(vx, vy, vz);

	{
		int debugIdx = (*vMP[0].mHandles.begin()).first;
		qout.output("Original pos: " + std::string(m_mesh[0].getVertex(debugIdx)->getPosition()));
		if (!vMP[0].mHandles.empty())
			qout.output("Handle pos: " + std::string((*vMP[0].mHandles.begin()).second));
		qout.output("Deformed pos: " + std::string(m_mesh[1].getVertex(debugIdx)->getPosition()));
	}

	double errorSum(0);
	for (int i = 0; i < m_mesh[0].getVerticesNum(); ++i)
	{
		errorSum += (m_mesh[0].getVertex_const(i)->getPosition() - m_mesh[1].getVertex_const(i)->getPosition()).length();
	}
	errorSum /= m_mesh[0].getVerticesNum() * m_mesh[0].getAvgEdgeLength();
	qout.output("Average position error: " + QString::number(errorSum));

	ui.glMeshWidget->update();
}

void QZGeometryWindow::filterExperimental()
{
	vector<double> vx, vy, vz;
	vMP[0].filterBySGW(vx, vy, vz);
	m_mesh[1].setVertexCoordinates(vx, vy, vz);

	double errorSum(0);
	for (int i = 0; i < m_mesh[0].getVerticesNum(); ++i)
	{
		errorSum += (m_mesh[0].getVertex_const(i)->getPosition() - m_mesh[1].getVertex_const(i)->getPosition()).length();
	}
	errorSum /= m_mesh[0].getVerticesNum() * m_mesh[0].getAvgEdgeLength();
	qout.output("Average position error: " + QString::number(errorSum));

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayNeighborVertices()
{
	int ring = (m_commonParameter > PARAMETER_SLIDER_CENTER) ? (m_commonParameter-PARAMETER_SLIDER_CENTER) : 1;

	int ref = vMP[0].getRefPointIndex();
	std::vector<int> vn = vMP[0].getMesh_const()->getNeighborVertexIndex(ref, ring);
//	std::vector<int> vn = vMP[0].getMesh()->getRingVertex(ref, ring);
	MeshFeatureList *mfl = new MeshFeatureList;

	for (auto iter = vn.begin(); iter != vn.end(); ++iter)
	{
		mfl->getFeatureVector()->push_back(new MeshFeature(*iter));
		mfl->setIDandName(FEATURE_NEIGHBORS, "Neighbors");
	}
	vMP[0].addProperty(mfl);

	vMP[0].setActiveFeatures(mfl->getFeatureVector());
	
	if (!ui.actionShowFeatures->isChecked())
		toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeHKS()
{
	double time_scale;
	if (m_commonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)m_commonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(m_commonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("Heat Kernel timescale: %f", time_scale));

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed())
		{
			DifferentialMeshProcessor& mp = vMP[i];
			mp.computeKernelSignature(time_scale, HEAT_KERNEL);
		}
	}
	
	displayHKS();

	current_operation = Compute_HKS;
}

void QZGeometryWindow::computeHK()
{
	double time_scale;
	if (m_commonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)m_commonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(m_commonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("Heat Kernel timescale: %f", time_scale));

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed())
		{
			DifferentialMeshProcessor& mp = vMP[i];
			int refPoint = mp.getRefPointIndex();
			mp.computeKernelDistanceSignature(time_scale, HEAT_KERNEL, refPoint);
		}
	}

	displayHK();

	current_operation = Compute_HK;
}

void QZGeometryWindow::displayHKS()
{
	displaySignature(SIGNATURE_HKS);
}

void QZGeometryWindow::displayHK()
{
	displaySignature(SIGNATURE_HK);
}

void QZGeometryWindow::repeatOperation()
{
	switch(current_operation)
	{
	case Compute_HKS:
		computeHKS();
		break;
	case Compute_HK:
		computeHK();
		break;
	case Compute_MHWS:
		computeMHWS();
		break;
	case Compute_MHW:
		computeMHW();
		break;
	case Compute_SGWS:
		computeSGWS();
		break;
	}
}

void QZGeometryWindow::computeHKSFeatures()
{
	vector<double> vTimes;
	vTimes.push_back(10);
	vTimes.push_back(30);
	vTimes.push_back(90);
	vTimes.push_back(270);

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i])
			vMP[i].computeKernelSignatureFeatures(vTimes, HEAT_KERNEL);
	}
	
	if (!ui.glMeshWidget->m_bShowFeatures)
		toggleShowFeatures();
}

void QZGeometryWindow::computeMHWFeatures()
{
	vector<double> vTimes;
	vTimes.push_back(10);
	vTimes.push_back(30);
	vTimes.push_back(90);
	vTimes.push_back(270);

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i])
			vMP[i].computeKernelSignatureFeatures(vTimes, MHW_KERNEL);
	}

	if (!ui.glMeshWidget->m_bShowFeatures)
		toggleShowFeatures();

}

void QZGeometryWindow::computeMHWS()
{
	double time_scale;
	if (m_commonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)m_commonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(m_commonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("MHW timescale: %f", time_scale));

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed())
		{
			DifferentialMeshProcessor& mp = vMP[i];
			mp.computeKernelSignature(time_scale, MHW_KERNEL);
		}
	}

	displayMHWS();

	current_operation = Compute_MHWS;
}

void QZGeometryWindow::displayMHWS()
{
	displaySignature(SIGNATURE_MHWS);
}

void QZGeometryWindow::computeSGWS()
{
	double time_scale;
	if (m_commonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)m_commonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(m_commonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("Spectral Graph Wavelet timescale: %f", time_scale));

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed())
		{
			DifferentialMeshProcessor& mp = vMP[i];
			mp.computeKernelSignature(time_scale, SGW_KERNEL);
		}
	}

	displaySGWS();

	current_operation = Compute_SGWS;
}

void QZGeometryWindow::displaySGWS()
{
	displaySignature(SIGNATURE_SGWS);
}

void QZGeometryWindow::displaySignature( int signatureID )
{
	for (int i = 0; i < 2; ++i)
	{
		DifferentialMeshProcessor& mp = vMP[i];
		MeshProperty* vs = mp.retrievePropertyByID(signatureID);
		if (vs != NULL)
			vRS[i].normalizeSignatureFrom(dynamic_cast<MeshFunction*>(vs)->getMeshFunction_const());
//			vRS[i].logNormalizeSignatureFrom(dynamic_cast<MeshFunction*>(vs)->getMeshFunction_const());
	}

	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();

	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeMHW()
{
	double time_scale;
	if (m_commonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)m_commonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(m_commonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("MHW timescale: %f", time_scale));

	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed())
		{
			DifferentialMeshProcessor& mp = vMP[i];
			int refPoint = mp.getRefPointIndex();
			mp.computeKernelDistanceSignature(time_scale, MHW_KERNEL, refPoint);
		}
	}

	displayMHW();

	current_operation = Compute_MHW;
}

void QZGeometryWindow::displayMHW()
{
	displaySignature(SIGNATURE_MHW);
}

void QZGeometryWindow::setTaskRegistration()
{
	g_task = TASK_REGISTRATION;
	ui.actionTaskRegistration->setChecked(true);
	ui.actionTaskEditing->setChecked(false);
}

void QZGeometryWindow::setTaskEditing()
{
	g_task = TASK_EDITING;
	ui.actionTaskRegistration->setChecked(false);
	ui.actionTaskEditing->setChecked(true);
}

void QZGeometryWindow::registerAutomatic()
{
	this->buildHierarchy();
	this->detectFeatures();
	this->matchFeatures();
}

void QZGeometryWindow::buildHierarchy()
{
	int nLevel = g_configMgr.getConfigValueInt("HIERARCHY_LEVEL");
	double ratio = g_configMgr.getConfigValueDouble("CONTRACTION_RATIO");
	std::string log_filename = g_configMgr.getConfigValue("HIERARCHY_OUTPUT_FILE");
	ofstream ostr(log_filename.c_str(), ios::trunc);

	qout.output("-- Build hierarchy --");
	shapeMatcher.constructPyramid(nLevel, ratio, ostr);
	qout.output("Mesh hierarchy constructed!");

	ostr.close();
}

void QZGeometryWindow::detectFeatures()
{
	qout.output("-- Detect initial features --");
	shapeMatcher.detectFeatures(0, 2, 4, DiffusionShapeMatcher::DEFAULT_FEATURE_TIMESCALE, DiffusionShapeMatcher::DEFAULT_T_MULTIPLIER, DiffusionShapeMatcher::DEFAULT_EXTREAMA_THRESH);
	shapeMatcher.detectFeatures(1, 2, 4, DiffusionShapeMatcher::DEFAULT_FEATURE_TIMESCALE, DiffusionShapeMatcher::DEFAULT_T_MULTIPLIER, DiffusionShapeMatcher::DEFAULT_EXTREAMA_THRESH);

	qout.output("Multi-scale mesh features detected!");
	qout.output(QString().sprintf("Mesh1 features#: %d; Mesh2 features#: %d", shapeMatcher.getSparseFeatures(0).size(), shapeMatcher.getSparseFeatures(1).size()));

	const vector<HKSFeature>& vf1 = shapeMatcher.vFeatures[0], &vf2 = shapeMatcher.vFeatures[1];

	int count_possible = 0;
	for (auto iter = vf1.begin(); iter != vf1.end(); ++iter)
	{
		for (auto iter2 = vf2.begin(); iter2 != vf2.end(); ++iter2)
		{
			if (std::abs(iter->m_index - iter2->m_index) <= 1 || m_mesh[1].isInNeighborRing(iter->m_index, iter2->m_index, 1))
			{
				count_possible++;
				break;
			}
		}
	}
	qout.output(QString().sprintf("-- Valid detections: %d", count_possible));

	if (!ui.glMeshWidget->m_bShowFeatures)
		toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::matchFeatures()
{
	bool force_matching = (1 == g_configMgr.getConfigValueInt("FORCE_MATCHING"));
	
	bool use_tensor = (g_configMgr.getConfigValueInt("USE_TENSOR_MATCHING") == 1);
	double matching_thresh_1 = g_configMgr.getConfigValueDouble("MATCHING_THRESH_1");
	double matching_thresh_2 = g_configMgr.getConfigValueDouble("MATCHING_THRESH_2");
	double tensor_matching_timescasle = g_configMgr.getConfigValueDouble("TENSOR_MATCHING_TIMESCALE");
	std::string log_filename = g_configMgr.getConfigValue("MATCH_OUTPUT_FILE");

	qout.output("-- Match initial features --");

	ofstream ofstr(log_filename.c_str(), ios::trunc);
	CStopWatch timer;
	timer.startTimer();

	if (use_tensor)
	{
		shapeMatcher.matchFeaturesTensor(ofstr, tensor_matching_timescasle, matching_thresh_2);
	}
	else
	{
		shapeMatcher.matchFeatures(ofstr, matching_thresh_1);
	}

	timer.stopTimer();
	ofstr.close();
	qout.output(QString().sprintf("Initial features matched! Matched#:%d. Time elapsed:%f", shapeMatcher.getMatchedFeaturesResults(shapeMatcher.getAlreadyMatchedLevel()).size(), timer.getElapsedTime()));

	if (force_matching)
	{
		string string_override = "";
		if (m_mesh[0].getMeshName() == "horse0")
		{
			string_override = g_configMgr.getConfigValue("HORSE0_FEATURE_OVERRIDE");
		}
		
		if (string_override != "")
		{
			vector<int> idx_override = splitStringToInt(string_override);
			vector<MatchPair> vmp;
			for (auto iter = begin(idx_override); iter != end(idx_override); ++iter)
			{
				vmp.push_back(MatchPair(*iter, *iter));
			}
			shapeMatcher.forceInitialAnchors(vmp);
			
			qout.output("!!Matched anchors manually assigned!!");
		}
	}

	if (!ui.glMeshWidget->m_bDrawMatching)
		toggleDrawMatching();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::registerStep()
{
	string log_filename = g_configMgr.getConfigValue("REGISTER_OUTPUT_FILE");

	qout.output(QString().sprintf("-- Register level %d --", shapeMatcher.getAlreadyRegisteredLevel() - 1));
	
	ofstream ofstr;	
	if (shapeMatcher.getAlreadyRegisteredLevel() == shapeMatcher.getTotalRegistrationLevels())
		ofstr.open(log_filename.c_str(), ios::trunc);
	else ofstr.open(log_filename.c_str(), ios::app);

	double time_elapsed = time_call([&]{
		shapeMatcher.refineRegister2(ofstr);
	}) / 1000.0;

	int level = shapeMatcher.getAlreadyRegisteredLevel();
	const vector<MatchPair>& vf = shapeMatcher.getMatchedFeaturesResults(shapeMatcher.getAlreadyMatchedLevel());
	const vector<MatchPair>& vr = shapeMatcher.getRegistrationResults(shapeMatcher.getAlreadyRegisteredLevel());
	
	qout.output(QString().sprintf("Registration level %d finished! Time elapsed:%f\n-Features Matced:%d; Registered:%d",
		                        level, time_elapsed, vf.size(), vr.size()));
	qout.output(QString().sprintf("Registered ratio: %f (%d/%d)", 
		                        double(vr.size())/shapeMatcher.getMesh(0, level)->getVerticesNum(),
								vr.size(), shapeMatcher.getMesh(0, level)->getVerticesNum()));
	/* ---- evaluation ---- */
	double vError[3];
	for(int k = 0; k < 2; ++k) { 
		vError[k] = DiffusionShapeMatcher::evaluateDistortion(vr, &m_mesh[0], &m_mesh[1], shapeMatcher.m_randPairs, 200 * k);
	}

	qout.output(QString().sprintf("Registration error: %f, %f", vError[0], vError[1]));

	if (!ui.glMeshWidget->m_bDrawRegistration)
		toggleDrawRegistration();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::registerFull()
{
	while(shapeMatcher.getAlreadyMatchedLevel() > 0)
		registerStep();
}

void QZGeometryWindow::showFiner()
{
	if (ui.glMeshWidget->m_nMeshLevel > 0)
		ui.glMeshWidget->m_nMeshLevel--;

	qout.output("Display mesh level " + QString::number(ui.glMeshWidget->m_nMeshLevel));
	ui.glMeshWidget->update();
}

void QZGeometryWindow::showCoarser()
{
	if (ui.glMeshWidget->m_nMeshLevel < shapeMatcher.getPyramidLevels()-1)
		ui.glMeshWidget->m_nMeshLevel++;

	qout.output("Display mesh level " + QString::number(ui.glMeshWidget->m_nMeshLevel));
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeBiharmonic()
{
	for (int i = 0; i < 2; ++i)
	{
		if (mesh_valid[i] && vMP[i].isLaplacianDecomposed())
		{
			DifferentialMeshProcessor& mp = vMP[i];
			int refPoint = mp.getRefPointIndex();
			mp.computeBiharmonicDistanceSignature(refPoint);
		}
	}

	displayBiharmonic();
}

void QZGeometryWindow::displayBiharmonic()
{
	displaySignature(SIGNATURE_BIHARMONIC_DISTANCE);
}

void QZGeometryWindow::evalDistance()
{
	CStopWatch timer;
	timer.startTimer();
	Concurrency::parallel_invoke(
		[&](){ cout << "Error geodesic1: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 0) << endl; },
//		[&](){ cout << "Error geodesic2: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl; },
		[&](){ cout << "Error biharmonic1: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_BIHARMONIC, std::vector<double>(), shapeMatcher.m_randPairs, 0) << endl; },
		[&](){ cout << "Error biharmonic2: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_BIHARMONIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl; }
//		[&](){cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 30.), shapeMatcher.m_randPairs, 0) << endl;},
//		cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 90.), shapeMatcher.m_randPairs, 0) << endl;
//		cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 270.), shapeMatcher.m_randPairs, 0) << endl;
	);
//	cout << "Error geodesic1: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl;
	timer.stopTimer();
	cout << "Eval Dist time (ppl): " << timer.getElapsedTime() << endl;

//  timer.startTimer();
//  cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 0) << endl;
// 	cout << "Error biharmonic1: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_BIHARMONIC, std::vector<double>(), shapeMatcher.m_randPairs, 0) << endl;
// 	cout << "Error biharmonic2: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_BIHARMONIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl;
//  timer.stopTimer();
//  cout << "Eval Dist time: " << timer.getElapsedTime() << endl;
}

void QZGeometryWindow::calculateSingleLaplacian( int obj, LaplacianType laplacianType /*= CotFormula*/ )
{
	if (!mesh_valid[obj]) return;
	DifferentialMeshProcessor& mp = vMP[obj];
	const CMesh& mesh = m_mesh[obj];

	if (!mp.vMHB[laplacianType].empty()) 
	{
		mp.setActiveMHB(laplacianType);
		return;
	}

	CStopWatch timer;
	timer.startTimer();

	string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "output/" + mp.getMesh_const()->getMeshName() + ".mhb." + s_idx;
	ifstream ifs(pathMHB.c_str());
	if (ifs && LOAD_MHB_CACHE)	// MHB cache available for the current mesh
	{
		mp.readMHB(pathMHB, laplacianType);
		ifs.close();
	}
	else // need to compute Laplacian and to cache
	{
		int nEig = min(DEFAULT_EIGEN_SIZE, mesh.getVerticesNum()-1);

		mp.decomposeLaplacian(nEig, laplacianType);
		mp.writeMHB(pathMHB, laplacianType);
		qout.output("MHB saved to " + pathMHB);

		std::string pathEVL = "output/" + mp.getMesh_const()->getMeshName() + ".evl";	//dump eigenvalues
		mp.getMHB().dumpEigenValues(pathEVL);
	}

	timer.stopTimer();
	qout.output(QString().sprintf("Laplacian %d decomposed in %f(s)", obj, timer.getElapsedTime()));
	qout.output(QString().sprintf("Min EigVal: %f, Max EigVal: %f", mp.getMHB().m_func.front().m_val, mp.getMHB().m_func.back().m_val));
}

void QZGeometryWindow::computeLaplacian1()
{
	LaplacianType laplacianType = Umbrella;

	Concurrency::parallel_for(0, num_meshes, [&](int obj)
	{
		if (!vMP[obj].vMeshLaplacian[laplacianType].m_bMatrixBuilt)
			cout << "Time to construct Umbrella Laplacian: " << time_call([&](){
				vMP[obj].vMeshLaplacian[laplacianType].constructFromMesh1(&m_mesh[obj]);
			}) / 1000. << "(s)" << endl;
	});

	calculateLaplacians(laplacianType);
}

void QZGeometryWindow::computeLaplacian2()
{
	LaplacianType laplacianType = CotFormula;

	Concurrency::parallel_for(0, num_meshes, [&](int obj)
	{
		if (!vMP[obj].vMeshLaplacian[laplacianType].m_bMatrixBuilt)
			cout << "Time to construct CotFormula Laplacian: " << time_call([&](){
				vMP[obj].vMeshLaplacian[laplacianType].constructFromMesh2(&m_mesh[obj]);
			}) / 1000. << "(s)" << endl;
	});

	calculateLaplacians(laplacianType);
}

void QZGeometryWindow::computeLaplacian3()
{
	LaplacianType laplacianType = Anisotropic1;

	Concurrency::parallel_for(0, num_meshes, [&](int obj)
	{		
		if (!vMP[obj].vMeshLaplacian[laplacianType].m_bMatrixBuilt)
		{
			const CMesh* tm = &m_mesh[obj];
			cout << "Time to construct Anisotropic_1 kernel: " << time_call([&](){
				double para1 = 2 * tm->getAvgEdgeLength() * tm->getAvgEdgeLength();
				double para2 = para1;
				vMP[obj].vMeshLaplacian[laplacianType].constructFromMesh3(&m_mesh[obj], 1, para1, para2);
			}) / 1000. << "(s)" << endl;
		}	
//		cout << "Kernel Sparsity: " << vMP[obj]vMeshLaplacian[Anisotropic1].vSS.size() / double(m_size*m_size) << endl;
	});

	calculateLaplacians(laplacianType);
}

void QZGeometryWindow::computeLaplacian4()
{
	LaplacianType laplacianType = Anisotropic2;
	Concurrency::parallel_for(0, num_meshes, [&](int obj)
	{
		if (!vMP[obj].vMeshLaplacian[laplacianType].m_bMatrixBuilt)
		{
			const CMesh* tm = &m_mesh[obj];
			cout << "Time to construct Anisotropic_2 kernel: " << time_call([&](){
				double para1 = 2 * tm->getAvgEdgeLength() * tm->getAvgEdgeLength();
				double para2 = tm->getAvgEdgeLength() / 2;
				vMP[obj].vMeshLaplacian[laplacianType].constructFromMesh4(&m_mesh[obj], 1, para1, para2);
			}) / 1000. << "(s)" << endl;
//		cout << "Kernel Sparsity: " << vMP[obj]vMeshLaplacian[laplacianType].vSS.size() / double(m_size*m_size) << endl;
		}
	});

	calculateLaplacians(laplacianType);
}

void QZGeometryWindow::calculateLaplacians( LaplacianType laplacianType /*= CotFormula*/ )
{
	if (LOAD_MHB_CACHE && num_meshes == 2)	// special acceleration condition
	{
		DifferentialMeshProcessor& mp1 = vMP[0], &mp2 = vMP[1];
		string s_idx = "0";
		s_idx[0] += (int)laplacianType;
		std::string pathMHB1 = "output/" + m_mesh[0].getMeshName() + ".mhb." + s_idx;
		std::string pathMHB2 = "output/" + m_mesh[1].getMeshName() + ".mhb." + s_idx;

		ifstream ifs1(pathMHB1.c_str()), ifs2(pathMHB2.c_str());

		if (ifs1 && ifs2)	// both MHB can be read from disk
		{
			CStopWatch timer;
			timer.startTimer();

			Concurrency::parallel_invoke( [&](){ mp1.readMHB(pathMHB1, laplacianType); },
										  [&](){ mp2.readMHB(pathMHB2, laplacianType); }
										);
			timer.stopTimer();
			qout.output(QString().sprintf("2 Laplacian decomposed in %f(s)", timer.getElapsedTime()));

			return;
		}		
	}

	for (int obj = 0; obj < num_meshes; ++obj)
	{
		calculateSingleLaplacian(obj, laplacianType);
	}

	for (int l = 0; l < LaplacianEnd; ++l)
		m_actionLaplacians[l]->setChecked(false);
	m_actionLaplacians[laplacianType]->setChecked(true);
}

void QZGeometryWindow::computeSimilarityMap()
{
	Concurrency::parallel_for(0, num_meshes, [&](int obj){
		vMP[obj].computeSimilarityMap(vMP[obj].getRefPointIndex());
	});

	displaySimilarityMap();
}

void QZGeometryWindow::displaySimilarityMap()
{
	displaySignature(SIGNATURE_SIMILARITY_MAP);
}

void QZGeometryWindow::computeSimilarityMap2()
{
	Concurrency::parallel_for(0, num_meshes, [&](int obj){
		vMP[obj].computeSimilarityMap2(vMP[obj].getRefPointIndex());
	});

	displaySimilarityMap();
}

void QZGeometryWindow::computeLaplacian5()
{
	LaplacianType laplacianType = IsoApproximate;
	Concurrency::parallel_for(0, num_meshes, [&](int obj)
	{
		if (!vMP[obj].vMeshLaplacian[laplacianType].m_bMatrixBuilt)
		{
			const CMesh* tm = &m_mesh[obj];
			cout << "Time to construct mesh laplacian: " << time_call([&](){
				vMP[obj].vMeshLaplacian[laplacianType].constructFromMesh5(&m_mesh[obj]);
			}) / 1000. << "(s)" << endl;
			//		cout << "Kernel Sparsity: " << vMP[obj]vMeshLaplacian[laplacianType].vSS.size() / double(m_size*m_size) << endl;
		}
	});

	calculateLaplacians(laplacianType);
}

