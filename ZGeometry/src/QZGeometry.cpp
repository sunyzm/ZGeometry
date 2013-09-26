#include "QZGeometry.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <set>
#include <stdexcept>
#include <ppl.h>
#include <boost/lexical_cast.hpp>
#include <QtWidgets/QMessageBox>
#include <QFileDialog>
#include <QTime>
#include <QProcess>
#include <ZUtil/ZUtil.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/DenseMatrix.h>
#include "global.h"

using namespace std;
using ZUtil::Int2String;
using ZUtil::logic_assert;
using ZUtil::runtime_assert;
using ZGeom::MatlabEngineWrapper;

int     QZGeometryWindow::DEFAULT_EIGEN_SIZE      = 300;
int     QZGeometryWindow::DEFAULT_DEFORM_RING     = 5 ;
int     QZGeometryWindow::LOAD_MHB_CACHE          = 0;
double  QZGeometryWindow::MIN_HK_TIMESCALE        = 1;
double  QZGeometryWindow::DEFUALT_HK_TIMESCALE    = 40.0;
double  QZGeometryWindow::MAX_HK_TIMESCALE        = 2000.0;
double  QZGeometryWindow::PARAMETER_SLIDER_CENTER = 50;
double  QZGeometryWindow::DR_THRESH_INCREMENT     = 0.00001;

QZGeometryWindow::QZGeometryWindow(QWidget *parent,  Qt::WindowFlags flags) 
	: QMainWindow(parent, flags), mEngineWrapper(256)
{
	/* read in configuration parameters from g_configMgr */
	g_configMgr.getConfigValueInt("LOAD_MHB_CACHE", LOAD_MHB_CACHE);
	g_configMgr.getConfigValueDouble("PARAMETER_SLIDER_CENTER", PARAMETER_SLIDER_CENTER);
	g_configMgr.getConfigValueDouble("DEFUALT_HK_TIMESCALE", DEFUALT_HK_TIMESCALE);
	g_configMgr.getConfigValueInt("DEFAULT_EIGEN_SIZE", DEFAULT_EIGEN_SIZE);

	mMeshCount = 0;
	mObjInFocus = -1;
	mCommonParameter = PARAMETER_SLIDER_CENTER;
	current_operation = None;
	deformType = Simple;
	refMove.xMove = refMove.yMove = refMove.zMove = 0;

	/* setup ui and connections */
	ui.setupUi(this);
	ui.centralWidget->setLayout(ui.mainLayout);
	ui.glMeshWidget->setup(&mProcessors, &mRenderManagers, &mShapeMatcher, &mShapeEditor);
	this->makeConnections();
	
	ui.spinBoxParameter->setMinimum(0);
	ui.spinBoxParameter->setMaximum(2 * PARAMETER_SLIDER_CENTER);
	ui.horizontalSliderParamter->setMinimum(0);
	ui.horizontalSliderParamter->setMaximum(2 * PARAMETER_SLIDER_CENTER);
	ui.horizontalSliderParamter->setSliderPosition(PARAMETER_SLIDER_CENTER);

	setDisplayMesh();
	setEditModeMove();

	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);    
}

QZGeometryWindow::~QZGeometryWindow()
{
	for (CMesh* m : mMeshes) delete m;
	for (DifferentialMeshProcessor* p : mProcessors) delete p;
	for (RenderSettings* rs : mRenderManagers) delete rs;

	for (QAction* a : m_actionDisplaySignatures) delete a;
	for (QAction* a : m_actionComputeSimilarities) delete a;
	for (QAction* a : m_actionComputeLaplacians) delete a;

	delete simlaritySignalMapper;
	delete laplacianSignalMapper;
	delete signatureSignalMapper;
}

void QZGeometryWindow::makeConnections()
{	
	/*  actionComputeLaplacians  */
	int laplacianTypeCount = MeshLaplacian::LaplacianTypeCount;
	m_actionComputeLaplacians.resize(laplacianTypeCount);
	laplacianSignalMapper = new QSignalMapper(this);
	for (int t = 0; t < laplacianTypeCount; ++t) {
		m_actionComputeLaplacians[t] = new QAction(QString("Laplacian type ") + QString::number(t), this);
		ui.menuComputeLaplacian->addAction(m_actionComputeLaplacians[t]);
		laplacianSignalMapper->setMapping(m_actionComputeLaplacians[t], t);
		QObject::connect(m_actionComputeLaplacians[t], SIGNAL(triggered()), laplacianSignalMapper, SLOT(map()));
	}
	QObject::connect(laplacianSignalMapper, SIGNAL(mapped(int)), this, SLOT(computeLaplacian(int)));

	/*  actionComputeSimilarities  */
	m_actionComputeSimilarities.resize(SIM_TYPE_COUNT);
	simlaritySignalMapper = new QSignalMapper(this);
	for (int t = 0; t < SIM_TYPE_COUNT; ++t) {
		m_actionComputeSimilarities[t] = new QAction(QString("Similarity type ") + QString::number(t), this);
		ui.menuComputeSimilarityMap->addAction(m_actionComputeSimilarities[t]);
		simlaritySignalMapper->setMapping(m_actionComputeSimilarities[t], t);
		QObject::connect(m_actionComputeSimilarities[t], SIGNAL(triggered()), simlaritySignalMapper, SLOT(map()));
	}
	QObject::connect(simlaritySignalMapper, SIGNAL(mapped(int)), this, SLOT(computeSimilarityMap(int)));
	
	/*  actionDisplaySignatures  */
	signatureSignalMapper = new QSignalMapper(this);
	QObject::connect(signatureSignalMapper, SIGNAL(mapped(int)), this, SLOT(displaySignature(int)));

	////////	file	////////
	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	QObject::connect(ui.actionSaveSignature, SIGNAL(triggered()), this, SLOT(saveSignature()));
	QObject::connect(ui.actionAddMesh, SIGNAL(triggered()), this, SLOT(addMesh()));
	QObject::connect(ui.actionSaveMatching, SIGNAL(triggered()), this, SLOT(saveMatchingResult()));
	QObject::connect(ui.actionLoadMatching, SIGNAL(triggered()), this, SLOT(loadMatchingResult()));

	////////	compute	////////
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(computeEigenfunction()));
	QObject::connect(ui.actionMeanCurvature, SIGNAL(triggered()), this, SLOT(computeCurvatureMean()));
	QObject::connect(ui.actionGaussCurvature, SIGNAL(triggered()), this, SLOT(computeCurvatureGauss()));
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
	QObject::connect(ui.actionRevert, SIGNAL(triggered()), this, SLOT(revert()));
	QObject::connect(ui.actionReconstructMHB, SIGNAL(triggered()), this, SLOT(reconstructMHB()));
	QObject::connect(ui.actionDeformDifferential, SIGNAL(triggered()), this, SLOT(deformDifferential()));
	QObject::connect(ui.actionDeformSimple, SIGNAL(triggered()), this, SLOT(deformSimple()));
	QObject::connect(ui.actionDeformSGW, SIGNAL(triggered()), this, SLOT(deformSGW()));
	QObject::connect(ui.actionDeformLaplace, SIGNAL(triggered()), this, SLOT(deformLaplace()));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
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
	QObject::connect(ui.actionRegisterTest, SIGNAL(triggered()), this, SLOT(registerTest()));
}

bool QZGeometryWindow::initialize(const std::string& mesh_list_name)
{
	qout.output("******** Welcome ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);
	qout.output('*', 24, OUT_CONSOLE);

	try	{
		CStopWatch timer;
		timer.startTimer();
		mEngineWrapper.open();     
		timer.stopTimer("-- Matlab loading time: ", " --");

		switch (g_task)
		{
		case TASK_VIEWING:
			g_configMgr.getConfigValueInt("NUM_PRELOAD_MESHES", mMeshCount);
			if (mMeshCount <= 0) return true;
			break;
		case TASK_REGISTRATION:
			mMeshCount = 2;
			break;
		case TASK_EDITING:
			mMeshCount = 1;
			break;
		default:
			break;
		}

		loadInitialMeshes(mesh_list_name); 

		/* compute and decompose mesh Laplacians */
#if 0		
		int init_laplacian_type = MeshLaplacian::CotFormula;
		g_configMgr.getConfigValueInt("INIT_LAPLACIAN_TYPE", init_laplacian_type);
		if (init_laplacian_type >= 0 && init_laplacian_type < MeshLaplacian::LaplacianTypeCount) { 
			timer.startTimer();
			computeLaplacian(init_laplacian_type);			
			timer.stopTimer("-- Time to decompose initial Laplacians: ", " --");		
		} else {
			std::cout << "Unrecognized Laplacian type; Nothing is done.";
		}
#endif
		computeLaplacian(MeshLaplacian::CotFormula);
		computeLaplacian(MeshLaplacian::SymCot);

		if (g_task == TASK_REGISTRATION) registerPreprocess();
		if (g_task == TASK_EDITING) mShapeEditor.init(mProcessors[0]);

	} catch (std::exception* e) {
		std::cerr << "!!Fatal error in initialization:\n" << e->what() << std::endl;
		return false;
	}

	return true;
}

void QZGeometryWindow::loadInitialMeshes(const std::string& mesh_list_name)
{
	ZUtil::runtime_assert (ZUtil::fileExist(mesh_list_name),  "Cannot open file mesh list file!");
	ifstream meshfiles(mesh_list_name);
	
	std::vector<std::string> vMeshFiles;
	while (!meshfiles.eof()) {
		std::string meshFileName;
		getline(meshfiles, meshFileName);
		if (meshFileName == "") continue;
		if (meshFileName[0] == '#') continue;
		else if (!ZUtil::fileExist(meshFileName))
			throw std::runtime_error("Cannot open file " + meshFileName);		
		vMeshFiles.push_back(meshFileName);
	}
	meshfiles.close();
	if (vMeshFiles.size() < mMeshCount)
		throw std::runtime_error("Not enough meshes in mesh list!");

	allocateStorage(mMeshCount);

	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
		CMesh& mesh = *mMeshes[obj];
		mesh.Load(vMeshFiles[obj]);
		mesh.scaleAreaToVertexNum();
		mesh.gatherStatistics();        
	});	

	for (int obj = 0; obj < mMeshCount; ++obj) {
		CMesh& mesh = *mMeshes[obj];
		Vector3D center = mesh.getCenter();
		Vector3D bbox = mesh.getBoundingBox();
		qout.output(QString().sprintf("Load mesh: %s; Size: %d", mesh.getMeshName().c_str(), mesh.vertCount()), OUT_TERMINAL);
		qout.output(QString().sprintf("Center: (%f, %f, %f)\nDimension: (%f, %f, %f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_TERMINAL);	
		
		mProcessors[obj]->init(&mesh, &mEngineWrapper);
		mRenderManagers[obj]->mesh_color = preset_mesh_colors[obj%2];
	}

	/* ---- update mesh-dependent ui ---- */
	if (mMeshCount >= 1) {
		ui.spinBox1->setMinimum(0);
		ui.spinBox1->setMaximum(mMeshes[0]->vertCount() - 1);
		ui.horizontalSlider1->setMinimum(0);
		ui.horizontalSlider1->setMaximum(mMeshes[0]->vertCount() - 1);
		ui.spinBox1->setValue(0);	
	}
	if (mMeshCount >= 2) {		
		ui.spinBox2->setMinimum(0);
		ui.spinBox2->setMaximum(mMeshes[1]->vertCount()-1);
		ui.horizontalSlider2->setMinimum(0);
		ui.horizontalSlider2->setMaximum(mMeshes[1]->vertCount()-1);
		ui.spinBox2->setValue(0);
	}
	ui.glMeshWidget->fieldView(mMeshes[0]->getCenter(), mMeshes[0]->getBoundingBox());

	mRenderManagers[0]->selected = true;
	mObjInFocus = 0;
}

void QZGeometryWindow::registerPreprocess()
{
	if (g_task != TASK_REGISTRATION || mMeshCount != 2) return;

	computeFunctionMaps(40);
	mShapeMatcher.initialize(mProcessors[0], mProcessors[1], mEngineWrapper.getEngine());
	std::string rand_data_file = g_configMgr.getConfigValue("RAND_DATA_FILE");
	mShapeMatcher.readInRandPair(rand_data_file);

	// ground truth 
	if (mMeshes[0]->getMeshName() == "march1_1_partial") {
		cout << "Ground truth available!" << endl;
		string mapFile = "./models/map1.txt";
		mShapeMatcher.loadGroundTruth(mapFile);
	}
	else if(mMeshes[0]->getMeshSize() == mMeshes[1]->getMeshSize()) {
		mShapeMatcher.autoGroundTruth();
	}

	mShapeMatcher.setRegistrationLevels(1);
	registerTest();
//	evalDistance();
}

void QZGeometryWindow::keyPressEvent( QKeyEvent *event )
{
	switch (event->key())
	{
	case Qt::Key_1:
		if (event->modifiers() & Qt::ControlModifier) {
			ui.boxObjSelect->setCurrentIndex(ui.boxObjSelect->findText("1"));
			selectObject(ui.boxObjSelect->findText("1"));
		}
		break;

	case Qt::Key_2:
		if (event->modifiers() & Qt::ControlModifier) {
			ui.boxObjSelect->setCurrentIndex(ui.boxObjSelect->findText("2"));
			selectObject(ui.boxObjSelect->findText("2"));
		}
		break;

	case Qt::Key_0:
		if (event->modifiers() & Qt::ControlModifier) {
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

	case Qt::Key_J:
		// a temporary hack
		if (mProcessors[0]->getActiveFeatures()->id == FEATURE_DEMO) {
			mProcessors[0]->setActiveFeaturesByID(FEATURE_DEMO2);
			mProcessors[1]->setActiveFeaturesByID(FEATURE_DEMO2);
			mShapeMatcher.swapMP();
		} else if (mProcessors[0]->getActiveFeatures()->id == FEATURE_DEMO2) {
			mProcessors[0]->setActiveFeaturesByID(FEATURE_DEMO);
			mProcessors[1]->setActiveFeaturesByID(FEATURE_DEMO);
			mShapeMatcher.swapMP();
		}
		
		ui.glMeshWidget->update();
		break;

	case Qt::Key_L:
		if (event->modifiers() & Qt::AltModifier)
			toggleShowMatchingLines();
		break;

	case Qt::Key_R:
		if (event->modifiers() & Qt::AltModifier) {
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
		if (mRenderManagers[0]->displayType == RenderSettings::Mesh)
			setDisplayWireframe();
		else if (mRenderManagers[0]->displayType == RenderSettings::Wireframe)
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
		mProcessors[0]->setConstrainWeight(mProcessors[0]->getConstrainWeight()*0.5);
		qout.output("Constrain weight: " + QString::number(mProcessors[0]->getConstrainWeight()));
		break;

	case Qt::Key_Equal:
		mProcessors[0]->setConstrainWeight(mProcessors[0]->getConstrainWeight()*2);
		qout.output("Constrain weight: " + QString::number(mProcessors[0]->getConstrainWeight()));
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

	for (int i = 0; i < 2; ++i) {
		mProcessors[i]->computeKernelSignatureFeatures(vTimes, SGW_KERNEL);
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
	vector<double> vTimes(timescales, timescales + nScales);

	mProcessors[0]->computeSGW(vTimes, &transferFunc1, true, &transferScalingFunc1);

	timer.stopTimer("Time for compute SGW: ");
}

void QZGeometryWindow::deformSimple()
{
	mShapeEditor.deformSimple();

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
	for (auto handle :mProcessors[0]->getHandles())
	{
		vHandle.push_back(handle.first);
		vHandlePos.push_back(handle.second);

		vector<int> vNeighbor = mProcessors[0]->getMesh_const()->getNeighborVertexIndex(handle.first, DEFAULT_DEFORM_RING);
		sFreeIdx.insert(vNeighbor.begin(), vNeighbor.end());
	}
	vFree.insert(vFree.begin(), sFreeIdx.begin(), sFreeIdx.end());

	try {
		mProcessors[0]->deform(vHandle, vHandlePos, vFree, vNewPos, Laplace);
		mMeshes[1]->setVertexCoordinates(vFree, vNewPos);
		mMeshes[1]->setVertexCoordinates(vHandle, vHandlePos);
	} catch (runtime_error* e) {
		qout.output(e->what(), OUT_MSGBOX);
	}
	
	deformType = Laplace;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::deformSGW()
{
	/*
	if (!mProcessors[0]->isSGWComputed())
		this->computeSGW();

	vector<Vector3D> vHandlePos;
	vector<int> vHandle;
	vector<Vector3D> vNewPos;
	vector<int> vFree;

	std::set<int> sFreeIdx;
	for (auto handle : mProcessors[0]->getHandles()) {
		vHandle.push_back(handle.first);
		vHandlePos.push_back(handle.second);

		vector<int> vNeighbor = mProcessors[0]->getMesh_const()->getNeighborVertexIndex(handle.first, DEFAULT_DEFORM_RING);
		sFreeIdx.insert(vNeighbor.begin(), vNeighbor.end());
	}
	vFree.insert(vFree.begin(), sFreeIdx.begin(), sFreeIdx.end());	

	try	{
		mProcessors[0]->deform(vHandle, vHandlePos, vFree, vNewPos, SGW);
		mMeshes[1]->setVertexCoordinates(vFree, vNewPos);
		mMeshes[1]->setVertexCoordinates(vHandle, vHandlePos);
	} catch (runtime_error* e) {
		qout.output(e->what(), OUT_MSGBOX);
	}
	*/

	mShapeEditor.deformSpectralWavelet();

	deformType = SGW;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::selectObject( int index )
{
	QString text = ui.boxObjSelect->itemText(index);	
	for (auto iter = mRenderManagers.begin(); iter != mRenderManagers.end(); ++iter) (*iter)->selected = false;	
	
	if (text == "1") 
	{			 
		if (mRenderManagers.size() >= 1) mRenderManagers[0]->selected = true;
		mObjInFocus = 0;
	}
	else if (text == "2")
	{
		if (mRenderManagers.size() >= 2) mRenderManagers[1]->selected = true;
		mObjInFocus = 1;
	}
	else if (text == "All")
	{
		for (auto iter = mRenderManagers.begin(); iter != mRenderManagers.end(); ++iter) (*iter)->selected = true;	 
		mObjInFocus = 0;
	}
	else if (text == "None")
	{
		mObjInFocus = -1;
	}

	qout.output("Selected object(s): " + text);
}

void QZGeometryWindow::setRefPoint1( int vn )
{
	if (mMeshCount < 1) return;
	mProcessors[0]->setRefPointIndex(vn);
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	updateReferenceMove(0);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setRefPoint2( int vn )
{
	if (mMeshCount < 2) return;

	mProcessors[1]->setRefPointIndex(vn);
	refMove.xMove = refMove.yMove = refMove.zMove = 0;
	updateReferenceMove(1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setCommonParameter( int p )
{
	mCommonParameter = p;

	if (current_operation == Compute_HKS || current_operation == Compute_HK 
		|| current_operation == Compute_MHWS || current_operation == Compute_MHW
		|| current_operation == Compute_SGWS || current_operation == Compute_SGW)
	{
		double time_scale;
		if (mCommonParameter <= PARAMETER_SLIDER_CENTER) 
			time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
		else 
			time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(mCommonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 
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

	for ( auto iter = begin(mRenderManagers); iter != end(mRenderManagers); ++iter)
	{
		(*iter)->displayType = RenderSettings::PointCloud;
		(*iter)->glPolygonMode = GL_POINT;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayWireframe()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(true);
	ui.actionDisplayMesh->setChecked(false);
	
	for ( auto iter = begin(mRenderManagers); iter != end(mRenderManagers); ++iter)
	{
		(*iter)->displayType = RenderSettings::Wireframe;
		(*iter)->glPolygonMode = GL_LINE;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayMesh()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(true);

	for ( auto iter = begin(mRenderManagers); iter != end(mRenderManagers); ++iter)
	{
		(*iter)->displayType = RenderSettings::Mesh;
		(*iter)->glPolygonMode = GL_FILL;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayExperimental()
{
	DifferentialMeshProcessor& mp = *mProcessors[0];

	vector<double> vExp;
	mp.computeExperimentalWavelet(vExp, 30); 

	mRenderManagers[0]->normalizeSignatureFrom(vExp);

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

void QZGeometryWindow::computeCurvatureMean()
{
	for (int i = 0; i < mMeshCount; ++i) {
		DifferentialMeshProcessor& mp = *mProcessors[i];
		vector<double> vCurvature;
		mp.computeCurvature(vCurvature, 0);
		auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
		qout.output(QString().sprintf("Min curvature: %d  Max curvature: %d", *mm.first, *mm.second));

		MeshFunction *mf = new MeshFunction(mp.getMesh_const()->getMeshSize());
		mf->copyValues(vCurvature);
		mf->setIDandName(SIGNATURE_MEAN_CURVATURE, "Mean Curvature");
		mp.replaceProperty(mf);			
	}

	displaySignature(SIGNATURE_MEAN_CURVATURE);
	qout.output("Show Mean curvature");
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::computeCurvatureGauss()
{
	for (int i = 0; i < mMeshCount; ++i) {		
		DifferentialMeshProcessor& mp = *mProcessors[i];
		vector<double> vCurvature;
		mp.computeCurvature(vCurvature, 1);
		auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
		qout.output(QString().sprintf("Min curvature: %d  Max curvature: %d", *mm.first, *mm.second));
			
		MeshFunction *mf = new MeshFunction(mp.getMesh_const()->getMeshSize());
		mf->copyValues(vCurvature);
		mf->setIDandName(SIGNATURE_GAUSS_CURVATURE, "Gauss Curvature");
		mp.replaceProperty(mf);			
	}

	displaySignature(SIGNATURE_GAUSS_CURVATURE);
	qout.output("Show Gauss curvature");

	updateDisplaySignatureMenu();
}

void QZGeometryWindow::displayDiffPosition()
{
	runtime_assert(mMeshes[0]->vertCount() == mMeshes[1]->vertCount());
	int size = mMeshes[0]->vertCount();
	vector<double> vDiff;
	vDiff.resize(size);

	for (int i = 0; i < mMeshes[0]->vertCount(); ++i)
	{
		vDiff[i] = (mMeshes[0]->getVertex(i)->getPosition() - mMeshes[1]->getVertex(i)->getPosition()).length() / mMeshes[0]->getAvgEdgeLength();
	}

	mRenderManagers[0]->normalizeSignatureFrom(vDiff);
	
	if (!ui.glMeshWidget->m_bShowSignature)
		toggleShowSignature();
	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::updateReferenceMove( int obj )
{
	DifferentialMeshProcessor& mp = *mProcessors[obj]; 

	double unitMove = (mp.getMesh_const()->getBoundingBox().x + mp.getMesh_const()->getBoundingBox().y + mp.getMesh_const()->getBoundingBox().z)/300.0;
	Vector3D originalPos = mp.getMesh_const()->getVertex(mp.getRefPointIndex())->getPosition();
	
	mp.setRefPointPosition(originalPos.x + unitMove * refMove.xMove,
						   originalPos.y + unitMove * refMove.yMove,
						   originalPos.z + unitMove * refMove.zMove);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::clone()
{
	if (mMeshCount != 1) {
		std::cout << "One and only one mesh should exist for cloning" << std::endl;
		return;
	}

	allocateStorage(2);
	mMeshes[1]->cloneFrom(*mMeshes[0]);
	mMeshes[1]->gatherStatistics();
	mProcessors[1]->init(mMeshes[1], &mEngineWrapper);
	mRenderManagers[1]->mesh_color = preset_mesh_colors[1];

	qout.output(QString().sprintf("Mesh %s constructed! Size: %d", mMeshes[1]->getMeshName().c_str(), mMeshes[1]->vertCount()));
	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::reconstructMHB()
{
	double ratio = min((double)mCommonParameter/PARAMETER_SLIDER_CENTER, 1.0);
	int nEig = mProcessors[0]->getMHB(MeshLaplacian::CotFormula).eigVecCount() * ratio;
	mShapeEditor.manifoldHarmonicsReconstruct(nEig);
	std::cout << "Reconstruct with " << nEig << " eigenvectors" << std::endl;

	ui.glMeshWidget->update();
}

void QZGeometryWindow::reconstructSGW()
{
	if (!mProcessors[0]->isSGWComputed())
		this->computeSGW();
	
	vector<double> vx, vy, vz;
	CStopWatch timer;
	timer.startTimer();
	try {
		mProcessors[0]->reconstructBySGW(vx, vy, vz, true);
	} catch (logic_error* e) {
		qout.output(e->what(), OUT_MSGBOX);
		return;
	}
	timer.stopTimer();
	qout.output(QString("SGW reconstruct time: ") + QString::number(timer.getElapsedTime()));

	mMeshes[1]->setVertexCoordinates(vx, vy, vz);

	{
		int debugIdx = mProcessors[0]->getHandles().begin()->first;
		qout.output("Original pos: " + std::string(mMeshes[0]->getVertex(debugIdx)->getPosition()));
		if (!mProcessors[0]->getHandles().empty())
			qout.output("Handle pos: " + std::string(mProcessors[0]->getHandles().begin()->second));
		qout.output("Deformed pos: " + std::string(mMeshes[1]->getVertex(debugIdx)->getPosition()));
	}

	double errorSum(0);
	for (int i = 0; i < mMeshes[0]->vertCount(); ++i)
	{
		errorSum += (mMeshes[0]->getVertex(i)->getPosition() - mMeshes[1]->getVertex(i)->getPosition()).length();
	}
	errorSum /= mMeshes[0]->vertCount() * mMeshes[0]->getAvgEdgeLength();
	qout.output("Average position error: " + QString::number(errorSum));

	ui.glMeshWidget->update();
}

void QZGeometryWindow::filterExperimental()
{
	vector<double> vx, vy, vz;
	mProcessors[0]->filterBySGW(vx, vy, vz);
	mMeshes[1]->setVertexCoordinates(vx, vy, vz);

	double errorSum(0);
	for (int i = 0; i < mMeshes[0]->vertCount(); ++i)
	{
		errorSum += (mMeshes[0]->getVertex(i)->getPosition() - mMeshes[1]->getVertex(i)->getPosition()).length();
	}
	errorSum /= mMeshes[0]->vertCount() * mMeshes[0]->getAvgEdgeLength();
	qout.output("Average position error: " + QString::number(errorSum));

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayNeighborVertices()
{
	int ring = (mCommonParameter > PARAMETER_SLIDER_CENTER) ? (mCommonParameter-PARAMETER_SLIDER_CENTER) : 1;

	int ref = mProcessors[0]->getRefPointIndex();
	std::vector<int> vn = mProcessors[0]->getMesh_const()->getNeighborVertexIndex(ref, ring);
//	std::vector<int> vn = vMP[0].getMesh()->getRingVertex(ref, ring);
	MeshFeatureList *mfl = new MeshFeatureList;

	for (auto iter = vn.begin(); iter != vn.end(); ++iter) {
		mfl->getFeatureVector()->push_back(new MeshFeature(*iter));
		mfl->setIDandName(FEATURE_NEIGHBORS, "Neighbors");
	}
	mProcessors[0]->addProperty(mfl);

	mProcessors[0]->setActiveFeaturesByID(FEATURE_NEIGHBORS);
	
	if (!ui.actionShowFeatures->isChecked())
		toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeEigenfunction()
{
	int select_eig = (mCommonParameter - PARAMETER_SLIDER_CENTER >= 0) ? (mCommonParameter - PARAMETER_SLIDER_CENTER + 1) : 1;

	for (int i = 0; i < mMeshCount; ++i) {
		DifferentialMeshProcessor& mp = *mProcessors[i];
		MeshFunction *mf = new MeshFunction(mp.getMesh_const()->getMeshSize());
		mf->copyValues(mp.getMHB(MeshLaplacian::CotFormula).getEigVec(select_eig).toStdVector());
		mf->setIDandName(SIGNATURE_EIG_FUNC, "Eigen_Function");
		mp.replaceProperty(mf);			
	}

	for (int i = 0; i < mMeshCount; ++i) {
		double *data = mProcessors[i]->getMHB(MeshLaplacian::SymCot).getEigVec(select_eig).c_ptr();
		int count = mProcessors[i]->getMesh()->getMeshSize();
		std::string varName = "eig" + boost::lexical_cast<std::string>(i) 
							  + "_" + boost::lexical_cast<std::string>(select_eig);
		mEngineWrapper.addVariable(data, count, 1, false, varName);
	}

	displaySignature(SIGNATURE_EIG_FUNC);
	current_operation = Compute_EIG_FUNC;
	qout.output("Show eigenfunction" + Int2String(select_eig));
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::computeHKS()
{
	double time_scale;
	if (mCommonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(mCommonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	for (int i = 0; i < mMeshCount; ++i) {
		DifferentialMeshProcessor& mp = *mProcessors[i];
		int meshSize = mp.getMesh_const()->getMeshSize();
		MeshFunction *mf = new MeshFunction(meshSize);
		mf->setIDandName(SIGNATURE_HKS, "HKS");

		const ManifoldHarmonics& mhb = mp.getMHB(MeshLaplacian::CotFormula);
		std::vector<double>& values = mf->getMeshFunction();
		Concurrency::parallel_for (0, meshSize, [&](int k) {
			values[k] = mhb.heatKernel(k, k, time_scale);
		});

		mp.replaceProperty(mf);	
	}
	
	displaySignature(SIGNATURE_HKS);
	current_operation = Compute_HKS;
	qout.output(QString().sprintf("HKS with timescale: %f", time_scale));
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::computeHK()
{
	double time_scale;
	if (mCommonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(mCommonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("Heat Kernel timescale: %f", time_scale));

	for (int i = 0; i < mMeshCount; ++i) {
		DifferentialMeshProcessor& mp = *mProcessors[i];
		int refPoint = mp.getRefPointIndex();
		mp.computeKernelDistanceSignature(time_scale, HEAT_KERNEL, refPoint);
	}

	displaySignature(SIGNATURE_HK);

	current_operation = Compute_HK;
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::repeatOperation()
{
	switch(current_operation)
	{
	case Compute_EIG_FUNC:
		computeEigenfunction();
		break;
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

	for (int i = 0; i < mMeshCount; ++i) {
		mProcessors[i]->computeKernelSignatureFeatures(vTimes, HEAT_KERNEL);
	}
	
	if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures();
}

void QZGeometryWindow::computeMHWFeatures()
{
	vector<double> vTimes;
	vTimes.push_back(10);
	vTimes.push_back(30);
	vTimes.push_back(90);
	vTimes.push_back(270);

	for (int i = 0; i < mMeshCount; ++i)
	{	
		mProcessors[i]->computeKernelSignatureFeatures(vTimes, MHW_KERNEL);
	}

	if (!ui.glMeshWidget->m_bShowFeatures)
		toggleShowFeatures();

}

void QZGeometryWindow::computeMHWS()
{
	double time_scale;
	if (mCommonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(mCommonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("MHW timescale: %f", time_scale));

	for (int i = 0; i < mMeshCount; ++i)
	{
		DifferentialMeshProcessor& mp = *mProcessors[i];
		mp.computeKernelSignature(time_scale, MHW_KERNEL);
	}

	displaySignature(SIGNATURE_MHWS);

	current_operation = Compute_MHWS;
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::computeSGWS()
{
	double time_scale;
	if (mCommonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(mCommonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("Spectral Graph Wavelet timescale: %f", time_scale));

	for (int i = 0; i < 2; ++i)
	{
		DifferentialMeshProcessor& mp = *mProcessors[i];
		mp.computeKernelSignature(time_scale, SGW_KERNEL);
	}

	displaySignature(SIGNATURE_SGWS);

	current_operation = Compute_SGWS;
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::displaySignature( int signatureID )
{
	for (int i = 0; i < mMeshCount; ++i)
	{
		DifferentialMeshProcessor& mp = *mProcessors[i];
		MeshProperty* vs = mp.retrievePropertyByID(signatureID);
		if (vs != NULL)
		{
			mRenderManagers[i]->normalizeSignatureFrom(dynamic_cast<MeshFunction*>(vs)->getMeshFunction_const());
//			vRS[i].logNormalizeSignatureFrom(dynamic_cast<MeshFunction*>(vs)->getMeshFunction_const());
			qout.output(QString().sprintf("Sig Min: %f; Sig Max: %f", mRenderManagers[i]->sigMin, mRenderManagers[i]->sigMax), OUT_CONSOLE);
		}
	}

	if (!ui.glMeshWidget->m_bShowSignature) toggleShowSignature();	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeMHW()
{
	double time_scale;
	if (mCommonParameter <= PARAMETER_SLIDER_CENTER) 
		time_scale = std::exp(std::log(DEFUALT_HK_TIMESCALE / MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)PARAMETER_SLIDER_CENTER) + std::log(MIN_HK_TIMESCALE));
	else 
		time_scale = std::exp(std::log(MAX_HK_TIMESCALE / DEFUALT_HK_TIMESCALE) * ((double)(mCommonParameter-PARAMETER_SLIDER_CENTER) / (double)PARAMETER_SLIDER_CENTER) + std::log(DEFUALT_HK_TIMESCALE)); 

	qout.output(QString().sprintf("MHW timescale: %f", time_scale));

	for (int i = 0; i < mMeshCount; ++i)
	{
		DifferentialMeshProcessor& mp = *mProcessors[i];
		int refPoint = mp.getRefPointIndex();
		mp.computeKernelDistanceSignature(time_scale, MHW_KERNEL, refPoint);
	}

	displaySignature(SIGNATURE_MHW);

	current_operation = Compute_MHW;
	updateDisplaySignatureMenu();
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
	mShapeMatcher.constructPyramid(nLevel, ratio, ostr);
	qout.output("Mesh hierarchy constructed!");

	ostr.close();
}

void QZGeometryWindow::detectFeatures()
{
	double featureDetectionBaseTimescale = g_configMgr.getConfigValueDouble("FEATURE_DETECTION_BASE_TIMESCALE");
	double featureDetectionTMultiplier = g_configMgr.getConfigValueDouble("FEATURE_DETECTION_T_MULTIPLIER");
	int numDetectScales = g_configMgr.getConfigValueInt("FEATURE_DETECTION_NUM_SCALES");
	double featureDetectionExtremaThresh = g_configMgr.getConfigValueDouble("FEATURE_DETECTION_EXTREMA_THRESH");
	int detectRing = g_configMgr.getConfigValueInt("FEATURE_DETECTION_RING");

	qout.output("-- Detect initial features --");
	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
		mShapeMatcher.detectFeatures(obj, detectRing, numDetectScales, featureDetectionBaseTimescale, 
									featureDetectionTMultiplier, featureDetectionExtremaThresh);
	});
	qout.output("Multi-scale mesh features detected!");
	std::cout << "Mesh1 features #: " << mShapeMatcher.getSparseFeatures(0).size() 
		      << "; Mesh 2 features #: " << mShapeMatcher.getSparseFeatures(1).size() << std::endl;

#if 0	
	if (mShapeMatcher.hasGroundTruth()) 
	{
		vector<HKSFeature>& vf1 = mShapeMatcher.getSparseFeatures(0), &vf2 = mShapeMatcher.getSparseFeatures(1);
		std::map<int, int> feature_count;
		int count_strict = 0; 
		int count_one_neighbor = 0;

		int count_possible2 = 0;

		double sim1(0.);

		for (auto iter1 = vf1.begin(); iter1 != vf1.end(); ) {
			bool candFound = false;
			for (auto iter2 = vf2.begin(); iter2 != vf2.end(); ++iter2) {
// 				double sim = shapeMatcher.calPointHksDissimilarity(&vMP[0], &vMP[1], iter1->m_index, iter2->m_index, vTimes, 1);
// 				if (sim < 0.20) {candFound = true; break;}

				if (mMeshes[1]->isInNeighborRing(iter2->m_index, iter1->m_index, 2)) {
					candFound = true; break;
				}
			}
			if (!candFound) {
				iter1 = vf1.erase(iter1);
				continue;
			}
			else ++iter1;
		}

		for (auto iter2 = vf2.begin(); iter2 != vf2.end(); ) {
			bool candFound = false;
			for (auto iter1 = vf1.begin(); iter1 != vf1.end(); ++iter1) {
// 				double sim = shapeMatcher.calPointHksDissimilarity(&vMP[0], &vMP[1], iter1->m_index, iter2->m_index, vTimes, 1);
// 				if (sim < 0.20) {candFound = true; break;}
				if (mMeshes[0]->isInNeighborRing(iter1->m_index, iter2->m_index, 2)) {
					candFound = true; break;
				}
			}
			if (!candFound) {
				iter2 = vf2.erase(iter2);
				continue;
			}
			else ++iter2;
		}

		for (auto iter1 = vf1.begin(); iter1 != vf1.end(); ++iter1) {
			for (auto iter2 = vf2.begin(); iter2 != vf2.end(); ++iter2) {
				if (mMeshes[1]->isInNeighborRing(iter1->m_index, iter2->m_index, 2)) {
					count_one_neighbor++;
					if (feature_count.find(iter1->m_index) == feature_count.end())
						feature_count.insert(make_pair(iter1->m_index, 1));
					else feature_count[iter1->m_index] += 1;

					if (iter1->m_index == iter2->m_index) {
						count_strict++;
					}

// 					double sim = shapeMatcher.calPointHksDissimilarity(&vMP[0], &vMP[1], iter1->m_index, iter2->m_index, vTimes, 1);
// 					cout << "corresponded sim: " << sim << endl;
// 					sim1 += sim;
				}
			}
		}
		cout << "v1: " << vf1.size() << "  v2: " << vf2.size() << endl;
		cout << "-- Potential matches: " << count_strict << "/" << count_one_neighbor << endl;
//		cout << "Average similarity of matched: " << sim1 / double(count_possible1) << endl;
	}
#endif	

	if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::matchFeatures()
{
	int force_matching = g_configMgr.getConfigValueInt("FORCE_MATCHING");
	string string_override = "";
	bool alreadyMached = false;

#if 0	
	if (force_matching == 1)
	{
		if (mMeshes[0]->getMeshName() == "horse0")
		{
			string_override = g_configMgr.getConfigValue("HORSE0_FEATURE_OVERRIDE");
		}
		else if (mMeshes[0]->getMeshName() == "eight")
		{
			string_override = g_configMgr.getConfigValue("EIGHT_FEATURE_OVERRIDE");
		}
		if (string_override != "")
		{
			vector<int> idx_override = ZUtil::splitStringToInt(string_override);
			if (!idx_override.empty())
			{
				vector<MatchPair> vmp;
				for (auto iter = begin(idx_override); iter != end(idx_override); ++iter)
				{
					vmp.push_back(MatchPair(*iter, *iter));
				}
				mShapeMatcher.forceInitialAnchors(vmp);
				qout.output("!!Matched anchors manually assigned!!");
				alreadyMached = true;
			}
		}
	}
	else if (force_matching == 2)
	{
		if (mMeshes[0]->getMeshName() == "horse0")
		{
			string_override = g_configMgr.getConfigValue("MATCHING_HORSE_FILE");
		}
		else if (mMeshes[0]->getMeshName() == "eight")
		{
			string_override = g_configMgr.getConfigValue("MATCHING_EIGHT_FILE");
		}
		if (string_override != "" && mShapeMatcher.loadInitialFeaturePairs(string_override))
		{
			qout.output("!!Matched anchors manually assigned!!");
			alreadyMached = true;
		}
	}
#endif
	
	int matching_method = g_configMgr.getConfigValueInt("FEATURE_MATCHING_METHOD");
	std::string log_filename = g_configMgr.getConfigValue("MATCH_OUTPUT_FILE");

	qout.output("-- Match initial features --");
	ofstream ofstr(log_filename.c_str(), ios::trunc);
	CStopWatch timer;
	timer.startTimer();
	if (matching_method == 2)	//tensor
	{
		std::cout << "Now do tensor matching!" << endl;
		double matching_thresh_2 = g_configMgr.getConfigValueDouble("MATCHING_THRESH_2");
		double tensor_matching_timescasle = g_configMgr.getConfigValueDouble("TENSOR_MATCHING_TIMESCALE");

		const vector<HKSFeature>& vftFine1 = mShapeMatcher.getSparseFeatures(0);
		const vector<HKSFeature>& vftFine2 = mShapeMatcher.getSparseFeatures(1);
		vector<int> vFeatures1, vFeatures2;
		for_each(vftFine1.begin(), vftFine1.end(), [&](const HKSFeature& f){vFeatures1.push_back(f.m_index);});
		for_each(vftFine2.begin(), vftFine2.end(), [&](const HKSFeature& f){vFeatures2.push_back(f.m_index);});

		vector<MatchPair> vPairs;
		double vPara[] = {40, 0.8, 400};
		double matchScore;
		matchScore = ShapeMatcher::TensorGraphMatching6(mEngineWrapper.getEngine(), mProcessors[0], mProcessors[1], vftFine1, vftFine2, vPairs, tensor_matching_timescasle, matching_thresh_2, /*verbose=*/true);
		//matchScore = DiffusionShapeMatcher::TensorMatchingExt(m_ep, &vMP[0], &vMP[1], vFeatures1, vFeatures2, vPairs, 0, vPara, cout, true);

		if (1 == g_configMgr.getConfigValueInt("GROUND_TRUTH_AVAILABLE"))
		{
			vector<double> vTimes;
			vTimes.push_back(20); vTimes.push_back(40); vTimes.push_back(80); vTimes.push_back(160); vTimes.push_back(320);
			for (auto iter = vPairs.begin(); iter != vPairs.end(); )
			{
				if (!mMeshes[1]->isInNeighborRing(iter->m_idx1, iter->m_idx2, 2))
					iter->m_note = -1;

				double dissim = mShapeMatcher.calPointHksDissimilarity(mProcessors[0], mProcessors[1], iter->m_idx1, iter->m_idx2, vTimes, 1);
				if (dissim > 0.18)
				{
					iter = vPairs.erase(iter);
					continue;
				}
				else ++iter;
			}
		}			

		std::cout << "Tensor match score: " << matchScore << endl;
		mShapeMatcher.forceInitialAnchors(vPairs);

		//shapeMatcher.matchFeaturesTensor(ofstr, tensor_matching_timescasle, matching_thresh_2);
	}
	else if (matching_method == 1)	// traditional pair-based matching
	{
		double matching_thresh_1 = g_configMgr.getConfigValueDouble("MATCHING_THRESH_1");
		mShapeMatcher.matchFeatures(ofstr, matching_thresh_1);
	}
	else // point based
	{
		mShapeMatcher.matchFeatureSimple();
	}
	timer.stopTimer();
	ofstr.close();
	qout.output(QString().sprintf("Initial features matched! Matched#:%d. Time elapsed:%f", mShapeMatcher.getMatchedFeaturesResults(mShapeMatcher.getAlreadyMatchedLevel()).size(), timer.getElapsedTime()));


	const std::vector<MatchPair>& result = mShapeMatcher.getMatchedFeaturesResults(mShapeMatcher.getAlreadyMatchedLevel());

	mShapeMatcher.evaluateWithGroundTruth(result);

	if (!ui.glMeshWidget->m_bDrawMatching) toggleDrawMatching();
}

void QZGeometryWindow::registerStep()
{
	string log_filename = g_configMgr.getConfigValue("REGISTER_OUTPUT_FILE");

	qout.output(QString().sprintf("-- Register level %d --", mShapeMatcher.getAlreadyRegisteredLevel() - 1));
	
	ofstream ofstr;	
	if (mShapeMatcher.getAlreadyRegisteredLevel() == mShapeMatcher.getTotalRegistrationLevels())
		ofstr.open(log_filename.c_str(), ios::trunc);
	else
		ofstr.open(log_filename.c_str(), ios::app);

	int regMethod = g_configMgr.getConfigValueInt("REGISTRATION_METHOD");

	double time_elapsed = time_call([&]{
		if (regMethod == 1) mShapeMatcher.refineRegister(ofstr);
		else if (regMethod == 2) mShapeMatcher.refineRegister2(ofstr);
	}) / 1000.0;

	int level = mShapeMatcher.getAlreadyRegisteredLevel();
	const vector<MatchPair>& vf = mShapeMatcher.getMatchedFeaturesResults(mShapeMatcher.getAlreadyMatchedLevel());
	const vector<MatchPair>& vr = mShapeMatcher.getRegistrationResults(mShapeMatcher.getAlreadyRegisteredLevel());

	qout.output(QString().sprintf("Registration level %d finished! Time elapsed:%f\n-Features Matched:%d; Registered:%d",
		level, time_elapsed, vf.size(), vr.size()));
	qout.output(QString().sprintf("Registered ratio: %f (%d/%d)", 
		double(vr.size())/mShapeMatcher.getMesh(0, level)->vertCount(),
		vr.size(), mShapeMatcher.getMesh(0, level)->vertCount()));
	/* ---- evaluation ---- */
	if (mShapeMatcher.hasGroundTruth())
	{
		cout << "Features - ";
		mShapeMatcher.evaluateWithGroundTruth(vf);
		cout << "Registrations - ";
		mShapeMatcher.evaluateWithGroundTruth(vr);
	}

// 	double vError[3];
// 	for(int k = 0; k < 2; ++k)
// 	{ 
// 		vError[k] = DiffusionShapeMatcher::evaluateDistortion(vr, &m_mesh[0], &m_mesh[1], shapeMatcher.m_randPairs, 200 * k);
// 	}
// 
// 	qout.output(QString().sprintf("Registration error: %f, %f", vError[0], vError[1]));

	if (!ui.glMeshWidget->m_bDrawRegistration)
		toggleDrawRegistration();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::registerFull()
{
	while(mShapeMatcher.getAlreadyMatchedLevel() > 0)
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
	if (ui.glMeshWidget->m_nMeshLevel < mShapeMatcher.getPyramidLevels()-1)
		ui.glMeshWidget->m_nMeshLevel++;

	qout.output("Display mesh level " + QString::number(ui.glMeshWidget->m_nMeshLevel));
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeBiharmonic()
{
	for (int i = 0; i < mMeshCount; ++i)
	{
		DifferentialMeshProcessor& mp = *mProcessors[i];
		int refPoint = mp.getRefPointIndex();
		mp.computeBiharmonicDistanceSignature(refPoint);
	}

	displaySignature(SIGNATURE_BIHARMONIC_DISTANCE);
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::evalDistance()
{
	if (mMeshCount < 2) return;

	CStopWatch timer;
	timer.startTimer();
	Concurrency::parallel_invoke(
		[&](){ cout << "Error geodesic1: " << ShapeMatcher::evaluateDistance(*mProcessors[0], *mProcessors[1], DISTANCE_GEODESIC, std::vector<double>(), mShapeMatcher.m_randPairs, 0) << endl; },
//		[&](){ cout << "Error geodesic2: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl; },
		[&](){ cout << "Error biharmonic1: " << ShapeMatcher::evaluateDistance(*mProcessors[0], *mProcessors[1], DISTANCE_BIHARMONIC, std::vector<double>(), mShapeMatcher.m_randPairs, 0) << endl; },
		[&](){ cout << "Error biharmonic2: " << ShapeMatcher::evaluateDistance(*mProcessors[0], *mProcessors[1], DISTANCE_BIHARMONIC, std::vector<double>(), mShapeMatcher.m_randPairs, 500) << endl; }
//		[&](){cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 30.), shapeMatcher.m_randPairs, 0) << endl;},
//		cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 90.), shapeMatcher.m_randPairs, 0) << endl;
//		cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 270.), shapeMatcher.m_randPairs, 0) << endl;
	);
//	cout << "Error geodesic1: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl;
	timer.stopTimer();
	cout << "Eval Dist time (ppl): " << timer.getElapsedTime() << endl;

}

void QZGeometryWindow::decomposeSingleLaplacian( int obj, int nEigVec, MeshLaplacian::LaplacianType laplacianType /*= CotFormula*/ )
{
	DifferentialMeshProcessor& mp = *mProcessors[obj];
	const CMesh& mesh = *mMeshes[obj];
	assert(nEigVec < mesh.getMeshSize());

	if (!mp.getMHB(laplacianType).empty()) return;
	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "output/" + mp.getMesh_const()->getMeshName() + ".mhb." + s_idx;
	
	if (LOAD_MHB_CACHE && ZUtil::fileExist(pathMHB))	// MHB cache available for the current mesh
	{
		ifstream ifs(pathMHB.c_str());
		mp.loadMHB(pathMHB, laplacianType);
		ifs.close();
	}
	else // need to compute Laplacian and to cache
	{
		mp.decomposeLaplacian(nEigVec, laplacianType);
		mp.saveMHB(pathMHB, laplacianType);
		qout.output("MHB saved to " + pathMHB);
	}

	std::cout << "Min EigVal: " << mp.getMHB(MeshLaplacian::CotFormula).getEigVals().front() 
			  << "; Max EigVal: " << mp.getMHB(MeshLaplacian::CotFormula).getEigVals().back() << std::endl;
}

void QZGeometryWindow::decomposeLaplacians( MeshLaplacian::LaplacianType laplacianType /*= CotFormula*/ )
{
	int totalToDecompose = 0;
	int nEigVec = DEFAULT_EIGEN_SIZE;

	for (int obj = 0; obj < mMeshCount; ++obj) {
		if (!mProcessors[obj]->hasLaplacian(laplacianType))
			throw std::logic_error("Laplacian type not valid!");
		if (laplacianRequireDecompose(obj, nEigVec, laplacianType)) 
			++totalToDecompose;
	}
	std::cout << totalToDecompose << " mesh Laplacians require explicit decomposition" << std::endl;
		
	if (totalToDecompose <= 1) {        
		Concurrency::parallel_for(0, mMeshCount, [&](int obj){
			decomposeSingleLaplacian(obj, nEigVec, laplacianType);    
		});
	} else {    // if both need explicit decomposition, then must run in sequence in Matlab
		for(int obj = 0; obj < mMeshCount; ++obj) {
			decomposeSingleLaplacian(obj, nEigVec, laplacianType);    
		}
	}

	for (int l = 0; l < MeshLaplacian::LaplacianTypeCount; ++l)
		m_actionComputeLaplacians[l]->setChecked(false);
	m_actionComputeLaplacians[laplacianType]->setChecked(true);
}

void QZGeometryWindow::saveSignature()
{
	if (mRenderManagers[0]->vOriginalSignature.empty())
	{
		qout.output("No signature available", OUT_MSGBOX);
		return;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Signature to File"),
		"./output/signature.txt",
		tr("Text Files (*.txt *.dat)"));
	vector<double> vSig = mRenderManagers[0]->vOriginalSignature;

	ZUtil::vector2file<double>(fileName.toStdString(), vSig);
}

void QZGeometryWindow::addMesh()
{
	QStringList filenames =  QFileDialog::getOpenFileNames(this, "Select one or more mesh files to open",
														   "../../Data/", "Meshes (*.obj *.off *.ply)");
	int cur_obj = mMeshCount;
	allocateStorage(++mMeshCount);

	CStopWatch timer;
	timer.startTimer();
	CMesh& mesh = *mMeshes[cur_obj];
	mesh.Load(filenames.begin()->toStdString());
	mesh.scaleAreaToVertexNum();
	mesh.gatherStatistics();
	timer.stopTimer();
	cout << "Time to load meshes: " << timer.getElapsedTime() << "s" << endl;

	Vector3D center = mesh.getCenter(), bbox = mesh.getBoundingBox();
	qout.output(QString().sprintf("Load mesh: %s; Size: %d", mesh.getMeshName().c_str(), mesh.vertCount()), OUT_CONSOLE);
	qout.output(QString().sprintf("Center: (%f,%f,%f)\nDimension: (%f,%f,%f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_CONSOLE);

	mProcessors[cur_obj]->init(&mesh, &mEngineWrapper);

	mRenderManagers[cur_obj]->selected = true;
	mRenderManagers[cur_obj]->mesh_color = preset_mesh_colors[cur_obj%2];

	if (cur_obj == 0)
	{
		ui.glMeshWidget->fieldView(mMeshes[0]->getCenter(), mMeshes[0]->getBoundingBox());
		ui.spinBox1->setMinimum(0);
		ui.spinBox1->setMaximum(mMeshes[0]->vertCount() - 1);
		ui.horizontalSlider1->setMinimum(0);
		ui.horizontalSlider1->setMaximum(mMeshes[0]->vertCount() - 1);
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::updateDisplaySignatureMenu()
{
	if (mObjInFocus < 0) return;

	const std::vector<MeshProperty*> vProperties = mProcessors[mObjInFocus]->properties();
	vector<MeshFunction*> vSigFunctions;
	for (MeshProperty* pp : vProperties) {
		if (pp->id > SIGNATURE_ID && pp->id < SIGNATURE_ID_COUNT)
			vSigFunctions.push_back(dynamic_cast<MeshFunction*>(pp));
	}
	
	QList<QAction*> signatureActions = ui.menuSignature->actions();
	for (QAction* qa : m_actionDisplaySignatures) {
		if (vSigFunctions.end() == find_if(vSigFunctions.begin(), vSigFunctions.end(), [&](MeshFunction* mf){ return mf->name == qa->text().toStdString();}))
		{
			ui.menuSignature->removeAction(qa);
			delete qa;			
		}
	}

	for (MeshFunction* pmf : vSigFunctions) {
		if (m_actionDisplaySignatures.end() == find_if(m_actionDisplaySignatures.begin(), m_actionDisplaySignatures.end(), [&](QAction* pa){ return pa->text().toStdString() == pmf->name;}))
		{
			QAction* newDisplayAction = new QAction(pmf->name.c_str(), this);
			m_actionDisplaySignatures.push_back(newDisplayAction);
			ui.menuSignature->addAction(m_actionDisplaySignatures.back());
			signatureSignalMapper->setMapping(newDisplayAction, (int)pmf->id);
			QObject::connect(newDisplayAction, SIGNAL(triggered()), signatureSignalMapper, SLOT(map()));
		}	
	}
}

void QZGeometryWindow::computeSimilarityMap( int simType )
{
	switch(simType)
	{
	case SIM_TYPE_1:
		Concurrency::parallel_for(0, mMeshCount, [&](int obj){
			mProcessors[obj]->computeSimilarityMap1(mProcessors[obj]->getRefPointIndex());
		});
		break;

	case SIM_TYPE_2:
		Concurrency::parallel_for(0, mMeshCount, [&](int obj){
			mProcessors[obj]->computeSimilarityMap2(mProcessors[obj]->getRefPointIndex());
		});
		break;

	case SIM_TYPE_3:
		Concurrency::parallel_for(0, mMeshCount, [&](int obj){
			mProcessors[obj]->computeSimilarityMap3(mProcessors[obj]->getRefPointIndex());
		});
		break;
	}

	displaySignature(SIGNATURE_SIMILARITY_MAP);
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::constructLaplacians( MeshLaplacian::LaplacianType laplacianType )
{
	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
		mProcessors[obj]->constructLaplacian(laplacianType);
	});
}

void QZGeometryWindow::computeLaplacian( int lapType )
{
	MeshLaplacian::LaplacianType laplacianType = (MeshLaplacian::LaplacianType)lapType;
	constructLaplacians(laplacianType);
	decomposeLaplacians(laplacianType);
}

void QZGeometryWindow::saveMatchingResult()
{
	if (g_task != TASK_REGISTRATION) return;

	const vector<MatchPair>& vPairs = mShapeMatcher.getInitialMatchedFeaturePairs();
	if (vPairs.empty()) {
		qout.output("No matching result available", OUT_MSGBOX);
		return;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Matching Result to File"),
													"./output/matching.txt",
													tr("Text Files (*.txt *.dat)"));
	
	ofstream ofs(fileName.toStdString().c_str(), ios::trunc);

	ofs << vPairs.size() << endl;
	for (auto iter = vPairs.begin(); iter != vPairs.end(); ++iter) {
		ofs << iter->m_idx1 << ' ' << iter->m_idx2 << ' ' << iter->m_score << endl;
	}

	ofs.close();
}

void QZGeometryWindow::loadMatchingResult()
{
	if (g_task != TASK_REGISTRATION) return;

	QString filename =  QFileDialog::getOpenFileName(this, "Select one or more mesh files to open",
													 "./output/", tr("Text Files (*.txt *.dat)"));

	mShapeMatcher.loadInitialFeaturePairs(filename.toStdString());

	if (ui.glMeshWidget->m_bDrawMatching)
		ui.glMeshWidget->update();
	else toggleDrawMatching();
}

void QZGeometryWindow::registerTest()
{
	//shapeMatcher.registerTesting1();
	//shapeMatcher.regsiterTesting2();
	//shapeMatcher.dataTesting1();
	//shapeMatcher.sparseMatchingTesting();
	//shapeMatcher.localCorrespondenceTesting();

// 	shapeMatcher.generateExampleMatching(20);
// 	if (!ui.glMeshWidget->m_bDrawMatching)
// 		toggleDrawMatching();

	ui.glMeshWidget->update();
}

bool QZGeometryWindow::laplacianRequireDecompose( int obj, int nEigVec, MeshLaplacian::LaplacianType laplacianType ) const
{
	const DifferentialMeshProcessor& mp = *mProcessors[obj];
	const CMesh& mesh = *mMeshes[obj];
	
	if (!mp.getMHB(laplacianType).empty()) return false; // already decomposed     
	if (!LOAD_MHB_CACHE) return true;    

	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "output/" + mp.getMesh_const()->getMeshName() + ".mhb." + s_idx;

	if (!ZUtil::fileExist(pathMHB)) return true;

	ifstream ifs(pathMHB.c_str(), ios::binary);
	int nEig, nSize;
	ifs.read((char*)&nEig, sizeof(int));
	ifs.read((char*)&nSize, sizeof(int));
	ifs.close();

	if (nEig != nEigVec || nSize != mesh.getMeshSize()) return true;

	return false;
}

void QZGeometryWindow::allocateStorage( int newMeshCount )
{
	int existingMeshCount = mMeshes.size();

	assert(newMeshCount > existingMeshCount);

	for (int k = 0; k < newMeshCount - existingMeshCount; ++k) {
		mMeshes.push_back(new CMesh());
		mProcessors.push_back(new DifferentialMeshProcessor());
		mRenderManagers.push_back(new RenderSettings());
	}

	mMeshCount = newMeshCount;
}

void QZGeometryWindow::computeFunctionMaps( int num )
{
	ZGeom::DenseMatrixd funcMap1(num, num), funcMap2(num, num);
	const MeshLaplacian &lap1 = mProcessors[0]->getMeshLaplacian(MeshLaplacian::CotFormula);
	const MeshLaplacian &lap2 = mProcessors[1]->getMeshLaplacian(MeshLaplacian::CotFormula);
	const ManifoldHarmonics& mhb1 = mProcessors[0]->getMHB(MeshLaplacian::CotFormula);
	const ManifoldHarmonics& mhb2 = mProcessors[1]->getMHB(MeshLaplacian::CotFormula);
	ZGeom::SparseMatrixCSR<double, int> csrMat1, csrMat2;
	lap1.getW().convertToCSR(csrMat1, ZGeom::MAT_FULL);
	lap2.getW().convertToCSR(csrMat2, ZGeom::MAT_FULL);

	for (int i = 0; i < num; ++i) {
		const ZGeom::VecNd& eig1i = mhb1.getEigVec(i);
		const ZGeom::VecNd& eig2i = mhb2.getEigVec(i);
		for (int j = 0; j < num; ++j) {
			const ZGeom::VecNd& eig1j = mhb1.getEigVec(j);
			const ZGeom::VecNd& eig2j = mhb2.getEigVec(j);
			funcMap1(i,j) = ZGeom::innerProductSym(eig2i, csrMat1, eig1j);
			funcMap2(i,j) = ZGeom::innerProductSym(eig1i, csrMat2, eig2j);
		}
	}
	
	funcMap1.print("output/funcmap1.txt");
	funcMap2.print("output/funcmap2.txt");

	lap1.getLS().print("output/Ls1.txt");
	lap1.getW().print("output/W1.txt");
	lap2.getLS().print("output/Ls2.txt");
	lap2.getW().print("output/W2.txt");

	mhb1.printEigVals("output/eigvals1.txt");
	mhb2.printEigVals("output/eigvals2.txt");
}

void QZGeometryWindow::verifyAreas() const
{
	for (int obj = 0; obj < 2; ++obj) {
		double areaSum(0);
		for (int i = 0; i < mMeshes[obj]->faceCount(); ++i) {
			areaSum += mMeshes[obj]->calFaceArea(i);
		}
		double weightSum(0);
		const MeshLaplacian& laplacian = mProcessors[obj]->getMeshLaplacian(MeshLaplacian::CotFormula);
		for (int i = 0; i < mMeshes[obj]->vertCount(); ++i) {
			weightSum += laplacian.getW().getElemByIndex(i);
		}
		std::cout << "Vert count: " << mMeshes[obj]->vertCount() << std::endl;
		std::cout << "Total surface area: " << areaSum << std::endl;
		std::cout << "Total vert weight: " << weightSum << std::endl;
	}
}

void QZGeometryWindow::revert()
{
	mShapeEditor.revert();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::deformDifferential()
{
	mShapeEditor.deformDifferential();
	ui.glMeshWidget->update();
}
