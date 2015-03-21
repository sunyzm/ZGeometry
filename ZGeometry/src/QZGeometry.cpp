#include "QZGeometry.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <set>
#include <unordered_set>
#include <random>
#include <stdexcept>
#include <ppl.h>
#include <boost/lexical_cast.hpp>
#include <Shellapi.h>
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>
#include <QTime>
#include <QImage>
#include <ZGeom/util.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/DenseMatrix.h>
#include "heat_diffusion.h"
#include "hole_fairing.h"


using std::vector;
using ZGeom::Colorf;
using ZGeom::logic_assert;
using ZGeom::runtime_assert;
using ZGeom::MatlabEngineWrapper;
using ZGeom::ColorSignature;
using ZGeom::MeshRegion;

QZGeometryWindow::QZGeometryWindow(QWidget *parent,  Qt::WindowFlags flags) : QMainWindow(parent, flags)
{
	mMeshCount				= 0;
	mObjInFocus				= -1;
	mCommonParameter		= gSettings.PARAMETER_SLIDER_CENTER;
	mLastOperation			= None;
	mDeformType				= DEFORM_Simple;
    mSignatureMode          = ZGeom::SM_Normalized;
	mActiveLalacian			= Umbrella;
	mColorMapType			= ZGeom::CM_JET;
	mDiffMax				= 2.0;
	mCurrentBasisScale		= 0;

	/* setup ui and connections */
	ui.setupUi(this);
	this->makeConnections();
	
	// comboBoxSigMode
	ui.comboBoxSigMode->clear();
    for (int i = 0; i < ZGeom::SM_CountSigModes; ++i)
		ui.comboBoxSigMode->addItem(QString(StrSignatureModeNames[i].c_str()));

	// comboBoxLaplacian
	ui.comboBoxLaplacian->clear();
	for (int i = 0; i < 4; ++i)
		ui.comboBoxLaplacian->addItem(StrLaplacianTypes[i].c_str());
	ui.comboBoxLaplacian->setCurrentIndex(ui.comboBoxLaplacian->findText("CotFormula"));

	// toolbar
	ui.spinBoxParameter->setMinimum(0);
	ui.spinBoxParameter->setMaximum(2 * gSettings.PARAMETER_SLIDER_CENTER);
	ui.horizontalSliderParamter->setMinimum(0);
	ui.horizontalSliderParamter->setMaximum(2 * gSettings.PARAMETER_SLIDER_CENTER);
	ui.horizontalSliderParamter->setSliderPosition(gSettings.PARAMETER_SLIDER_CENTER);

	// status bar
	mStatusLabel.setParent(ui.statusBar);
	ui.statusBar->addPermanentWidget(&mStatusLabel);
	qout.setLabel(&mStatusLabel);
	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);    	
	mStatusLabel.setText("Editing");

	setDisplayMesh();
	setEditModeMove();
}

QZGeometryWindow::~QZGeometryWindow()
{
	for (QAction* a : m_actionDisplaySignatures) delete a;
	for (QAction* a : m_actionComputeLaplacians) delete a;
	for (QAction* a : m_actionDisplayFeatures) delete a;
    for (QAction* a : m_actionDisplayLines) delete a;

	delete m_laplacianSignalMapper;
	delete m_signatureSignalMapper;
	delete m_featureSignalMapper;
    delete m_linesSignalMapper;
}

void QZGeometryWindow::makeConnections()
{
	QObject::connect(ui.actionAboutQt, SIGNAL(triggered()), this, SIGNAL(displayQtVersion()));

	/*  actionComputeLaplacians  */
	int laplacianTypeCount = LaplacianTypeCount;
	m_actionComputeLaplacians.resize(laplacianTypeCount);
	m_laplacianSignalMapper = new QSignalMapper(this);
	for (int t = 0; t < laplacianTypeCount; ++t) {
		m_actionComputeLaplacians[t] = new QAction(QString("Laplacian type ") + QString::number(t), this);
		ui.menuComputeLaplacian->addAction(m_actionComputeLaplacians[t]);
		m_laplacianSignalMapper->setMapping(m_actionComputeLaplacians[t], t);
		QObject::connect(m_actionComputeLaplacians[t], SIGNAL(triggered()), m_laplacianSignalMapper, SLOT(map()));
	}
	QObject::connect(m_laplacianSignalMapper, SIGNAL(mapped(int)), this, SLOT(computeLaplacian(int)));

	/*  actionDisplaySignatures  */
	m_signatureSignalMapper = new QSignalMapper(this);
	QObject::connect(m_signatureSignalMapper, SIGNAL(mapped(QString)), this, SLOT(displaySignature(QString)));

	/* actionDisplayFeatures */
	m_featureSignalMapper = new QSignalMapper(this);
	QObject::connect(m_featureSignalMapper, SIGNAL(mapped(QString)), this, SLOT(displayFeature(QString)));

    /* actionDisplayLines */
    m_linesSignalMapper = new QSignalMapper(this);
    QObject::connect(m_linesSignalMapper, SIGNAL(mapped(QString)), this, SLOT(displayLine(QString)));

	////////    Toolbar Controls    ////////
	QObject::connect(ui.spinBox1, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked1(int)));
	QObject::connect(ui.horizontalSlider1, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked1(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked1(int)), this, SLOT(setMesh1RefPoint(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked1(int)), ui.horizontalSlider1, SLOT(setValue(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked1(int)), ui.spinBox1, SLOT(setValue(int)));

	QObject::connect(ui.spinBox2, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked2(int)));
	QObject::connect(ui.horizontalSlider2, SIGNAL(valueChanged(int)), ui.glMeshWidget, SIGNAL(vertexPicked2(int)));
	QObject::connect(ui.glMeshWidget, SIGNAL(vertexPicked2(int)), this, SLOT(setMesh2RefPoint(int)));
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

	////////    Tabbed Controls    ////////
	QObject::connect(ui.sliderApprox1, SIGNAL(valueChanged(int)), this, SLOT(continuousApprox1(int)));
	QObject::connect(ui.sliderApprox2, SIGNAL(valueChanged(int)), this, SLOT(continuousApprox2(int)));
	QObject::connect(ui.sliderApprox3, SIGNAL(valueChanged(int)), this, SLOT(continuousApprox3(int)));
	QObject::connect(ui.sliderApprox4, SIGNAL(valueChanged(int)), this, SLOT(continuousApprox4(int)));
	QObject::connect(ui.sliderPointSize, SIGNAL(valueChanged(int)), this, SLOT(setFeaturePointSize(int)));
	QObject::connect(ui.sliderEditBasis, SIGNAL(valueChanged(int)), this, SLOT(displayBasis(int)));
	QObject::connect(ui.comboBoxSigMode, SIGNAL(activated(const QString&)), this, SLOT(setSignatureMode(const QString&)));
	QObject::connect(ui.sliderSigMin, SIGNAL(valueChanged(int)), this, SLOT(updateSignatureMin(int)));
	QObject::connect(ui.sliderSigMax, SIGNAL(valueChanged(int)), this, SLOT(updateSignatureMax(int)));
	QObject::connect(ui.comboBoxLaplacian, SIGNAL(activated(const QString&)), this, SLOT(setLaplacianType(const QString&)));

	////////    Menus	////////
	////  File  ////
	QObject::connect(ui.actionExit, SIGNAL(triggered()), this, SLOT(close()));
	QObject::connect(ui.actionSaveSignature, SIGNAL(triggered()), this, SLOT(saveSignature()));
	QObject::connect(ui.actionAddMesh, SIGNAL(triggered()), this, SLOT(addMesh()));
	QObject::connect(ui.actionSaveMatching, SIGNAL(triggered()), this, SLOT(saveMatchingResult()));
	QObject::connect(ui.actionLoadMatching, SIGNAL(triggered()), this, SLOT(loadMatchingResult()));

	////  Compute  ////
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(computeEigenfunction()));
	QObject::connect(ui.actionComputeBasis, SIGNAL(triggered()), this, SLOT(computeEditBasis()));
    QObject::connect(ui.actionComputeCurvatures, SIGNAL(triggered()), this, SLOT(computeCurvatures()));
	QObject::connect(ui.actionComputeHK, SIGNAL(triggered()), this, SLOT(computeHK()));
	QObject::connect(ui.actionComputeHKS, SIGNAL(triggered()), this, SLOT(computeHKS()));
	QObject::connect(ui.actionComputeHKSFeatures, SIGNAL(triggered()), this, SLOT(computeHKSFeatures()));
	QObject::connect(ui.actionComputeBiharmonic, SIGNAL(triggered()), this, SLOT(computeBiharmonic()));
	QObject::connect(ui.actionComputeGeodesics, SIGNAL(triggered()), this, SLOT(computeGeodesics()));
	QObject::connect(ui.actionComputeHeatTransfer, SIGNAL(triggered()), this, SLOT(computeHeatTransfer()));
	QObject::connect(ui.actionComputeVertNormals, SIGNAL(triggered()), this, SLOT(computeVertNormals()));
	QObject::connect(ui.actionComputeFaceNormals, SIGNAL(triggered()), this, SLOT(computeFaceNormals()));
    QObject::connect(ui.actionComputeHoleNeighbors, SIGNAL(triggered()), this, SLOT(computeHoleNeighbors()));

	////  Edit  ////
	QObject::connect(ui.actionClearHandles, SIGNAL(triggered()), this, SLOT(clearHandles()));
    QObject::connect(ui.actionListAttributes, SIGNAL(triggered()), this, SLOT(listMeshAttributes()));
	QObject::connect(ui.actionRevertCoordinate, SIGNAL(triggered()), this, SLOT(revertCoord()));
	QObject::connect(ui.actionNextCoordinate, SIGNAL(triggered()), this, SLOT(switchToNextCoordinate()));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
	QObject::connect(ui.actionAddNoise, SIGNAL(triggered()), this, SLOT(addNoise()));
	QObject::connect(ui.actionReconstructMHB, SIGNAL(triggered()), this, SLOT(reconstructMHB()));
	QObject::connect(ui.actionDeformSimple, SIGNAL(triggered()), this, SLOT(deformSimple()));
	QObject::connect(ui.actionDeformLaplace, SIGNAL(triggered()), this, SLOT(deformLaplace()));
	QObject::connect(ui.actionDeformLaplace2, SIGNAL(triggered()), this, SLOT(deformLaplace2()));
	QObject::connect(ui.actionDeformBiLaplace, SIGNAL(triggered()), this, SLOT(deformBiLaplace()));
	QObject::connect(ui.actionDeformMixedLaplace, SIGNAL(triggered()), this, SLOT(deformMixedLaplace()));
	QObject::connect(ui.actionDiffusionFlow, SIGNAL(triggered()), this, SLOT(diffusionFlow()));
	QObject::connect(ui.actionRunTests, SIGNAL(triggered()), this, SLOT(runTests()));
    QObject::connect(ui.actionFillHoles, SIGNAL(triggered()), this, SLOT(fillHoles()));
    QObject::connect(ui.actionFairFilledHoles, SIGNAL(triggered()), this, SLOT(holeFairingAll()));
    QObject::connect(ui.actionHoleFairingLS, SIGNAL(triggered()), this, SLOT(holeFairingLS()));
    QObject::connect(ui.actionHoleFairingFourierOMP, SIGNAL(triggered()), this, SLOT(holeFairingFourierOMP()));
    QObject::connect(ui.actionFourierLARS, SIGNAL(triggered()), this, SLOT(holeFairingLARS()));
    QObject::connect(ui.actionHoleEstimateCurvature, SIGNAL(triggered()), this, SLOT(holeEstimateCurvature()));
   
    QObject::connect(ui.actionSwitchMesh, SIGNAL(triggered()), this, SLOT(switchToNextMesh()));
    QObject::connect(ui.actionGenerateHoles, SIGNAL(triggered()), this, SLOT(generateHoles()));
    QObject::connect(ui.actionGenerateRingHoles, SIGNAL(triggered()), this, SLOT(generateRingHoles()));
    QObject::connect(ui.actionAutoGenHoles, SIGNAL(triggered()), this, SLOT(autoGenerateHoles()));
    QObject::connect(ui.actionDegradeHoles, SIGNAL(triggered()), this, SLOT(degradeHoles()));
    QObject::connect(ui.actionCutHoles, SIGNAL(triggered()), this, SLOT(cutHoles()));
    QObject::connect(ui.actionCutToSelected, SIGNAL(triggered()), this, SLOT(cutToSelected()));
    QObject::connect(ui.actionTriangulateHoles, SIGNAL(triggered()), this, SLOT(triangulateHoles()));
    QObject::connect(ui.actionRefineHoles, SIGNAL(triggered()), this, SLOT(refineHoles()));
    QObject::connect(ui.actionHoleFairingLeastSquares, SIGNAL(triggered()), this, SLOT(fairHoleLeastSquares()));

	////  Display  ////
	QObject::connect(ui.actionDisplayMesh, SIGNAL(triggered()), this, SLOT(setDisplayMesh()));
    QObject::connect(ui.actionChangeShadeMode, SIGNAL(triggered()), ui.glMeshWidget, SLOT(changeShadeMode()));
	QObject::connect(ui.actionDisplayWireframe, SIGNAL(triggered()), this, SLOT(setDisplayWireframe()));
	QObject::connect(ui.actionDisplayPointCloud, SIGNAL(triggered()), this, SLOT(setDisplayPointCloud()));
	QObject::connect(ui.actionDisplayNeighbors, SIGNAL(triggered()), this, SLOT(displayNeighborVertices()));
	QObject::connect(ui.actionShowFeatures, SIGNAL(triggered(bool)), this, SLOT(toggleShowFeatures(bool)));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered(bool)), this, SLOT(toggleShowRefPoint(bool)));
	QObject::connect(ui.actionShowSignature, SIGNAL(triggered(bool)), this, SLOT(toggleShowSignature(bool)));
	QObject::connect(ui.actionShowWireframeOverlay, SIGNAL(triggered(bool)), this, SLOT(toggleShowWireframeOverlay(bool)));
    QObject::connect(ui.actionShowHoles, SIGNAL(triggered(bool)), this, SLOT(toggleShowHoles(bool)));
    QObject::connect(ui.actionShowSurrounding, SIGNAL(triggered(bool)), this, SLOT(toggleShowSurrounding(bool)));
    QObject::connect(ui.actionShowBbox, SIGNAL(triggered(bool)), this, SLOT(toggleShowBoundingBox(bool)));
	QObject::connect(ui.actionShowColorLegend, SIGNAL(triggered(bool)), this, SLOT(toggleShowColorLegend(bool)));
	QObject::connect(ui.actionShowVectors, SIGNAL(triggered(bool)), this, SLOT(toggleShowLines(bool)));
	QObject::connect(ui.actionDrawMatching, SIGNAL(triggered(bool)), this, SLOT(toggleDrawMatching(bool)));
	QObject::connect(ui.actionShowMatchingLines, SIGNAL(triggered(bool)), this, SLOT(toggleShowMatchingLines(bool)));
	QObject::connect(ui.actionDrawRegistration, SIGNAL(triggered(bool)), this, SLOT(toggleDrawRegistration(bool)));	
	QObject::connect(ui.actionDiffPosition, SIGNAL(triggered()), this, SLOT(displayDiffPosition()));
	
    QObject::connect(ui.actionCapture, SIGNAL(triggered()), this, SLOT(captureGL()));
	QObject::connect(ui.actionCaptureAs, SIGNAL(triggered()), this, SLOT(captureGLAs()));

	////  Task  ////
	QObject::connect(ui.actionTaskRegistration, SIGNAL(triggered()), this, SLOT(setTaskRegistration()));
	QObject::connect(ui.actionTaskEditing, SIGNAL(triggered()), this, SLOT(setTaskEditing()));

	////  Register  ////
	QObject::connect(ui.actionRegisterAutomatic, SIGNAL(triggered()), this, SLOT(registerAutomatic()));
	QObject::connect(ui.actionDetectFeatures, SIGNAL(triggered()), this, SLOT(detectFeatures()));
	QObject::connect(ui.actionMatchFeatures, SIGNAL(triggered()), this, SLOT(matchFeatures()));
	QObject::connect(ui.actionRegisterStep, SIGNAL(triggered()), this, SLOT(registerStep()));
	QObject::connect(ui.actionRegisterFull, SIGNAL(triggered()), this, SLOT(registerFull()));
	QObject::connect(ui.actionRegisterTest, SIGNAL(triggered()), this, SLOT(registerTest()));

	////  Tools  ////
	QObject::connect(ui.actionExploreScreenshots, SIGNAL(triggered()), this, SLOT(openSreenshotLocation()));
	QObject::connect(ui.actionExploreOutput, SIGNAL(triggered()), this, SLOT(openOutputLocation()));
	

	////////    mShapeEditor    ////////
	QObject::connect(&mShapeEditor, SIGNAL(approxStepsChanged(int, int)), this, SLOT(resizeApproxSlider(int, int)));
    QObject::connect(&mShapeEditor, SIGNAL(meshSignatureAdded()), this, SLOT(updateMenuDisplaySignature()));
    QObject::connect(&mShapeEditor, SIGNAL(meshPointFeatureChanged()), this, SLOT(updateMenuDisplayFeatures()));
    QObject::connect(&mShapeEditor, SIGNAL(meshLineFeatureChanged()), this, SLOT(updateMenuDisplayLines()));
    QObject::connect(&mShapeEditor, SIGNAL(showSignature(QString)), this, SLOT(displaySignature(QString)));
	QObject::connect(&mShapeEditor, SIGNAL(coordinateSelected(int, int)), this, SLOT(visualizeCompression(int, int)));
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

	case Qt::Key_B:
		break;

	case Qt::Key_C:
		mColorMapType = ZGeom::ColorMapType(int(mColorMapType + 1) % (int)ZGeom::CM_COUNT);
		break;

	case Qt::Key_D:
		setEditModeDrag();
		break;

	case Qt::Key_E:
		if (event->modifiers() & Qt::AltModifier) {
			openOutputLocation();
		} else {
			if (mDeformType == DEFORM_Simple)
				deformSimple();
			else if (mDeformType == DEFORM_Laplace)
				deformLaplace();
		}
		break;

	case Qt::Key_F:
		if (event->modifiers() & Qt::AltModifier)
			toggleShowFeatures();
		break;

	case Qt::Key_G:
		break;

	case Qt::Key_J:
		break;

	case Qt::Key_L:
		if (event->modifiers() & Qt::AltModifier)
			toggleShowMatchingLines();
		else if (event->modifiers() & Qt::ControlModifier)
			listMeshAttributes();
		break;

	case Qt::Key_M:
		setEditModeMove();
		break;

	case Qt::Key_N:
        if (event->modifiers() & Qt::AltModifier) {            
            if (event->modifiers() & Qt::ShiftModifier)
                switchToPreviousMesh();
            else 
                switchToNextMesh();
        }
        else { 
            if (event->modifiers() & Qt::ShiftModifier)
                switchToPrevCoordinate();
            else
                switchToNextCoordinate(); 
        }
		break;

	case Qt::Key_P:
		setEditModePick();
		if (!ui.glMeshWidget->m_bShowRefPoint)
			toggleShowRefPoint();
		break;

	case Qt::Key_Q:
		if (event->modifiers() & Qt::AltModifier && event->modifiers() & Qt::ShiftModifier)
			captureGLAs();
		else if (event->modifiers() & Qt::AltModifier)
			captureGL();
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
		if (mRenderManagers[0].displayType == RenderSettings::Mesh)
			setDisplayWireframe();
		else if (mRenderManagers[0].displayType == RenderSettings::Wireframe)
			setDisplayPointCloud();
		else 
            setDisplayMesh();
		break;
	}	
	
	case Qt::Key_X:
		if (event->modifiers() & Qt::ShiftModifier)
			;
		else ;
		break;

	case Qt::Key_Y:
		if (event->modifiers() & Qt::ShiftModifier)
			;
		else ;
		break;

	case Qt::Key_Z:
		if (event->modifiers() & Qt::ShiftModifier)
			;
		else;
		break;

	case Qt::Key_BracketLeft:
		showFiner();
		break;

	case Qt::Key_BracketRight:
		showCoarser();
		break;

	case Qt::Key_Minus:
		mDiffMax -= 0.1;
		qout.output(QString().sprintf("diff max = %f", mDiffMax), OUT_STATUS);
		visualizeCompression(mSelectedApprox, mCoordIdx);
		break;

	case Qt::Key_Equal:
		mDiffMax += 0.1;
		qout.output(QString().sprintf("diff max = %f", mDiffMax), OUT_STATUS);
		visualizeCompression(mSelectedApprox, mCoordIdx);
		break;

	default: QWidget::keyPressEvent(event);
	}
}

bool QZGeometryWindow::initialize(const std::string& mesh_list_name)
{
	qout.output("******** Welcome ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);
	qout.output('*', 24, OUT_CONSOLE);

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

    std::vector<MeshHelper*> vpHelper;
    for (MeshHelper& mh : mMeshHelper) vpHelper.push_back(&mh);
    ui.glMeshWidget->setup(vpHelper, mRenderManagers, &mShapeMatcher);

    /* compute and decompose mesh Laplacians */
    //computeLaplacian(Umbrella);
    //computeLaplacian(NormalizedUmbrella);	
    //computeLaplacian(CotFormula);
    //computeLaplacian(SymCot);
    //computeLaplacian(Anisotropic1); 	
    //computeLaplacian(Anisotropic2);
    //setLaplacianType("Anisotropic1");

	if (g_task == TASK_REGISTRATION) {
		registerPreprocess();
	}
	if (g_task == TASK_EDITING) {
		mShapeEditor.init(mMeshHelper[0]);
		mShapeEditor.runTests();
	}

    updateMenuDisplaySignature();
    displaySignature(CMesh::StrAttrColorSigDefault.c_str());
	return true;
}

void QZGeometryWindow::loadInitialMeshes(const std::string& mesh_list_name)
{
	runtime_assert (fileExist(mesh_list_name),  "Cannot open file mesh list file!");
	std::ifstream meshfiles(mesh_list_name);	
	vector<std::string> vMeshFiles;
	while (!meshfiles.eof()) {
		std::string meshFileName;
		getline(meshfiles, meshFileName);
		for (auto iter = meshFileName.begin(); iter != meshFileName.end(); ) {
			if (*iter == ' ' || *iter == '\t') iter = meshFileName.erase(iter);
			else ++iter;
		}
		if (meshFileName == "") continue;
		if (meshFileName[0] == '#') continue;
        runtime_assert(fileExist(meshFileName), "Cannot open file " + meshFileName);
		vMeshFiles.push_back(meshFileName);
	}
	meshfiles.close();
	if (vMeshFiles.size() < mMeshCount)
		throw std::runtime_error("Not enough meshes in mesh list!");

	allocateStorage(mMeshCount);
	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
        loadMesh(vMeshFiles[obj], obj); 
	});	

	/* ---- update mesh-dependent ui ---- */
	if (mMeshCount >= 1) {
		ui.spinBox1->setMinimum(0);
        ui.spinBox1->setMaximum(getMesh(0)->vertCount() - 1);
		ui.horizontalSlider1->setMinimum(0);
		ui.horizontalSlider1->setMaximum(getMesh(0)->vertCount() - 1);
		ui.spinBox1->setValue(0);	
	}
	if (mMeshCount >= 2) {		
		ui.spinBox2->setMinimum(0);
		ui.spinBox2->setMaximum(getMesh(1)->vertCount()-1);
		ui.horizontalSlider2->setMinimum(0);
		ui.horizontalSlider2->setMaximum(getMesh(1)->vertCount()-1);
		ui.spinBox2->setValue(0);
	}
	ui.glMeshWidget->fieldView(getMesh(0)->getCenter(), getMesh(0)->getBoundingBox());

	mRenderManagers[0].selected = true;
	mObjInFocus = 0;
}

void QZGeometryWindow::loadMesh(std::string mesh_filename, int obj)
{
    assert(obj < mMeshCount);
    CMesh *newMesh = new CMesh();
    newMesh->load(mesh_filename);
    newMesh->scaleToUnitBox();
    ZGeom::gatherMeshStatistics(*newMesh);
    Colorf meshColor(ZGeom::MeshPresetColors[obj % 3]);
    newMesh->setDefaultColor(meshColor);
    newMesh->initNamedCoordinates();

    auto center = newMesh->getCenter();
    auto bbox = newMesh->getBoundingBox();
    qout.output(QString().sprintf("Load mesh: %s; Size: %d", newMesh->getMeshName().c_str(), newMesh->vertCount()), OUT_TERMINAL);
    qout.output(QString().sprintf("Center: (%f, %f, %f)\nDimension: (%f, %f, %f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_TERMINAL);
    if (newMesh->hasBoundary()) std::cout << "Original mesh has holes!\n";

    mMeshHelper[obj].init(newMesh);
    mRenderManagers[obj].mActiveColorSignatureName = CMesh::StrAttrColorSigDefault;
}

void QZGeometryWindow::addMesh()
{
	QStringList filenames =  QFileDialog::getOpenFileNames(this, "Select one or more mesh files to open",
														   "../../Data/", "Meshes (*.obj *.off *.ply)");
	int cur_obj = mMeshCount;
	allocateStorage(++mMeshCount);

    CMesh *newMesh = new CMesh();
    newMesh->load(filenames.front().toStdString());
    newMesh->scaleToUnitBox();
    ZGeom::gatherMeshStatistics(*newMesh);
    auto center = newMesh->getCenter();
    auto bbox = newMesh->getBoundingBox();
    qout.output(QString().sprintf("Load mesh: %s; Size: %d", newMesh->getMeshName().c_str(), newMesh->vertCount()), OUT_TERMINAL);
    qout.output(QString().sprintf("Center: (%f, %f, %f)\nDimension: (%f, %f, %f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_TERMINAL);
    Colorf meshColor(ZGeom::MeshPresetColors[cur_obj % 3]);
    newMesh->setDefaultColor(meshColor);
    newMesh->initNamedCoordinates();

	mMeshHelper[cur_obj].init(newMesh);
    mRenderManagers[cur_obj].mActiveColorSignatureName = CMesh::StrAttrColorSigDefault;
	mRenderManagers[cur_obj].selected = true;

	if (cur_obj == 0) {
		ui.glMeshWidget->fieldView(getMesh(0)->getCenter(), getMesh(0)->getBoundingBox());
		ui.spinBox1->setMinimum(0);
		ui.spinBox1->setMaximum(getMesh(0)->vertCount() - 1);
		ui.horizontalSlider1->setMinimum(0);
		ui.horizontalSlider1->setMaximum(getMesh(0)->vertCount() - 1);
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::clone()
{
	if (mMeshCount != 1) {
		std::cout << "One and only one mesh should exist for cloning" << std::endl;
		return;
	}

	allocateStorage(2);
    CMesh *newMesh = new CMesh(*getMesh(0));
	mMeshHelper[1].init(newMesh);

    qout.output(QString().sprintf("Mesh %s constructed! Size: %d", getMesh(1)->getMeshName().c_str(), getMesh(1)->vertCount()));
	ui.glMeshWidget->update();
}

void QZGeometryWindow::deformSimple()
{
	mShapeEditor.deformSimple();

	mDeformType = DEFORM_Simple;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::deformLaplace()
{
	mShapeEditor.deformLaplacian();
	mDeformType = DEFORM_Laplace;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::deformBiLaplace()
{
	mShapeEditor.deformBiLaplacian();
	mDeformType = DEFORM_BiLaplace;
	ui.glMeshWidget->update();
	setEditModeMove();
}


void QZGeometryWindow::deformMixedLaplace()
{
	double ks = 1.0, kb = 1.0;
	mShapeEditor.deformMixedLaplacian(ks, kb);
	mDeformType = DEFORM_Shell;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::selectObject( int index )
{
	QString text = ui.boxObjSelect->itemText(index);	
	for (auto iter = mRenderManagers.begin(); iter != mRenderManagers.end(); ++iter) iter->selected = false;	
	
	if (text == "1") {			 
		if (mRenderManagers.size() >= 1) mRenderManagers[0].selected = true;
		mObjInFocus = 0;
	}
	else if (text == "2") {
		if (mRenderManagers.size() >= 2) mRenderManagers[1].selected = true;
		mObjInFocus = 1;
	}
	else if (text == "All") {
		for (auto& rs : mRenderManagers) rs.selected = true;	 
		mObjInFocus = 0;
	}
	else if (text == "None") {
		mObjInFocus = -1;
	}

	qout.output("Selected object(s): " + text);
}

void QZGeometryWindow::setMesh1RefPoint( int vn )
{
	if (mMeshCount < 1) return;
	mMeshHelper[0].setRefPointIndex(vn);
	updateReferenceMove(0);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setMesh2RefPoint( int vn )
{
	if (mMeshCount < 2) return;

	mMeshHelper[1].setRefPointIndex(vn);
	updateReferenceMove(1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setCommonParameter( int p )
{
	mCommonParameter = p;
	int sliderCenter = ui.horizontalSliderParamter->maximum()/2;

	if (mLastOperation == Compute_HKS || mLastOperation == Compute_HK 
		|| mLastOperation == Compute_MHWS || mLastOperation == Compute_MHW
		|| mLastOperation == Compute_SGWS || mLastOperation == Compute_SGW)
	{
		double time_scale;
		if (mCommonParameter <= sliderCenter) 
			time_scale = std::exp(std::log(gSettings.DEFUALT_HK_TIMESCALE / gSettings.MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)sliderCenter) + std::log(gSettings.MIN_HK_TIMESCALE));
		else 
			time_scale = std::exp(std::log(gSettings.MAX_HK_TIMESCALE / gSettings.DEFUALT_HK_TIMESCALE) * (double(mCommonParameter - sliderCenter) / sliderCenter) + std::log(gSettings.DEFUALT_HK_TIMESCALE)); 
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

void QZGeometryWindow::setDisplayPointCloud()
{
	ui.actionDisplayPointCloud->setChecked(true);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(false);

	for ( auto& rm : mRenderManagers ) {
		rm.displayType = RenderSettings::PointCloud;
		rm.glPolygonMode = GL_POINT;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayWireframe()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(true);
	ui.actionDisplayMesh->setChecked(false);
	
	for ( auto& rm : mRenderManagers ) {
		rm.displayType = RenderSettings::Wireframe;
		rm.glPolygonMode = GL_LINE;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayMesh()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(true);

	for ( auto& rs : mRenderManagers) {
		rs.displayType = RenderSettings::Mesh;
		rs.glPolygonMode = GL_FILL;
	}

	ui.glMeshWidget->update();
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

void QZGeometryWindow::toggleShowWireframeOverlay(bool show /*= false*/)
{
	bool toShow = !ui.glMeshWidget->m_bShowWireframeOverlay;
	ui.glMeshWidget->m_bShowWireframeOverlay = toShow;
	ui.actionShowWireframeOverlay->setChecked(toShow);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowBoundingBox(bool show /*= false*/)
{
	bool toShow = !ui.glMeshWidget->m_bShowBoundingBox;
	ui.glMeshWidget->m_bShowBoundingBox = toShow;
	ui.actionShowBbox->setChecked(toShow);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowLines( bool show /*= false*/ )
{
	bool bToShow = !ui.glMeshWidget->m_bShowLines;
	ui.glMeshWidget->m_bShowLines = bToShow;
	ui.actionShowVectors->setChecked(bToShow);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowHoles(bool show /*= false*/)
{
    bool bToShow = !ui.glMeshWidget->m_bShowHoles;
    ui.glMeshWidget->m_bShowHoles = bToShow;
    ui.actionShowHoles->setChecked(bToShow);

    ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowSurrounding(bool show /* = false */)
{
    bool bToShow = !ui.glMeshWidget->m_bShowSurrounding;
    ui.glMeshWidget->m_bShowSurrounding = bToShow;
    ui.actionShowSurrounding->setChecked(bToShow);

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

void QZGeometryWindow::computeCurvatures()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
        ZGeom::calMeshAttrMeanGaussCurvatures(*getMesh(obj));
        vector<double> vCM = ZGeom::getMeshMeanCurvatures(*getMesh(obj));
        vector<double> vCG = ZGeom::getMeshGaussCurvatures(*getMesh(obj));
//         vector<double> vCP1 = getMesh(obj)->calPrincipalCurvature(1);
//         vector<double> vCP2 = getMesh(obj)->calPrincipalCurvature(2);

        ColorSignature colorCM(vCM, mColorMapType);
        ColorSignature colorCG(vCG, mColorMapType);
//         ColorSignature colorCP1(vCP1, mColorMapType);
//         ColorSignature colorCP2(vCP2, mColorMapType);

        getMesh(obj)->addColorSigAttr("color_mean_curvature", colorCM);
        getMesh(obj)->addColorSigAttr("color_gauss_curvature", colorCG);
//         getMesh(obj)->addColorAttr("color_principal_curvature_1", colorCP1);
//         getMesh(obj)->addColorAttr("color_principal_curvature_2", colorCP2);


        auto mm1 = std::minmax_element(vCM.begin(), vCM.end());
        auto mm2 = std::minmax_element(vCG.begin(), vCG.end());
        qout.output(QString("- mean curvature -  min: %1, max: %2").arg(*mm1.first).arg(*mm1.second), OUT_TERMINAL);
        qout.output(QString("- Gauss curvature -  min: %1, max: %2").arg(QString::number(*mm2.first), QString::number(*mm2.second)), OUT_TERMINAL);
	}

	updateMenuDisplaySignature();
    displaySignature("color_mean_curvature");
	qout.output("Visualize mean curvature");	
}

void QZGeometryWindow::updateReferenceMove( int obj )
{
	MeshHelper& mp = mMeshHelper[obj]; 

	double unitMove = (mp.getMesh()->getBoundingBox().x + mp.getMesh()->getBoundingBox().y + mp.getMesh()->getBoundingBox().z)/300.0;
	ZGeom::Vec3d originalPos = mp.getMesh()->vert(mp.getRefPointIndex())->pos();
	
	mp.setRefPointPosition(originalPos.x, originalPos.y, originalPos.z);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::registerPreprocess()
{
	if (g_task != TASK_REGISTRATION || mMeshCount != 2) return;

	computeFunctionMaps(40);
	mShapeMatcher.initialize(&mMeshHelper[0], &mMeshHelper[1], g_engineWrapper.getEngine());
	std::string rand_data_file = g_configMgr.getConfigValue("RAND_DATA_FILE");
	mShapeMatcher.readInRandPair(rand_data_file);

	mShapeMatcher.setRegistrationLevels(1);
	registerTest();
//	evalDistance();
}

void QZGeometryWindow::reconstructMHB()
{
	int sliderCenter = ui.horizontalSliderParamter->maximum() / 2;
	double ratio = std::min((double)mCommonParameter/sliderCenter, 1.0);
	int nEig = mMeshHelper[0].getMHB(CotFormula).eigVecCount() * ratio;
	double avgLen = getMesh(0)->getAvgEdgeLength();

	mShapeEditor.fourierReconstruct(nEig);
    std::cout << "Reconstruct with " << nEig << " eigenvectors" << std::endl;

    ZGeom::VecNd vPosDiff = mShapeEditor.getOldMeshCoord().vertDifference(getMesh(0)->getVertCoordinates());
    for (double& v : vPosDiff) v /= avgLen;
    std::cout << "Avg Error as ratio of AEL: " << vPosDiff.mean() << std::endl;

	addColorSignature(0, vPosDiff.toStdVector(), StrAttrColorPosDiff);
	displaySignature(StrAttrColorPosDiff.c_str());
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayNeighborVertices()
{
	int sliderCenter = ui.horizontalSliderParamter->maximum()/2;
	int ring = (mCommonParameter > sliderCenter) ? (mCommonParameter - sliderCenter) : 1;

	int ref = mMeshHelper[0].getRefPointIndex();
	std::vector<int> vn = mMeshHelper[0].getMesh()->getVertNeighborVerts(ref, ring, false);
	MeshFeatureList *mfl = new MeshFeatureList;

	for (auto iter = vn.begin(); iter != vn.end(); ++iter) {
		mfl->addFeature(new MeshFeature(*iter));
	}
	getMesh(0)->addAttrMeshFeatures(*mfl, StrAttrFeatureNeighbors);

	if (!ui.actionShowFeatures->isChecked()) toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeHoleNeighbors()
{
    CMesh& mesh = *getMesh(0);
    MeshRegion* hole = nullptr;
    if (mesh.hasAttr(ZGeom::StrAttrMeshHoleRegions)) {
        auto &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
        if (!vHoles.empty()) hole = &vHoles[0];
    }    
    if (hole == nullptr && mesh.hasAttr(StrAttrManualHoles)) {
        auto &vHoles = mesh.getAttrValue<vector<MeshRegion>>(StrAttrManualHoles);
        if (!vHoles.empty()) hole = &vHoles[0];        
    }
    if (hole == nullptr) {
        std::cout << "Mesh has no hole!" << std::endl;
        return;
    }

    int ring = 3;
    bool ok;
    ring = QInputDialog::getInt(this, tr("Input surrounding ring"),
        tr("Ring:"), ring, 1, 50, 1, &ok);
    if (!ok) return;

    vector<int> surrounding_verts = ZGeom::meshRegionSurroundingVerts(mesh, *hole, ring);
    mesh.addAttrMeshFeatures(MeshFeatureList(surrounding_verts, ZGeom::ColorMagenta), "hole_ring_neighbor_verts");
    updateMenuDisplayFeatures();

    vector<int> surrounding_faces = ZGeom::getFaceIdxEncompassedByVerts(mesh, surrounding_verts);
    mesh.addAttr<vector<int>>(surrounding_faces, StrAttrHoleSurroundingFaces, AR_UNIFORM, AT_VEC_INT);
    
    ui.glMeshWidget->update();
}

void QZGeometryWindow::computeEigenfunction()
{
	LaplacianType lapType = mActiveLalacian;
	int sliderCenter = ui.horizontalSliderParamter->maximum() / 2;
	int select_eig = (mCommonParameter >= sliderCenter) ? (mCommonParameter - sliderCenter + 1) : 1;
	
	for (int i = 0; i < mMeshCount; ++i) {
		MeshHelper& mp = mMeshHelper[i];
		vector<double> eigVec = mp.getMHB(lapType).getEigVec(select_eig).toStdVector();
        getMesh(i)->addColorSigAttr(StrAttrColorEigenFunction, ColorSignature(eigVec, mColorMapType));
	}

	displaySignature(StrAttrColorEigenFunction.c_str());
	mLastOperation = Compute_Eig_Func;
	qout.output("Show eigenfunction" + Int2String(select_eig), OUT_CONSOLE);
	updateMenuDisplaySignature();
}

void QZGeometryWindow::computeHK()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int obj = 0; obj < mMeshCount; ++obj) {
		MeshHelper& mp = mMeshHelper[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();
		const ZGeom::EigenSystem &mhb = mp.getMHB(mActiveLalacian);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int vIdx) {
            values[vIdx] = ZGeom::calHK(mhb, refPoint, vIdx, time_scale);
		});

		addColorSignature(obj, values, StrAttrColorHK);
	}

	displaySignature(StrAttrColorHK.c_str());
	qout.output(QString().sprintf("HK with timescale: %f", time_scale));
	updateMenuDisplaySignature();
	mLastOperation = Compute_HK;
}

void QZGeometryWindow::computeHKS()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int i = 0; i < mMeshCount; ++i) {
		MeshHelper& mp = mMeshHelper[i];
		const int meshSize = mp.getMesh()->vertCount();
		const ZGeom::EigenSystem& mhb = mp.getMHB(mActiveLalacian);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int k) {
            values[k] = ZGeom::calHK(mhb, k, k, time_scale);
		});

		addColorSignature(i, values, StrAttrColorHKS);
	}
	
	displaySignature(StrAttrColorHKS.c_str());
	qout.output(QString().sprintf("HKS with timescale: %f", time_scale));
	updateMenuDisplaySignature();
	mLastOperation = Compute_HKS;
}

void QZGeometryWindow::computeBiharmonic()
{
	for (int obj = 0; obj < mMeshCount; ++obj)
	{
		MeshHelper& mp = mMeshHelper[obj];
		const int vertCount = getMesh(obj)->vertCount();
		const int refPoint = mp.getRefPointIndex();
        const ZGeom::EigenSystem& es = mp.getMHB(CotFormula);

		std::vector<double> vVals(vertCount);
		for (int vIdx = 0; vIdx < vertCount; ++vIdx) {
			vVals[vIdx] = ZGeom::calBiharmonicDist(es, refPoint, vIdx);
		}

		addColorSignature(obj, vVals, StrAttrColorBiharmonic);
	}

	displaySignature(StrAttrColorBiharmonic.c_str());
	updateMenuDisplaySignature();
	mLastOperation = Compute_Biharmonic;
}

void QZGeometryWindow::computeHKSFeatures()
{
    std::vector<double> vTimes{ 10, 30, 90, 270 };

	for (int i = 0; i < mMeshCount; ++i) {

	}
	
	if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures();
}

void QZGeometryWindow::repeatOperation()
{
	switch(mLastOperation)
	{
	case Compute_Eig_Func:
		computeEigenfunction(); break;

	case Compute_Edit_Basis:
		computeEditBasis(); break;
	
	case Compute_HKS:
		computeHKS(); break;

	case Compute_HK:
		computeHK(); break;

	case Compute_Biharmonic:
		computeBiharmonic(); break;

	case Compute_Heat:
		computeHeatTransfer(); break;
	}
}

void QZGeometryWindow::displayDiffPosition()
{
	runtime_assert(mMeshCount >= 2 && getMesh(0)->vertCount() == getMesh(1)->vertCount());
	int size = getMesh(0)->vertCount();
	std::vector<double> vDiff;
	vDiff.resize(size);

    for (int i = 0; i < getMesh(0)->vertCount(); ++i) {
        vDiff[i] = (getMesh(0)->vert(i)->pos() - getMesh(1)->vert(i)->pos()).length() / getMesh(0)->getAvgEdgeLength();
	}

	addColorSignature(0, vDiff, StrAttrColorPosDiff);
	displaySignature(StrAttrColorPosDiff.c_str());
	updateMenuDisplaySignature();
}

void QZGeometryWindow::displaySignature(QString sigName )
{
	for (int i = 0; i < mMeshCount; ++i) {
        if (getMesh(i)->hasAttr(sigName.toStdString())) {
            mRenderManagers[i].mActiveColorSignatureName = sigName.toStdString();            
        }        
	}
    
    for (QAction *qa : m_actionDisplaySignatures) {
        if (qa->text() == sigName) qa->setChecked(true);
        else qa->setChecked(false);        
    }

	if (ui.glMeshWidget->m_bShowSignature == false) toggleShowSignature(true);	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayFeature( QString pointFeatureName )
{
	for (int obj = 0; obj < mMeshCount; ++obj) 
    {
		if (!isMeshSelected(obj)) continue;
        auto &activePointFeatures = mRenderManagers[obj].mActivePointFeatures;
        std::string newPointFeature = pointFeatureName.toStdString();
        QAction *activeAction = nullptr;
        for (QAction *qa : m_actionDisplayFeatures) {
            if (qa->text() == pointFeatureName) {
                activeAction = qa;
                break;
            }
        }
        logic_assert(activeAction != nullptr, "Should have an active action");

        if (!setHas(activePointFeatures, newPointFeature)) {
            activePointFeatures.insert(newPointFeature);
            activeAction->setChecked(true);
            if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures(true);
        }
        else {
            activePointFeatures.erase(newPointFeature);
            activeAction->setChecked(false);            
        }
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayLine(QString lineFeatureName)
{
    for (int obj = 0; obj < mMeshCount; ++obj) {
        if (!isMeshSelected(obj)) continue;
        auto &activeLineFeatures = mRenderManagers[obj].mActiveLineNames;
        std::string newLineFeature = lineFeatureName.toStdString();
        QAction *activeAction = nullptr;
        for (QAction *qa : m_actionDisplayLines) {
            if (qa->text() == lineFeatureName) {
                activeAction = qa;
                break;
            }
        }
        logic_assert(activeAction != nullptr, "Should have an active action");

        if (!setHas(activeLineFeatures, newLineFeature)) {
            activeLineFeatures.insert(newLineFeature);
            activeAction->setChecked(true);
            if (!ui.glMeshWidget->m_bShowLines) toggleShowLines(true);
        }
        else {
            activeLineFeatures.erase(newLineFeature);
            activeAction->setChecked(false);
        }
    }

    ui.glMeshWidget->update();
}

void QZGeometryWindow::updateMenuDisplaySignature()
{
	int obj = (mObjInFocus <= 0 ? 0 : mObjInFocus);
    vector<AttrVertColors*> vColorAttributes = getMesh(obj)->getColorAttrList();
	for (QAction* qa : m_actionDisplaySignatures) {
    	ui.menuDisplaySignatures->removeAction(qa);
		delete qa;	
	}
    m_actionDisplaySignatures.clear();

	for (AttrVertColors* attr : vColorAttributes) {
        QString action_name = attr->attrName().c_str();
        QAction* newDisplayAction = new QAction(action_name, this);
        m_actionDisplaySignatures.push_back(newDisplayAction);
        ui.menuDisplaySignatures->addAction(newDisplayAction);
        newDisplayAction->setCheckable(true);
        m_signatureSignalMapper->setMapping(newDisplayAction, action_name);
        QObject::connect(newDisplayAction, SIGNAL(triggered()), m_signatureSignalMapper, SLOT(map()));
	}
}

void QZGeometryWindow::updateMenuDisplayFeatures()
{
	int obj = (mObjInFocus <= 0 ? 0 : mObjInFocus);
    std::vector<AttrMeshFeatures*> vFeatureAttr = getMesh(obj)->getMeshFeatureList();
	for (QAction* qa : m_actionDisplayFeatures) {
    	ui.menuDisplayFeatures->removeAction(qa);
		delete qa;		
	}
    m_actionDisplayFeatures.clear();

	for (AttrMeshFeatures* attr : vFeatureAttr) {
        QString action_name = attr->attrName().c_str();
        QAction* newDisplayAction = new QAction(action_name, this);
		m_actionDisplayFeatures.push_back(newDisplayAction);
        ui.menuDisplayFeatures->addAction(newDisplayAction);
        newDisplayAction->setCheckable(true);
        m_featureSignalMapper->setMapping(newDisplayAction, action_name);
		QObject::connect(newDisplayAction, SIGNAL(triggered()), m_featureSignalMapper, SLOT(map()));		
	}
}

void QZGeometryWindow::updateMenuDisplayLines()
{
    int obj = (mObjInFocus <= 0 ? 0 : mObjInFocus);
    vector<AttrMeshLines*> vLineAttr = getMesh(obj)->getMeshLineList();
    for (QAction *qa : m_actionDisplayLines) {
        ui.menuDisplayFeatures->removeAction(qa);
        delete qa;
    }
    m_actionDisplayLines.clear();

    for (AttrMeshLines* attr : vLineAttr) {
        QString action_name = attr->attrName().c_str();
        QAction* newDisplayAction = new QAction(action_name, this);
        m_actionDisplayLines.push_back(newDisplayAction);
        ui.menuDisplayLines->addAction(newDisplayAction);
        newDisplayAction->setCheckable(true);
        m_linesSignalMapper->setMapping(newDisplayAction, action_name);
        QObject::connect(newDisplayAction, SIGNAL(triggered()), m_linesSignalMapper, SLOT(map()));
    }
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
    //this->buildHierarchy();
	this->detectFeatures();
	this->matchFeatures();
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

				if (getMesh(1)->isInNeighborRing(iter2->m_index, iter1->m_index, 2)) {
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
				if (getMesh(0)->isInNeighborRing(iter1->m_index, iter2->m_index, 2)) {
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
				if (getMesh(1)->isInNeighborRing(iter1->m_index, iter2->m_index, 2)) {
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
	using namespace std;

	int force_matching = g_configMgr.getConfigValueInt("FORCE_MATCHING");
	string string_override = "";
	bool alreadyMached = false;

#if 0	
	if (force_matching == 1)
	{
		if (getMesh(0)->getMeshName() == "horse0")
		{
			string_override = g_configMgr.getConfigValue("HORSE0_FEATURE_OVERRIDE");
		}
		else if (getMesh(0)->getMeshName() == "eight")
		{
			string_override = g_configMgr.getConfigValue("EIGHT_FEATURE_OVERRIDE");
		}
		if (string_override != "")
		{
			vector<int> idx_override = splitStringToInt(string_override);
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
		if (getMesh(0)->getMeshName() == "horse0")
		{
			string_override = g_configMgr.getConfigValue("MATCHING_HORSE_FILE");
		}
		else if (getMesh(0)->getMeshName() == "eight")
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
	std::ofstream ofstr(log_filename.c_str(), std::ios::trunc);
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
		matchScore = ShapeMatcher::TensorGraphMatching6(g_engineWrapper.getEngine(), &mMeshHelper[0], &mMeshHelper[1], vftFine1, vftFine2, vPairs, tensor_matching_timescasle, matching_thresh_2, /*verbose=*/true);
		//matchScore = DiffusionShapeMatcher::TensorMatchingExt(m_ep, &vMP[0], &vMP[1], vFeatures1, vFeatures2, vPairs, 0, vPara, cout, true);

		if (1 == g_configMgr.getConfigValueInt("GROUND_TRUTH_AVAILABLE")) {
			vector<double> vTimes;
			vTimes.push_back(20); vTimes.push_back(40); vTimes.push_back(80); vTimes.push_back(160); vTimes.push_back(320);
			for (auto iter = vPairs.begin(); iter != vPairs.end(); ) {
				if (!getMesh(1)->isInNeighborRing(iter->m_idx1, iter->m_idx2, 2))
					iter->m_note = -1;

				double dissim = mShapeMatcher.calPointHksDissimilarity(&mMeshHelper[0], &mMeshHelper[1], iter->m_idx1, iter->m_idx2, vTimes, 1);
				if (dissim > 0.18) {
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
#if 0
	using namespace std;

	string log_filename = g_configMgr.getConfigValue("REGISTER_OUTPUT_FILE");

	qout.output(QString().sprintf("-- Register level %d --", mShapeMatcher.getAlreadyRegisteredLevel() - 1));
	
	ofstream ofstr;	
	if (mShapeMatcher.getAlreadyRegisteredLevel() == mShapeMatcher.getTotalRegistrationLevels())
		ofstr.open(log_filename.c_str(), ios::trunc);
	else
		ofstr.open(log_filename.c_str(), ios::app);

	int regMethod = g_configMgr.getConfigValueInt("REGISTRATION_METHOD");

	double time_elapsed = time_call_sec([&]{
		if (regMethod == 1) mShapeMatcher.refineRegister(ofstr);
		else if (regMethod == 2) mShapeMatcher.refineRegister2(ofstr);
	});

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

	if (!ui.glMeshWidget->m_bDrawRegistration)
		toggleDrawRegistration();
	ui.glMeshWidget->update();
#endif
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

void QZGeometryWindow::decomposeSingleLaplacian( int obj, int nEigVec, LaplacianType laplacianType /*= CotFormula*/ )
{
	MeshHelper& mp = mMeshHelper[obj];
	const CMesh& mesh = *getMesh(obj);
	const int vertCount = mesh.vertCount();
	if (nEigVec >= vertCount) nEigVec = vertCount - 1;

	if (!mp.getMHB(laplacianType).empty()) return;
	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "cache/" + mp.getMesh()->getMeshName() + ".mhb." + s_idx;
	
	if (gSettings.LOAD_MHB_CACHE && fileExist(pathMHB))	// MHB cache available for the current mesh
	{
		std::ifstream ifs(pathMHB.c_str());
		mp.loadMHB(pathMHB, laplacianType);
		ifs.close();
	}
	else // need to compute Laplacian and to cache
	{
		mp.decomposeLaplacian(nEigVec, laplacianType);
		mp.saveMHB(pathMHB, laplacianType);
		qout.output("MHB saved to " + pathMHB, OUT_CONSOLE);
	}

	std::cout << "Min EigVal: " << mp.getMHB(laplacianType).getAllEigVals().front() 
			  << "; Max EigVal: " << mp.getMHB(laplacianType).getAllEigVals().back() << std::endl;
}

void QZGeometryWindow::saveSignature()
{
	if (!getMesh(0)->hasAttr(StrAttrOriginalSignature)) {
		qout.output("No signature available", OUT_MSGBOX);
		return;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Signature to File"),
		"./output/signature.txt",
		tr("Text Files (*.txt *.dat)"));
    const std::vector<double>& vSig = getMesh(0)->getAttrValue<vector<double>>(StrAttrOriginalSignature);

	vector2file<double>(fileName.toStdString(), vSig);
}

void QZGeometryWindow::computeLaplacian( int lapType )
{
	LaplacianType laplacianType = (LaplacianType)lapType;
	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
		mMeshHelper[obj].constructLaplacian(laplacianType);
	});

	int totalToDecompose = 0;
	int nEigVec = gSettings.DEFAULT_EIGEN_SIZE;
	
	for (int obj = 0; obj < mMeshCount; ++obj) {
		if (!mMeshHelper[obj].hasLaplacian(laplacianType))
			throw std::logic_error("Laplacian type not valid!");
        int nEigen = (nEigVec == -1) ? (getMesh(obj)->vertCount() - 1) : nEigVec;
		if (laplacianRequireDecompose(obj, nEigen, laplacianType)) ++totalToDecompose;
	}
	std::cout << totalToDecompose << " mesh Laplacians require explicit decomposition" << std::endl;

	for(int obj = 0; obj < mMeshCount; ++obj) {
        int nEigen = (nEigVec == -1) ? (getMesh(obj)->vertCount() - 1) : nEigVec;
		decomposeSingleLaplacian(obj, nEigen, laplacianType);
	}

	for (int l = 0; l < LaplacianTypeCount; ++l)
		m_actionComputeLaplacians[l]->setChecked(false);
	m_actionComputeLaplacians[laplacianType]->setChecked(true);
}

void QZGeometryWindow::saveMatchingResult()
{
	if (g_task != TASK_REGISTRATION) return;
	const std::vector<MatchPair>& vPairs = mShapeMatcher.getInitialMatchedFeaturePairs();
	if (vPairs.empty()) {
		qout.output("No matching result available", OUT_MSGBOX);
		return;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Matching Result to File"),
													"./output/matching.txt",
													tr("Text Files (*.txt *.dat)"));	
	std::ofstream ofs(fileName.toStdString().c_str(), std::ios::trunc);

	ofs << vPairs.size() << std::endl;
	for (auto iter = vPairs.begin(); iter != vPairs.end(); ++iter) {
		ofs << iter->m_idx1 << ' ' << iter->m_idx2 << ' ' << iter->m_score << std::endl;
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
	//shapeMatcher.generateExampleMatching(20);
	//if (!ui.glMeshWidget->m_bDrawMatching)
	//	toggleDrawMatching();

	ui.glMeshWidget->update();
}

bool QZGeometryWindow::laplacianRequireDecompose( int obj, int nEigVec, LaplacianType laplacianType )
{
	const MeshHelper& mp = mMeshHelper[obj];
    CMesh& mesh = *getMesh(obj);
	
	if (!mp.getMHB(laplacianType).empty()) return false; // already decomposed     
	if (!gSettings.LOAD_MHB_CACHE) return true;    

	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "cache/" + mp.getMesh()->getMeshName() + ".mhb." + s_idx;
	if (!fileExist(pathMHB)) return true;

	std::ifstream ifs(pathMHB.c_str(), std::ios::binary);
	int nEig, nSize;
	ifs.read((char*)&nEig, sizeof(int));
	ifs.read((char*)&nSize, sizeof(int));
	ifs.close();

	if (nEig != nEigVec || nSize != mesh.vertCount()) return true;

	return false;
}

void QZGeometryWindow::allocateStorage( int newMeshCount )
{
	int existingMeshCount = (int)mMeshHelper.size();
	assert(newMeshCount > existingMeshCount);
	for (int k = 0; k < newMeshCount - existingMeshCount; ++k) {
        mMeshHelper.emplace_back();
        mRenderManagers.emplace_back();
	}
	mMeshCount = newMeshCount;
}

void QZGeometryWindow::computeFunctionMaps( int num )
{
	ZGeom::DenseMatrixd funcMap1(num, num), funcMap2(num, num);
	const MeshLaplacian &lap1 = mMeshHelper[0].getMeshLaplacian(CotFormula);
	const MeshLaplacian &lap2 = mMeshHelper[1].getMeshLaplacian(CotFormula);
	const ZGeom::EigenSystem& mhb1 = mMeshHelper[0].getMHB(CotFormula);
	const ZGeom::EigenSystem& mhb2 = mMeshHelper[1].getMHB(CotFormula);
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

void QZGeometryWindow::revertCoord()
{
    mMeshHelper[0].getMesh()->revertCoordinate();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::addColorSignature( int obj, const std::vector<double>& vVals, const std::string& sigName )
{
	auto iResult = minmax_element(vVals.begin(), vVals.end());
	double sMin = *(iResult.first); 
	double sMax = *(iResult.second);
	std::cout << "-- Signature Min = " << sMin << ", Signature Max = " << sMax << std::endl;

    std::vector<double>& vSig = getMesh(obj)->addAttrVertScalars(StrAttrOriginalSignature).attrValue();
	vSig = vVals;

    getMesh(obj)->addColorSigAttr(sigName, ColorSignature(vVals, mColorMapType));
}

double QZGeometryWindow::parameterFromSlider( double sDefault, double sMin, double sMax, bool verbose /*= false*/ )
{
	int sliderCenter = ui.horizontalSliderParamter->maximum() / 2;
	double para;
	if (mCommonParameter <= sliderCenter) 
		para = std::exp( std::log(sDefault / sMin) * (double(mCommonParameter) / sliderCenter) + std::log(sMin) );
	else 
		para = std::exp( std::log(sMax / sDefault) * (double(mCommonParameter - sliderCenter) / sliderCenter) + std::log(sDefault) ); 
	
	if (verbose) std::cout << "Parameter value: " << para << std::endl;
	return para;
}

void QZGeometryWindow::computeGeodesics()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		MeshHelper& mp = mMeshHelper[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();

		std::vector<double> values(meshSize);
		for (int vIdx = 0; vIdx < meshSize; ++vIdx) {
            values[vIdx] = ZGeom::calGeodesic(*getMesh(obj), refPoint, vIdx);
		}

		addColorSignature(obj, values, StrAttrColorGeodesics);
	}

	displaySignature(StrAttrColorGeodesics.c_str());
	updateMenuDisplaySignature();
	mLastOperation = Compute_Geodesics;
}

void QZGeometryWindow::computeHeatTransfer()
{
    double tMultiplier = 2.0;
	for (int obj = 0; obj < mMeshCount; ++obj) {
		MeshHelper *mp = &mMeshHelper[obj];
        const int vertCount = getMesh(obj)->vertCount();
		const int vSrc = mp->getRefPointIndex();		

        ZGeom::SparseSymMatVecSolver heat_solver;
        computeHeatDiffuseMatrix(*getMesh(0), tMultiplier, heat_solver);
        std::vector<double> vHeat = calHeat(*getMesh(0), vSrc, heat_solver);

		addColorSignature(obj, vHeat, StrAttrColorHeat);
	}

	displaySignature(StrAttrColorHeat.c_str());
	updateMenuDisplaySignature();
	mLastOperation = Compute_Heat;
}

void QZGeometryWindow::diffusionFlow()
{
	double tMultiplier = 0.1;
	mShapeEditor.meanCurvatureFlow(tMultiplier, 1);

	std::cout << "Diffusion flow done!" << std::endl;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::addNoise()
{
	double phi = 0.01;
    MeshCoordinates noisyCoord = mShapeEditor.getNoisyCoord(phi);
    mMeshHelper[0].getMesh()->addNamedCoordinate(noisyCoord, "noisy");
	std::cout << "Add Gauss noise with phi=" << phi << std::endl;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::updateSignature( ZGeom::SignatureMode smode )
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		const std::string& currentSig = mRenderManagers[obj].mActiveColorSignatureName;
        if (!getMesh(obj)->hasAttr(currentSig)) continue;
        ColorSignature colors = getMesh(obj)->getColorSignature(currentSig);
        colors.changeColorMap(mColorMapType);
        colors.changeSignatureMode(smode);
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeEditBasis()
{
	if (mShapeEditor.mEditBasis.empty()) {
		std::cout << "Edit basis unavailable!" << std::endl;
		return;
	}

	ui.sliderEditBasis->setMaximum(mShapeEditor.mEditBasis.size()-1);
	ui.sliderEditBasis->setValue(0);
	displayBasis(0);
}

void QZGeometryWindow::clearHandles()
{
	for (int obj = 0; obj < mMeshCount; ++obj) 
		mMeshHelper[obj].clearAllHandles();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToNextCoordinate()
{
    QString coord_name = QString(mMeshHelper[0].getMesh()->switchCoordinate().c_str());
    qout.outputStatus("coord: " + coord_name);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToPrevCoordinate()
{
    QString coord_name = QString(mMeshHelper[0].getMesh()->switchPrevCoordinate().c_str());
    qout.outputStatus("coord: " + coord_name);
    ui.glMeshWidget->update();
}



void QZGeometryWindow::deformLaplace2()
{
	mShapeEditor.deformLaplacian_v2();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::runTests()
{
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setFeaturePointSize( int v )
{
	double scale = std::pow(1.1, double(v-10));
	ui.glMeshWidget->zoomPointSize(scale);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::continuousApprox1( int level )
{
	qout.output("#Reconstruct Basis: " + boost::lexical_cast<std::string>(level), OUT_STATUS);
	mShapeEditor.continuousReconstruct(0, level-1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::continuousApprox2( int level )
{
	qout.output("#Reconstruct Basis: " + boost::lexical_cast<std::string>(level), OUT_STATUS);
	mShapeEditor.continuousReconstruct(1, level-1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::continuousApprox3( int level )
{
	qout.output("#Reconstruct Basis: " + boost::lexical_cast<std::string>(level), OUT_STATUS);
	mShapeEditor.continuousReconstruct(2, level-1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::continuousApprox4( int level )
{
	qout.output("#Reconstruct Basis: " + boost::lexical_cast<std::string>(level), OUT_STATUS);
	mShapeEditor.continuousReconstruct(3, level-1);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayBasis( int idx )
{
	if (mShapeEditor.mEditBasis.empty()) return;
	int select_basis = idx;

	for (int i = 0; i < 1; ++i) {
		MeshHelper& mp = mMeshHelper[i];
		std::vector<double> eigVec = mShapeEditor.mEditBasis[idx].toStdVector();
		addColorSignature(i, eigVec, StrAttrColorWaveletBasis);
	}

	displaySignature(StrAttrColorWaveletBasis.c_str());
	mLastOperation = Compute_Edit_Basis;
	qout.output("Show basis #" + Int2String(select_basis), OUT_STATUS);
	updateMenuDisplaySignature();

	if (mSignatureMode == ZGeom::SM_BandCurved) {
		ui.sliderSigMin->triggerAction(QAbstractSlider::SliderToMinimum);
		ui.sliderSigMax->triggerAction(QAbstractSlider::SliderToMaximum);
		updateSignatureMin(ui.sliderSigMin->minimum());
		updateSignatureMax(ui.sliderSigMax->maximum());
	}
}

void QZGeometryWindow::setSignatureMode( const QString& sigModeName )
{
    for (int i = 0; i < ZGeom::SM_CountSigModes; ++i)
        if (sigModeName == StrSignatureModeNames[i].c_str()) mSignatureMode = (ZGeom::SignatureMode)i;
	
    if (mSignatureMode == ZGeom::SM_BandCurved) {
		ui.sliderSigMin->setEnabled(true);
		ui.sliderSigMax->setEnabled(true);
		ui.sliderSigMin->setValue(ui.sliderSigMin->minimum());
		ui.sliderSigMax->setValue(ui.sliderSigMax->maximum());
		updateSignatureMin(ui.sliderSigMin->minimum());
	} else {
		ui.sliderSigMin->setEnabled(false);
		ui.sliderSigMax->setEnabled(false);
	}

	updateSignature(mSignatureMode);
}

void QZGeometryWindow::updateSignatureMin( int sMin )
{
    if (mSignatureMode != ZGeom::SM_BandCurved) return;
    if (!getMesh(0)->hasAttr(StrAttrOriginalSignature)) return;
	if (sMin >= ui.sliderSigMax->value()) return;

    std::vector<double>& vSig = getMesh(0)->getAttrValue<std::vector<double>>(StrAttrOriginalSignature);
	auto mmp = std::minmax_element(vSig.begin(), vSig.end());
	double vMin = *mmp.first, vMax = *mmp.second;

	double newVal = (double)sMin / (double)ui.sliderSigMin->maximum() * (vMax - vMin) + vMin;
	ui.labelSigMin->setText("Min: " + QString::number(newVal));
	
    updateSignature(ZGeom::SM_BandCurved);
}

void QZGeometryWindow::updateSignatureMax( int sMax )
{
    if (mSignatureMode != ZGeom::SM_BandCurved) return;
    if (!getMesh(0)->hasAttr(StrAttrOriginalSignature)) return;
	if (sMax <= ui.sliderSigMin->value()) return;

    std::vector<double>& vSig = getMesh(0)->getAttrValue<std::vector<double>>(StrAttrOriginalSignature);
	auto mmp = std::minmax_element(vSig.begin(), vSig.end());
	double vMin = *mmp.first, vMax = *mmp.second;

	double newVal = (double)sMax / (double)ui.sliderSigMax->maximum() * (vMax - vMin) + vMin;
	ui.labelSigMax->setText("Max: " + QString::number(newVal));

    updateSignature(ZGeom::SM_BandCurved);
}

void QZGeometryWindow::setLaplacianType( const QString& laplacianTypeName )
{
	qout.output("Select " + laplacianTypeName, OUT_STATUS);

	if (laplacianTypeName == "Umbrella") mActiveLalacian = Umbrella;
	else if (laplacianTypeName == "CotFormula") mActiveLalacian = CotFormula;
	else if (laplacianTypeName == "Anisotropic1") mActiveLalacian = Anisotropic1;
	else if (laplacianTypeName == "Anisotropic2") mActiveLalacian = Anisotropic2;
	else qout.output("Invalid Laplacian type selected!", OUT_MSGBOX);
}

void QZGeometryWindow::captureGL()
{
    QImage img = ui.glMeshWidget->getScreenShot();
	QString filename = "output/screenshots/" + QDateTime::currentDateTime().toString("MM-dd-yyyy_hh.mm.ss") + ".png";
	
	if (img.save(filename))
		qout.output("Screenshot saved to " +  filename + "\n", OUT_CONSOLE);
}

void QZGeometryWindow::captureGLAs()
{
    QImage img = ui.glMeshWidget->getScreenShot();
	QString defaultFilename = "output/screenshots/" + QDateTime::currentDateTime().toString("MM-dd-yyyy_hh.mm.ss") + ".png";
	QString filename = QFileDialog::getSaveFileName(this, tr("Save screenshot"),
													defaultFilename,
													tr("Images (*.png *.jpg)"));

	if (img.save(filename))
	 	qout.output("Screenshot saved to " +  filename + "\n", OUT_CONSOLE);
}

void QZGeometryWindow::openOutputLocation()
{
	ShellExecute(NULL, L"explore", L".\\output", NULL, NULL, SW_RESTORE);
}

void QZGeometryWindow::openSreenshotLocation()
{
	ShellExecute(NULL, L"explore", L".\\output\\screenshots", NULL, NULL, SW_RESTORE);
}

void QZGeometryWindow::resizeApproxSlider( int slider, int newSize )
{
	if (slider == 0) ui.sliderApprox1->setMaximum(newSize);
	else if (slider == 1) ui.sliderApprox2->setMaximum(newSize);
	else if (slider == 2) ui.sliderApprox3->setMaximum(newSize);
	else if (slider == 3) ui.sliderApprox4->setMaximum(newSize);
}

void QZGeometryWindow::visualizeCompression( int selectedApprox, int coordIdx )
{
	mSelectedApprox = selectedApprox;
	mCoordIdx = coordIdx;

    const int vertCount = getMesh(0)->vertCount();
	std::vector<double> vDiff;
	vDiff.resize(vertCount);

	const MeshCoordinates &oldCoord = mShapeEditor.getOldMeshCoord(),
		                  &newCoord = mShapeEditor.getApproximateCoordinate(selectedApprox, coordIdx);

	for (int i = 0; i < vertCount; ++i) {
		vDiff[i] = (oldCoord[i] - newCoord[i]).length();
	}

    getMesh(0)->addColorSigAttr(StrAttrColorPosDiff, ColorSignature(vDiff));

	displaySignature(StrAttrColorPosDiff.c_str());
	updateMenuDisplaySignature();
}

bool QZGeometryWindow::isMeshSelected( int obj )
{
	if (mObjInFocus == -1) return true;
	else return mObjInFocus == obj;
}

void QZGeometryWindow::listMeshAttributes()
{
    std::vector<std::string> attrList = getMesh(0)->getAttrNamesList();
	std::ostringstream ostr;
    CMesh* mesh = getMesh(0);
	ostr << "\nList all attributes:";
    for (size_t i = 0; i < attrList.size(); ++i) {
        std::string attr_name = attrList[i];
        ostr << "\n " << i << "." << attr_name;        
        AttrType attr_type = mesh->getAttrType(attr_name);
        
        if (attr_type == AT_STRING) {
            ostr << ": " << mesh->getAttrValue < std::string>(attr_name);
        }
        else if (attr_type == AT_INT) {
            ostr << ": " << mesh->getAttrValue<int>(attr_name);
        }
        else if (attr_type == AT_DBL) {
            ostr << ": " << mesh->getAttrValue<double>(attr_name);
        }
        else if (attr_type == AT_VEC3) {
            ostr << ": " << mesh->getAttrValue<ZGeom::Vec3d>(attr_name);
        }
    }
    ostr << std::flush;
	qout.output(ostr.str(), OUT_CONSOLE);
}

void QZGeometryWindow::computeVertNormals()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh* mesh = getMesh(obj);
		int vertCount = mesh->vertCount();
        ZGeom::calMeshAttrVertNormals(*mesh);
        auto vNormals = ZGeom::getMeshVertNormals(*mesh);
		MeshLineList mvl;
		for (int i = 0; i < vertCount; ++i)	{
			const ZGeom::Vec3d& vi = mesh->vertPos(i);            
			mvl.push_back(LineSegment(vi, vNormals[i], true));
		}
		
        mesh->addAttrLines(mvl, StrAttrLineVertNormal);
	}

    updateMenuDisplayLines();
    displayLine(StrAttrLineVertNormal.c_str());
}

void QZGeometryWindow::computeFaceNormals()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh* mesh = getMesh(obj);
		int faceCount = mesh->faceCount();
        mesh->calAttrFaceNormals();
		auto fNormals = mesh->getFaceNormals();
		MeshLineList mvl;
		for (int fIdx = 0; fIdx < faceCount; ++fIdx)	{
			ZGeom::Vec3d vc = mesh->getFace(fIdx)->calBarycenter();
			mvl.push_back(LineSegment(vc, fNormals[fIdx], true));
		}

        mesh->addAttrLines(mvl, StrAttrLineFaceNormal);
	}

    updateMenuDisplayLines();
    displayLine(StrAttrLineFaceNormal.c_str());
}

void QZGeometryWindow::fillHoles()
{
    bool skipExternal = false;
    mShapeEditor.fillHoles(skipExternal);

    MeshRegion hole;
    std::set<int> vIn, vBoundary;
    for (const ZGeom::MeshRegion& bv : mShapeEditor.filled_boundaries) {
        for (int vi : bv.vert_inside) vIn.insert(vi);
        for (int vi : bv.vert_on_boundary) vBoundary.insert(vi);
    }
    hole.vert_on_boundary = std::vector < int > {vBoundary.begin(), vBoundary.end()};
    hole.vert_inside = std::vector < int > {vIn.begin(), vIn.end()};
    std::set<int> fHole;
    for (int vi : hole.vert_inside) {
        auto nf = getMesh(0)->vert(vi)->getAdjacentFaces();
        for (const CFace* f : nf) fHole.insert(f->getFaceIndex());
    }
    hole.face_inside = std::vector < int > {fHole.begin(), fHole.end()};

    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_inside, ZGeom::ColorGreen), "hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_on_boundary, ZGeom::ColorRed), "hole_boundary_verts");
    updateMenuDisplayFeatures();
    displayFeature("hole_vertex");
    displayFeature("hole_boundary_verts");


    getMesh(0)->addAttr<vector<int>>(hole.face_inside, "hole_faces", AR_UNIFORM, AT_VEC_INT);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::holeFairingAll()
{
    mShapeEditor.holeFairingFourierOMP();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::holeFairingLS()
{
    mShapeEditor.holeFairingLeastSquare();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::holeFairingFourierOMP()
{
    mShapeEditor.holeFairingFourierOMP();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::holeFairingLARS()
{
    mShapeEditor.holeFairingFourierLARS();
    ui.glMeshWidget->update();
}


void QZGeometryWindow::holeEstimateCurvature()
{
    mShapeEditor.holeEstimateCurvature();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::holeEstimateNormals()
{
    mShapeEditor.holeEstimateNormals();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::generateHoles()
{
    int refIdx = mMeshHelper[0].getRefPointIndex();
    int holeVertCount = 30;

    bool ok;
    int i = QInputDialog::getInt(this, tr("Input hole size"),
        tr("Hole size:"), 25, 1, 10000, 1, &ok);
    if (ok) holeVertCount = i; 
    else return;
    
    MeshRegion generated_holes = ZGeom::generateRandomMeshHole(*getMesh(0), vector<int>{refIdx}, holeVertCount);
    getMesh(0)->addAttr<vector<MeshRegion>>(vector < MeshRegion > {generated_holes}, StrAttrManualHoles, AR_UNIFORM);    
    
    MeshRegion &hole = generated_holes;
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_inside, ZGeom::ColorGreen), "mesh_hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_on_boundary, ZGeom::ColorRed), "mesh_hole_boundary_verts");
    updateMenuDisplayFeatures();

    ui.glMeshWidget->update();
}

void QZGeometryWindow::generateRingHoles()
{
    int refIdx = mMeshHelper[0].getRefPointIndex();
    int ring = 5;

    bool ok;
    int i = QInputDialog::getInt(this, tr("Input hole rings"),
        tr("ring:"), ring, 1, 100, 1, &ok);
    if (ok) ring = i;
    else return;

    MeshRegion generated_holes = ZGeom::generateMeshRingHole(*getMesh(0), refIdx, ring);
    getMesh(0)->addAttr<vector<MeshRegion>>(vector < MeshRegion > {generated_holes}, StrAttrManualHoles, AR_UNIFORM);

    MeshRegion &hole = generated_holes;
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_inside, ZGeom::ColorGreen), "mesh_hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_on_boundary, ZGeom::ColorRed), "mesh_hole_boundary_verts");
    updateMenuDisplayFeatures();

    ui.glMeshWidget->update();
}

void QZGeometryWindow::autoGenerateHoles()
{
    bool ok;
    int N = getMesh(0)->vertCount();

    double missing_ratio = 0.2;
    double s = QInputDialog::getDouble(this, tr("Missing vertex ratio"),
        tr("missing_ratio:"), 0.2, 0.01, 0.75, 2, &ok);
    if (ok) missing_ratio = s;
    int holeVertCount = std::round(missing_ratio * (double)getMesh(0)->vertCount());

    int hole_count = 1;
    int i = QInputDialog::getInt(this, tr("Input number of holes"),
        tr("Hole size:"), 1, 0, 100, 1, &ok);
    if (ok) hole_count = i;
    if (hole_count == 0) hole_count = holeVertCount;


    std::vector<int> seedVerts(N);
    for (int i = 0; i < N; ++i) seedVerts[i] = i;
    std::random_shuffle(seedVerts.begin(), seedVerts.end());
    seedVerts = vector<int>{seedVerts.begin(), seedVerts.begin() + hole_count};

    MeshRegion generated_holes = ZGeom::generateRandomMeshHole(*getMesh(0), seedVerts, holeVertCount);
    getMesh(0)->addAttr<vector<MeshRegion>>(vector<MeshRegion>{generated_holes}, StrAttrManualHoles, AR_UNIFORM);
    
    MeshRegion &hole = generated_holes;
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_inside, ZGeom::ColorGreen), "hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_on_boundary, ZGeom::ColorRed), "hole_boundary_verts");
    updateMenuDisplayFeatures();

    ui.glMeshWidget->update();
}

void QZGeometryWindow::degradeHoles()
{
    double sigma = 0.02;
    bool ok;
    double s = QInputDialog::getDouble(this, tr("Input noise sigma"),
        tr("Noise Signal:"), 0.02, 0, 0.1, 3, &ok);
    if (ok) sigma = s;

    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    auto attrHoles = original_mesh->getAttr<vector<MeshRegion>>(StrAttrManualHoles);
    if (attrHoles == nullptr) {
        std::cout << "No holes selected to degrade" << std::endl;
    }
    vector<MeshRegion>& generated_holes = attrHoles->attrValue();
    mShapeEditor.generateNoise(generated_holes[0].vert_inside, sigma);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::inpaintHoles1()
{
    double eps = 1e-4;
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    auto attrHoles = original_mesh->getAttr<vector<MeshRegion>>(StrAttrManualHoles);
    if (attrHoles == nullptr) {
        std::cout << "No holes selected to inpaint" << std::endl;
    }
    vector<MeshRegion>& generated_holes = attrHoles->attrValue();

    mShapeEditor.inpaintHolesLARS(generated_holes[0].vert_inside, eps);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::inpaintHoles2()
{
    double lambda = 1e-3;
    double tol = 1e-3;
    bool ok;
    double s = QInputDialog::getDouble(this, tr("Input lambda"),
        tr("lambda:"), 1e-3, 0, 0.5, 4, &ok);
    if (ok) lambda = s;

    double eps = 1e-4;
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    auto attrHoles = original_mesh->getAttr<vector<MeshRegion>>(StrAttrManualHoles);
    if (attrHoles == nullptr) {
        std::cout << "No holes selected to inpaint" << std::endl;
    }
    vector<MeshRegion>& generated_holes = attrHoles->attrValue();

    mShapeEditor.inpaintHolesL1LS(generated_holes[0].vert_inside, lambda, tol);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::cutHoles()
{
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    if (!original_mesh->hasAttr(StrAttrManualHoles)) {
        std::cout << "No faces selected to cut" << std::endl;
        return;
    }
    vector<MeshRegion>& generated_holes = original_mesh->getAttrValue<vector<MeshRegion>>(StrAttrManualHoles);
    std::unique_ptr<CMesh> newMesh = std::move(ZGeom::cutFromMesh(*mMeshHelper[0].getMesh(), generated_holes[0].getInsideFaceIdx()));   
        
    newMesh->initNamedCoordinates();
    mMeshHelper[0].addMesh(std::move(newMesh), "selected_partial_mesh");

    mShapeEditor.init(mMeshHelper[0]);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::cutToSelected()
{
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    if (!original_mesh->hasAttr(StrAttrManualHoles)) {
        std::cout << "No faces selected to cut" << std::endl;
        return;
    }
    vector<MeshRegion>& generated_holes = original_mesh->getAttrValue<vector<MeshRegion>>(StrAttrManualHoles);
    std::unique_ptr<CMesh> newMesh = std::move(ZGeom::cutMeshTo(*mMeshHelper[0].getMesh(), generated_holes[0].getInsideFaceIdx()));    
    
    newMesh->initNamedCoordinates();
    mMeshHelper[0].addMesh(std::move(newMesh), "hole_cut_mesh");

    mShapeEditor.init(mMeshHelper[0]);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToNextMesh()
{
    mMeshHelper[0].nextMesh();
    mShapeEditor.init(mMeshHelper[0]);
    
    updateUI();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToPreviousMesh()
{
    mMeshHelper[0].prevMesh();
    mShapeEditor.init(mMeshHelper[0]);

    updateUI();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::detectHoles()
{
}


void QZGeometryWindow::triangulateHoles()
{
    CMesh* oldMesh = mMeshHelper[0].getMesh();
    if (ZGeom::getMeshBoundaryLoops(*oldMesh).empty()) {
        std::cout << "No holes found!" << std::endl;
        return;
    }

    std::unique_ptr<CMesh> newMesh(new CMesh(*oldMesh));
    ZGeom::triangulateMeshHoles(*newMesh);
    mMeshHelper[0].addMesh(std::move(newMesh), "hole_triangulated_mesh");

    ui.glMeshWidget->update();
}

void QZGeometryWindow::refineHoles()
{
    CMesh* oldMesh = mMeshHelper[0].getMesh();
    const auto& oldBoundaries = ZGeom::getMeshBoundaryLoops(*oldMesh);
    if (!oldBoundaries.empty()) {
        std::unique_ptr<CMesh> newMesh(new CMesh(*oldMesh));

        double lambda = 0.8;
        bool ok;
        double r = QInputDialog::getDouble(this, tr("Input refine coefficient"),
            tr("lambda:"), lambda, 0.1, 2, 2, &ok);
        if (ok) lambda = r;
        else return;

        ZGeom::refineMeshHoles(*newMesh, lambda);
        mMeshHelper[0].addMesh(std::move(newMesh), "hole_refined_mesh");
        std::cout << "Mesh hole refined!" << std::endl;

        //evaluateCurrentInpainting();
    }
    else if (oldMesh->hasAttr(StrAttrManualHoles)) {
        auto& manualHoles = oldMesh->getAttrValue<vector<MeshRegion>>(StrAttrManualHoles);
        if (!manualHoles.empty()) {
            std::unique_ptr<CMesh> newMesh(new CMesh(*oldMesh));
            newMesh->addAttr<vector<MeshRegion>>(manualHoles, ZGeom::StrAttrMeshHoleRegions, AR_UNIFORM);
            newMesh->removeAttr(StrAttrManualHoles);
            mMeshHelper[0].addMesh(std::move(newMesh), "hole_refined_mesh");
            std::cout << "Manual hole copied to new mesh!" << std::endl;
        }
    }
    else {
        std::cout << "No holes found! Do nothing." << std::endl;
        return;
    }

    ui.glMeshWidget->update();
}

void QZGeometryWindow::evaluateCurrentInpainting()
{
    /* compare inpainting result with the original generated mesh */
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    CMesh* cur_mesh = mMeshHelper[0].getMesh();
    if (original_mesh->hasAttr(StrAttrManualHoles) && cur_mesh->hasAttr(ZGeom::StrAttrMeshHoleRegions)) {
        vector<MeshRegion>& refinedHoles = mMeshHelper[0].getMesh()->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
        vector<MeshRegion>& manualHoles = mMeshHelper[0].getOriginalMesh()->getAttrValue<vector<MeshRegion>>(StrAttrManualHoles);
        if (refinedHoles.size() == 1 && manualHoles.size() == 1) {
            double rmse = distSubMesh(*cur_mesh, refinedHoles[0].getInsideFaceIdx(), *original_mesh, manualHoles[0].getInsideFaceIdx(), ZGeom::RMSE);
            std::cout << "Inpainting Error: " << rmse << std::endl;
        }
    }
}

void QZGeometryWindow::fairHoleLeastSquares()
{
    /* hole fairing LS */
    CMesh& mesh = *mMeshHelper[0].getMesh();
    const ZGeom::MeshRegion& hole_region = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions)[0];
    int anchor_ring = 3;
    double anchor_weight = 1.0;

    /* input parameters */
    bool ok;
    int r = QInputDialog::getInt(this, tr("Input surrounding ring"),
        tr("Ring:"), anchor_ring, 1, 50, 1, &ok);
    if (ok) anchor_ring = r;
    else return;
    double w = QInputDialog::getDouble(this, tr("Input anchor weight"),
        tr("Weight:"), anchor_weight, 0.1, 10000, 2, &ok);
    if (ok) anchor_weight = w;
    else return;

    MeshCoordinates coord_ls = least_square_fairing(mesh, hole_region, anchor_ring, anchor_weight);
    mesh.addNamedCoordinate(coord_ls, "ls_hole_fairing");
    
    evaluateCurrentInpainting();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::updateUI()
{
    updateMenuDisplaySignature();
    updateMenuDisplayFeatures();
    updateMenuDisplayLines();
    displaySignature(CMesh::StrAttrColorSigDefault.c_str());
}
