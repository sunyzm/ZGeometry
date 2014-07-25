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
#include <stdexcept>
#include <ppl.h>
#include <boost/lexical_cast.hpp>
#include <Shellapi.h>
#include <QMessageBox>
#include <QFileDialog>
#include <QTime>
#include <QImage>
#include <ZGeom/util.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/DenseMatrix.h>

using ZGeom::logic_assert;
using ZGeom::runtime_assert;
using ZGeom::MatlabEngineWrapper;


void normalizeSignature(const std::vector<double>& vOriginalSignature, std::vector<double>& vDisplaySignature)
{
	auto iResult = minmax_element(vOriginalSignature.begin(), vOriginalSignature.end());
	double sMin = *(iResult.first); 
	double sMax = *(iResult.second);

	vDisplaySignature.resize(vOriginalSignature.size());
	std::transform(vOriginalSignature.begin(), vOriginalSignature.end(), 
				   vDisplaySignature.begin(),
				   [=](double v){ return (v-sMin)/(sMax-sMin); });
}

void QZGeometryWindow::signaturesToColors(const std::vector<double>& vOriSig, std::vector<ZGeom::Colorf>& vColors, SignatureMode smode/* = SignatureMode::Normalized*/)
{
	size_t vSize = vOriSig.size();
	vColors.resize(vSize);
	std::vector<double> tmpSig = vOriSig;

	if (smode == SM_BandCurved) {
		auto mmp = std::minmax_element(tmpSig.begin(), tmpSig.end());
		double sMin = *mmp.first, sMax = *mmp.second;
		double bandMin = (double)ui.sliderSigMin->value() / ui.sliderSigMin->maximum() * (sMax - sMin) + sMin;
		double bandMax = (double)ui.sliderSigMax->value() / ui.sliderSigMax->maximum() * (sMax - sMin) + sMin;

		for (size_t a = 0; a < vSize; ++a) {
			if (tmpSig[a] < bandMin) tmpSig[a] = bandMin;
			if (tmpSig[a] > bandMax) tmpSig[a] = bandMax;
		}
	}
	
	if (smode == SM_AbsNormalized) {
		for (double& v : tmpSig) v = std::fabs(v);
	}

	if (smode == SM_LogNormalized) {
		double sMin = *(std::minmax_element(tmpSig.begin(), tmpSig.end()).first);
		if (sMin <= 0) return;
		for (double& v : tmpSig) v = std::log(v);		
	}

	if (smode == SM_PosNegPlot)
	{
		auto mmp = std::minmax_element(tmpSig.begin(), tmpSig.end());
		double sMin = *mmp.first, smax = *mmp.second;
		double absMax = std::max(std::fabs(sMin), std::fabs(smax));
		for (int i = 0; i < vSize; ++i) 
			vColors[i].posNegColor(tmpSig[i] / absMax, ZGeom::ColorOrange, ZGeom::ColorAzure);
	} 
	else if (smode == SM_Normalized || smode == SM_AbsNormalized || smode == SM_BandCurved ||
			 smode == SM_LogNormalized || smode == SM_MarkNegNormalized) 
	{
		std::vector<double> normalizedSig;
		normalizeSignature(tmpSig, normalizedSig);	
		for (int i = 0; i < vSize; ++i) {
			vColors[i].falseColor(normalizedSig[i], 1.f, mColorMapType);		
		}
	}	

	if (smode == SM_MarkNegNormalized) {
		for (int i = 0; i < vSize; ++i) {
			if (vOriSig[i] < 0) vColors[i].setAs(ZGeom::ColorBlack);
		}
	}
}

QZGeometryWindow::QZGeometryWindow(QWidget *parent,  Qt::WindowFlags flags) : QMainWindow(parent, flags)
{
	mMeshCount				= 0;
	mObjInFocus				= -1;
	mCommonParameter		= gSettings.PARAMETER_SLIDER_CENTER;
	mLastOperation			= None;
	mDeformType				= DEFORM_Simple;
	mSignatureMode			= SM_Normalized;
	mActiveLalacian			= CotFormula;
	mColorMapType			= ZGeom::CM_JET;
	mDiffMax				= 2.0;
	mCurrentBasisScale		= 0;


	/* setup ui and connections */
	ui.setupUi(this);
	ui.glMeshWidget->setup(&mProcessors, &mRenderManagers, &mShapeMatcher, &mShapeEditor);
	this->makeConnections();
	
	// comboBoxSigMode
	ui.comboBoxSigMode->clear();
	for (int i = 0; i < SM_CountSigModes; ++i)
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
	mStatusLabel1.setParent(ui.statusBar);
	ui.statusBar->addPermanentWidget(&mStatusLabel1);
	qout.setLabel(&mStatusLabel1);
	qout.setConsole(ui.consoleOutput);
	qout.setStatusBar(ui.statusBar);    	
	mStatusLabel1.setText("Editing");

	setDisplayMesh();
	setEditModeMove();
}

QZGeometryWindow::~QZGeometryWindow()
{
	for (CMesh* m : mMeshes) delete m;
	for (DifferentialMeshProcessor* p : mProcessors) delete p;
	for (RenderSettings* rs : mRenderManagers) delete rs;

	for (QAction* a : m_actionDisplaySignatures) delete a;
	for (QAction* a : m_actionComputeLaplacians) delete a;
	for (QAction* a : m_actionDisplayFeatures) delete a;

	delete m_laplacianSignalMapper;
	delete m_signatureSignalMapper;
	delete m_featureSignalMapper;
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
	QObject::connect(ui.actionComputeMeanCurvature, SIGNAL(triggered()), this, SLOT(computeCurvatureMean()));
	QObject::connect(ui.actionComputeGaussCurvature, SIGNAL(triggered()), this, SLOT(computeCurvatureGauss()));
	QObject::connect(ui.actionComputeHK, SIGNAL(triggered()), this, SLOT(computeHK()));
	QObject::connect(ui.actionComputeHKS, SIGNAL(triggered()), this, SLOT(computeHKS()));
	QObject::connect(ui.actionComputeHKSFeatures, SIGNAL(triggered()), this, SLOT(computeHKSFeatures()));
	QObject::connect(ui.actionComputeMHW, SIGNAL(triggered()), this, SLOT(computeMHW()));
	QObject::connect(ui.actionComputeMHWS, SIGNAL(triggered()), this, SLOT(computeMHWS()));
	QObject::connect(ui.actionComputeSGW, SIGNAL(triggered()), this, SLOT(computeSGW()));
	QObject::connect(ui.actionComputeBiharmonic, SIGNAL(triggered()), this, SLOT(computeBiharmonic()));
	QObject::connect(ui.actionComputeGeodesics, SIGNAL(triggered()), this, SLOT(computeGeodesics()));
	QObject::connect(ui.actionComputeHeatTransfer, SIGNAL(triggered()), this, SLOT(computeHeatTransfer()));
	QObject::connect(ui.actionComputeDictionaryAtom, SIGNAL(triggered()), this, SLOT(computeDictAtom()));
	QObject::connect(ui.actionComputeVertNormals, SIGNAL(triggered()), this, SLOT(computeVertNormals()));
	QObject::connect(ui.actionComputeFaceNormals, SIGNAL(triggered()), this, SLOT(computeFaceNormals()));

	////  Edit  ////
	QObject::connect(ui.actionClearHandles, SIGNAL(triggered()), this, SLOT(clearHandles()));
	QObject::connect(ui.actionRevertCoordinate, SIGNAL(triggered()), this, SLOT(revert()));
	QObject::connect(ui.actionNextCoordinate, SIGNAL(triggered()), this, SLOT(nextCoordinate()));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
	QObject::connect(ui.actionAddNoise, SIGNAL(triggered()), this, SLOT(addNoise()));
	QObject::connect(ui.actionReconstructMHB, SIGNAL(triggered()), this, SLOT(reconstructMHB()));
	QObject::connect(ui.actionReconstructSGW, SIGNAL(triggered()), this, SLOT(reconstructSGW()));
	QObject::connect(ui.actionDeformSimple, SIGNAL(triggered()), this, SLOT(deformSimple()));
	QObject::connect(ui.actionDeformLaplace, SIGNAL(triggered()), this, SLOT(deformLaplace()));
	QObject::connect(ui.actionDeformLaplace2, SIGNAL(triggered()), this, SLOT(deformLaplace2()));
	QObject::connect(ui.actionDeformBiLaplace, SIGNAL(triggered()), this, SLOT(deformBiLaplace()));
	QObject::connect(ui.actionDeformMixedLaplace, SIGNAL(triggered()), this, SLOT(deformMixedLaplace()));
	QObject::connect(ui.actionDeformSGW, SIGNAL(triggered()), this, SLOT(deformSGW()));
	QObject::connect(ui.actionDiffusionFlow, SIGNAL(triggered()), this, SLOT(diffusionFlow()));
	QObject::connect(ui.actionRunTests, SIGNAL(triggered()), this, SLOT(runTests()));
	 
	////  Display  ////
	QObject::connect(ui.actionDisplayMesh, SIGNAL(triggered()), this, SLOT(setDisplayMesh()));
	QObject::connect(ui.actionDisplayWireframe, SIGNAL(triggered()), this, SLOT(setDisplayWireframe()));
	QObject::connect(ui.actionDisplayPointCloud, SIGNAL(triggered()), this, SLOT(setDisplayPointCloud()));
	QObject::connect(ui.actionDisplayNeighbors, SIGNAL(triggered()), this, SLOT(displayNeighborVertices()));
	QObject::connect(ui.actionShowFeatures, SIGNAL(triggered(bool)), this, SLOT(toggleShowFeatures(bool)));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered(bool)), this, SLOT(toggleShowRefPoint(bool)));
	QObject::connect(ui.actionShowSignature, SIGNAL(triggered(bool)), this, SLOT(toggleShowSignature(bool)));
	QObject::connect(ui.actionShowWireframeOverlay, SIGNAL(triggered(bool)), this, SLOT(toggleShowWireframeOverlay(bool)));
	QObject::connect(ui.actionShowColorLegend, SIGNAL(triggered(bool)), this, SLOT(toggleShowColorLegend(bool)));
	QObject::connect(ui.actionShowVectors, SIGNAL(triggered(bool)), this, SLOT(toggleShowVectors(bool)));
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
	QObject::connect(ui.actionBuildHierarchy, SIGNAL(triggered()), this, SLOT(buildHierarchy()));
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
	QObject::connect(&mShapeEditor, SIGNAL(signatureComputed(QString)), this, SLOT(displaySignature(QString)));
	QObject::connect(&mShapeEditor, SIGNAL(signatureComputed(QString)), this, SLOT(updateDisplaySignatureMenu()));
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
		nextBasisScale();
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
			else if (mDeformType == DEFORM_SGW)
				deformSGW();
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
		nextCoordinate();
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
		if (mRenderManagers[0]->displayType == RenderSettings::Mesh)
			setDisplayWireframe();
		else if (mRenderManagers[0]->displayType == RenderSettings::Wireframe)
			setDisplayPointCloud();
		else setDisplayMesh();
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

	/* compute and decompose mesh Laplacians */
	//computeLaplacian(Umbrella);
	//computeLaplacian(NormalizedUmbrella);	
	computeLaplacian(CotFormula);
	//computeLaplacian(SymCot);
	//computeLaplacian(Anisotropic1); 	
	//computeLaplacian(Anisotropic2);
	//setLaplacianType("CotFormula");

	if (g_task == TASK_REGISTRATION) {
		registerPreprocess();
	}
	if (g_task == TASK_EDITING) {
		mShapeEditor.init(mProcessors[0]);
		mShapeEditor.runTests();
		if (mMeshes[0]->hasAttr(StrAttrFeatureSparseSGW))
			displayFeature(StrAttrFeatureSparseSGW.c_str());
	}

	return true;
}

void QZGeometryWindow::loadInitialMeshes(const std::string& mesh_list_name)
{
	runtime_assert (fileExist(mesh_list_name),  "Cannot open file mesh list file!");
	std::ifstream meshfiles(mesh_list_name);
	
	std::vector<std::string> vMeshFiles;
	while (!meshfiles.eof()) {
		std::string meshFileName;
		getline(meshfiles, meshFileName);
		if (meshFileName == "") continue;
		if (meshFileName[0] == '#') continue;
		else if (!fileExist(meshFileName))
			throw std::runtime_error("Cannot open file " + meshFileName);		
		vMeshFiles.push_back(meshFileName);
	}
	meshfiles.close();
	if (vMeshFiles.size() < mMeshCount)
		throw std::runtime_error("Not enough meshes in mesh list!");

	allocateStorage(mMeshCount);

	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
		CMesh& mesh = *mMeshes[obj];
		mesh.load(vMeshFiles[obj]);
		mesh.scaleAreaToVertexNum();
		mesh.gatherStatistics();        
	});	

	for (int obj = 0; obj < mMeshCount; ++obj) {
		CMesh& mesh = *mMeshes[obj];
		Vector3D center = mesh.getCenter();
		Vector3D bbox = mesh.getBoundingBox();
		qout.output(QString().sprintf("Load mesh: %s; Size: %d", mesh.getMeshName().c_str(), mesh.vertCount()), OUT_TERMINAL);
		qout.output(QString().sprintf("Center: (%f, %f, %f)\nDimension: (%f, %f, %f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_TERMINAL);	
		
		mProcessors[obj]->init(&mesh);
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

void QZGeometryWindow::addMesh()
{
	QStringList filenames =  QFileDialog::getOpenFileNames(this, "Select one or more mesh files to open",
														   "../../Data/", "Meshes (*.obj *.off *.ply)");
	int cur_obj = mMeshCount;
	allocateStorage(++mMeshCount);

	CStopWatch timer;
	timer.startTimer();
	CMesh& mesh = *mMeshes[cur_obj];
	mesh.load(filenames.begin()->toStdString());
	mesh.scaleAreaToVertexNum();
	mesh.gatherStatistics();
	timer.stopTimer();
	std::cout << "Time to load meshes: " << timer.getElapsedTime() << "s" << std::endl;

	Vector3D center = mesh.getCenter(), bbox = mesh.getBoundingBox();
	qout.output(QString().sprintf("Load mesh: %s; Size: %d", mesh.getMeshName().c_str(), mesh.vertCount()), OUT_CONSOLE);
	qout.output(QString().sprintf("Center: (%f,%f,%f)\nDimension: (%f,%f,%f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z), OUT_CONSOLE);

	mProcessors[cur_obj]->init(&mesh);

	mRenderManagers[cur_obj]->selected = true;
	mRenderManagers[cur_obj]->mesh_color = preset_mesh_colors[cur_obj%2];

	if (cur_obj == 0) {
		ui.glMeshWidget->fieldView(mMeshes[0]->getCenter(), mMeshes[0]->getBoundingBox());
		ui.spinBox1->setMinimum(0);
		ui.spinBox1->setMaximum(mMeshes[0]->vertCount() - 1);
		ui.horizontalSlider1->setMinimum(0);
		ui.horizontalSlider1->setMaximum(mMeshes[0]->vertCount() - 1);
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
	mMeshes[1]->cloneFrom(*mMeshes[0]);
	mMeshes[1]->gatherStatistics();
	mProcessors[1]->init(mMeshes[1]);
	mRenderManagers[1]->mesh_color = preset_mesh_colors[1];

	qout.output(QString().sprintf("Mesh %s constructed! Size: %d", mMeshes[1]->getMeshName().c_str(), mMeshes[1]->vertCount()));
	
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

void QZGeometryWindow::deformSGW()
{
	mShapeEditor.deformSpectralWavelet();
	mDeformType = DEFORM_SGW;
	ui.glMeshWidget->update();
	setEditModeMove();
}

void QZGeometryWindow::selectObject( int index )
{
	QString text = ui.boxObjSelect->itemText(index);	
	for (auto iter = mRenderManagers.begin(); iter != mRenderManagers.end(); ++iter) (*iter)->selected = false;	
	
	if (text == "1") {			 
		if (mRenderManagers.size() >= 1) mRenderManagers[0]->selected = true;
		mObjInFocus = 0;
	}
	else if (text == "2") {
		if (mRenderManagers.size() >= 2) mRenderManagers[1]->selected = true;
		mObjInFocus = 1;
	}
	else if (text == "All") {
		for (auto iter = mRenderManagers.begin(); iter != mRenderManagers.end(); ++iter) (*iter)->selected = true;	 
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
	mProcessors[0]->setRefPointIndex(vn);
	updateReferenceMove(0);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::setMesh2RefPoint( int vn )
{
	if (mMeshCount < 2) return;

	mProcessors[1]->setRefPointIndex(vn);
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

	for ( auto rm : mRenderManagers ) {
		rm->displayType = RenderSettings::PointCloud;
		rm->glPolygonMode = GL_POINT;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayWireframe()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(true);
	ui.actionDisplayMesh->setChecked(false);
	
	for ( auto rm : mRenderManagers ) {
		rm->displayType = RenderSettings::Wireframe;
		rm->glPolygonMode = GL_LINE;
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::setDisplayMesh()
{
	ui.actionDisplayPointCloud->setChecked(false);
	ui.actionDisplayWireframe->setChecked(false);
	ui.actionDisplayMesh->setChecked(true);

	for ( auto iter = begin(mRenderManagers); iter != end(mRenderManagers); ++iter) {
		(*iter)->displayType = RenderSettings::Mesh;
		(*iter)->glPolygonMode = GL_FILL;
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

void QZGeometryWindow::toggleShowVectors( bool show /*= false*/ )
{
	bool bToShow = !ui.glMeshWidget->m_bShowVectors;
	ui.glMeshWidget->m_bShowVectors = bToShow;
	ui.actionShowVectors->setChecked(bToShow);

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

void QZGeometryWindow::computeCurvatureMean()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		std::vector<double> vCurvature = mMeshes[obj]->getMeanCurvature();
		auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
		qout.output(QString().sprintf("Min curvature: %d  Max curvature: %d", *mm.first, *mm.second));
		addColorSignature(obj, vCurvature, StrAttrColorMeanCurvature);
	}

	updateDisplaySignatureMenu();
	displaySignature(StrAttrColorMeanCurvature.c_str());
	qout.output("Visualize mean curvature");	
}

void QZGeometryWindow::computeCurvatureGauss()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {		
		std::vector<double> vCurvature = mMeshes[obj]->getGaussCurvature();
		auto mm = std::minmax_element(vCurvature.begin(), vCurvature.end());
		qout.output(QString().sprintf("Min curvature: %d  Max curvature: %d", *mm.first, *mm.second));
		addColorSignature(obj, vCurvature, StrAttrColorGaussCurvature);
	}

	updateDisplaySignatureMenu();
	displaySignature(StrAttrColorGaussCurvature.c_str());
	qout.output("Visualize Gauss curvature");
}

void QZGeometryWindow::updateReferenceMove( int obj )
{
	DifferentialMeshProcessor& mp = *mProcessors[obj]; 

	double unitMove = (mp.getMesh_const()->getBoundingBox().x + mp.getMesh_const()->getBoundingBox().y + mp.getMesh_const()->getBoundingBox().z)/300.0;
	Vector3D originalPos = mp.getMesh_const()->getVertex(mp.getRefPointIndex())->getPosition();
	
	mp.setRefPointPosition(originalPos.x, originalPos.y, originalPos.z);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::registerPreprocess()
{
	if (g_task != TASK_REGISTRATION || mMeshCount != 2) return;

	computeFunctionMaps(40);
	mShapeMatcher.initialize(mProcessors[0], mProcessors[1], g_engineWrapper.getEngine());
	std::string rand_data_file = g_configMgr.getConfigValue("RAND_DATA_FILE");
	mShapeMatcher.readInRandPair(rand_data_file);

	// ground truth 
	if (mMeshes[0]->getMeshName() == "march1_1_partial") {
		std::cout << "Ground truth available!" << std::endl;
		std::string mapFile = "./models/map1.txt";
		mShapeMatcher.loadGroundTruth(mapFile);
	}
	else if(mMeshes[0]->vertCount() == mMeshes[1]->vertCount()) {
		mShapeMatcher.autoGroundTruth();
	}

	mShapeMatcher.setRegistrationLevels(1);
	registerTest();
//	evalDistance();
}

void QZGeometryWindow::reconstructMHB()
{
	int sliderCenter = ui.horizontalSliderParamter->maximum() / 2;
	double ratio = std::min((double)mCommonParameter/sliderCenter, 1.0);
	int nEig = mProcessors[0]->getMHB(CotFormula).eigVecCount() * ratio;
	double avgLen = mMeshes[0]->getAvgEdgeLength();

	mShapeEditor.fourierReconstruct(nEig);
	std::cout << "Reconstruct with " << nEig << " eigenvectors" << std::endl;

	std::vector<double> vPosDiff;
	mMeshes[0]->diffCoordinates(mShapeEditor.getOldMeshCoord(), vPosDiff);
	for (double& v : vPosDiff) v /= avgLen;
	std::cout << "Avg Error as ratio of AEL: " << ZGeom::vecMean(vPosDiff) << std::endl;

	addColorSignature(0, vPosDiff, StrAttrColorPosDiff);
	displaySignature(StrAttrColorPosDiff.c_str());
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayNeighborVertices()
{
	int sliderCenter = ui.horizontalSliderParamter->maximum()/2;
	int ring = (mCommonParameter > sliderCenter) ? (mCommonParameter - sliderCenter) : 1;

	int ref = mProcessors[0]->getRefPointIndex();
	std::vector<int> vn;
	mProcessors[0]->getMesh_const()->vertRingNeighborVerts(ref, ring, vn, false);
	MeshFeatureList *mfl = new MeshFeatureList;

	for (auto iter = vn.begin(); iter != vn.end(); ++iter) {
		mfl->addFeature(new MeshFeature(*iter));
	}
	mMeshes[0]->addAttrMeshFeatures(*mfl, StrAttrFeatureNeighbors);

	if (!ui.actionShowFeatures->isChecked()) toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeEigenfunction()
{
	LaplacianType lapType = mActiveLalacian;
	int sliderCenter = ui.horizontalSliderParamter->maximum() / 2;
	int select_eig = (mCommonParameter >= sliderCenter) ? (mCommonParameter - sliderCenter + 1) : 1;
	
	for (int i = 0; i < mMeshCount; ++i) {
		DifferentialMeshProcessor& mp = *mProcessors[i];
		std::vector<double> eigVec = mp.getMHB(lapType).getEigVec(select_eig).toStdVector();
		addColorSignature(i, eigVec, StrAttrColorEigenFunction);
	}

	displaySignature(StrAttrColorEigenFunction.c_str());
	mLastOperation = Compute_Eig_Func;
	qout.output("Show eigenfunction" + Int2String(select_eig), OUT_CONSOLE);
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::computeHK()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int obj = 0; obj < mMeshCount; ++obj) {
		DifferentialMeshProcessor& mp = *mProcessors[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();
		const ZGeom::EigenSystem &mhb = mp.getMHB(mActiveLalacian);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int vIdx) {
			values[vIdx] = mhb.heatKernel(refPoint, vIdx, time_scale);
		});

		addColorSignature(obj, values, StrAttrColorHK);
	}

	displaySignature(StrAttrColorHK.c_str());
	qout.output(QString().sprintf("HK with timescale: %f", time_scale));
	updateDisplaySignatureMenu();
	mLastOperation = Compute_HK;
}

void QZGeometryWindow::computeHKS()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int i = 0; i < mMeshCount; ++i) {
		DifferentialMeshProcessor& mp = *mProcessors[i];
		const int meshSize = mp.getMesh()->vertCount();
		const ZGeom::EigenSystem& mhb = mp.getMHB(mActiveLalacian);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int k) {
			values[k] = mhb.heatKernel(k, k, time_scale);
		});

		addColorSignature(i, values, StrAttrColorHKS);
	}
	
	displaySignature(StrAttrColorHKS.c_str());
	qout.output(QString().sprintf("HKS with timescale: %f", time_scale));
	updateDisplaySignatureMenu();
	mLastOperation = Compute_HKS;
}

void QZGeometryWindow::computeBiharmonic()
{
	for (int obj = 0; obj < mMeshCount; ++obj)
	{
		DifferentialMeshProcessor& mp = *mProcessors[obj];
		const int vertCount = mMeshes[obj]->vertCount();
		const int refPoint = mp.getRefPointIndex();

		std::vector<double> vVals(vertCount);
		for (int vIdx = 0; vIdx < vertCount; ++vIdx) {
			vVals[vIdx] = mp.calBiharmonic(refPoint, vIdx);
		}

		addColorSignature(obj, vVals, StrAttrColorBiharmonic);
	}

	displaySignature(StrAttrColorBiharmonic.c_str());
	updateDisplaySignatureMenu();
	mLastOperation = Compute_Biharmonic;
}

void QZGeometryWindow::computeHKSFeatures()
{
	std::vector<double> vTimes;
	vTimes.push_back(10);
	vTimes.push_back(30);
	vTimes.push_back(90);
	vTimes.push_back(270);

	for (int i = 0; i < mMeshCount; ++i) {
		mProcessors[i]->computeKernelSignatureFeatures(vTimes, HEAT_KERNEL);
	}
	
	if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures();
}

void QZGeometryWindow::computeMHW()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int obj = 0; obj < mMeshCount; ++obj) {
		DifferentialMeshProcessor& mp = *mProcessors[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();
		const ZGeom::EigenSystem& mhb = mp.getMHB(CotFormula);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int vIdx) {
			values[vIdx] = mp.calMHW(refPoint, vIdx, time_scale);
		});

		addColorSignature(obj, values, StrAttrColorMHW);
	}

	displaySignature(StrAttrColorMHW.c_str());
	qout.output(QString().sprintf("MHW timescale: %f", time_scale));
	updateDisplaySignatureMenu();
	mLastOperation = Compute_MHW;
}

void QZGeometryWindow::computeMHWS()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int obj = 0; obj < mMeshCount; ++obj) {
		DifferentialMeshProcessor& mp = *mProcessors[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const ZGeom::EigenSystem& mhb = mp.getMHB(CotFormula);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int vIdx) {
			values[vIdx] = mp.calMHW(vIdx, vIdx, time_scale);
		});

		addColorSignature(obj, values, StrAttrColorMHWS);
	}

	displaySignature(StrAttrColorMHWS.c_str());
	qout.output(QString().sprintf("MHWS timescale: %f", time_scale));
	updateDisplaySignatureMenu();
	mLastOperation = Compute_MHWS;
}


void QZGeometryWindow::computeSGW()
{
	double time_scale = parameterFromSlider(gSettings.DEFUALT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE, true);

	for (int obj = 0; obj < mMeshCount; ++obj) {
		DifferentialMeshProcessor& mp = *mProcessors[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();
		const ZGeom::EigenSystem& mhb = mp.getMHB(mActiveLalacian);

		std::vector<double> values(meshSize);
		Concurrency::parallel_for (0, meshSize, [&](int vIdx) {
			values[vIdx] = mp.calSGW(refPoint, vIdx, time_scale, mActiveLalacian);
		});

		addColorSignature(obj, values, StrAttrColorSGW);
	}

	displaySignature(StrAttrColorSGW.c_str());
	updateDisplaySignatureMenu();
	mLastOperation = Compute_SGW;
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

	case Compute_MHWS:
		computeMHWS(); break;

	case Compute_MHW:
		computeMHW(); break;

	case Compute_SGW:
		computeSGW(); break;

	case Compute_Heat:
		computeHeatTransfer(); break;
	}
}

void QZGeometryWindow::displaySignature(QString sigName )
{
	for (int i = 0; i < mMeshCount; ++i) {
		if (mMeshes[i]->hasAttr(sigName.toStdString()))
			mRenderManagers[i]->mActiveColorSignatureName = sigName.toStdString();
	}

	if (!ui.glMeshWidget->m_bShowSignature && sigName.toStdString() != StrAttrColorPosDiff) 
		toggleShowSignature();	
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayFeature( QString featureName )
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		if (!isMeshSelected(obj)) continue;
		mRenderManagers[obj]->mActiveFeatureName = featureName.toStdString();
	}

	if (!ui.glMeshWidget->m_bShowFeatures) 
		toggleShowFeatures(true);
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayDiffPosition()
{
	runtime_assert(mMeshCount >= 2 && mMeshes[0]->vertCount() == mMeshes[1]->vertCount());
	int size = mMeshes[0]->vertCount();
	std::vector<double> vDiff;
	vDiff.resize(size);

	for (int i = 0; i < mMeshes[0]->vertCount(); ++i) {
		vDiff[i] = (mMeshes[0]->getVertex(i)->getPosition() - mMeshes[1]->getVertex(i)->getPosition()).length() / mMeshes[0]->getAvgEdgeLength();
	}

	addColorSignature(0, vDiff, StrAttrColorPosDiff);
	displaySignature(StrAttrColorPosDiff.c_str());
	updateDisplaySignatureMenu();
}

void QZGeometryWindow::updateDisplaySignatureMenu()
{
	int obj = (mObjInFocus <= 0 ? 0 : mObjInFocus);
	std::vector<AttrVertColors*> vColorAttributes = mMeshes[obj]->getColorAttrList();

	for (QAction* qa : m_actionDisplaySignatures) {
		if (find_if(vColorAttributes.begin(), vColorAttributes.end(), [&](AttrVertColors* attr){ return attr->attrName() == qa->text().toStdString();}) 
			== vColorAttributes.end())
		{
				ui.menuDisplaySignatures->removeAction(qa);
				delete qa;	
		}
	}

	for (AttrVertColors* attr : vColorAttributes) {
		if (find_if(m_actionDisplaySignatures.begin(), m_actionDisplaySignatures.end(), [&](QAction* pa){ return pa->text().toStdString() == attr->attrName();})
			== m_actionDisplaySignatures.end())
		{
			QAction* newDisplayAction = new QAction(attr->attrName().c_str(), this);
			m_actionDisplaySignatures.push_back(newDisplayAction);
			ui.menuDisplaySignatures->addAction(m_actionDisplaySignatures.back());
			m_signatureSignalMapper->setMapping(newDisplayAction, attr->attrName().c_str());
			QObject::connect(newDisplayAction, SIGNAL(triggered()), m_signatureSignalMapper, SLOT(map()));
		}	
	}
}

void QZGeometryWindow::updateDisplayFeatureMenu()
{
	int obj = (mObjInFocus <= 0 ? 0 : mObjInFocus);
	std::vector<AttrMeshFeatures*> vFeatureAttr = mMeshes[obj]->getMeshFeatureList();

	for (QAction* qa : m_actionDisplayFeatures) {
		if (find_if(vFeatureAttr.begin(), vFeatureAttr.end(), [&](AttrMeshFeatures* attr){ return attr->attrName() == qa->text().toStdString();}) 
			== vFeatureAttr.end())
		{
			ui.menuDisplayFeatures->removeAction(qa);
			delete qa;	
		}
	}

	for (AttrMeshFeatures* attr : vFeatureAttr) {
		if (find_if(m_actionDisplaySignatures.begin(), m_actionDisplaySignatures.end(), [&](QAction* pa){ return pa->text().toStdString() == attr->attrName();})
			== m_actionDisplaySignatures.end())
		{
			QAction* newDisplayAction = new QAction(attr->attrName().c_str(), this);
			m_actionDisplayFeatures.push_back(newDisplayAction);
			ui.menuDisplayFeatures->addAction(m_actionDisplayFeatures.back());
			m_featureSignalMapper->setMapping(newDisplayAction, attr->attrName().c_str());
			QObject::connect(newDisplayAction, SIGNAL(triggered()), m_featureSignalMapper, SLOT(map()));
		}	
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
	this->buildHierarchy();
	this->detectFeatures();
	this->matchFeatures();
}

void QZGeometryWindow::buildHierarchy()
{
	int nLevel = g_configMgr.getConfigValueInt("HIERARCHY_LEVEL");
	double ratio = g_configMgr.getConfigValueDouble("CONTRACTION_RATIO");
	std::string log_filename = g_configMgr.getConfigValue("HIERARCHY_OUTPUT_FILE");
	std::ofstream ostr(log_filename.c_str(), std::ios::trunc);

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
	using namespace std;

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
		matchScore = ShapeMatcher::TensorGraphMatching6(g_engineWrapper.getEngine(), mProcessors[0], mProcessors[1], vftFine1, vftFine2, vPairs, tensor_matching_timescasle, matching_thresh_2, /*verbose=*/true);
		//matchScore = DiffusionShapeMatcher::TensorMatchingExt(m_ep, &vMP[0], &vMP[1], vFeatures1, vFeatures2, vPairs, 0, vPara, cout, true);

		if (1 == g_configMgr.getConfigValueInt("GROUND_TRUTH_AVAILABLE")) {
			vector<double> vTimes;
			vTimes.push_back(20); vTimes.push_back(40); vTimes.push_back(80); vTimes.push_back(160); vTimes.push_back(320);
			for (auto iter = vPairs.begin(); iter != vPairs.end(); ) {
				if (!mMeshes[1]->isInNeighborRing(iter->m_idx1, iter->m_idx2, 2))
					iter->m_note = -1;

				double dissim = mShapeMatcher.calPointHksDissimilarity(mProcessors[0], mProcessors[1], iter->m_idx1, iter->m_idx2, vTimes, 1);
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
	using namespace std;

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

void QZGeometryWindow::evalDistance()
{
	using namespace std;

	if (mMeshCount < 2) return;

	CStopWatch timer;
	timer.startTimer();
	Concurrency::parallel_invoke(
		[&](){ std::cout << "Error geodesic1: " << ShapeMatcher::evaluateDistance(*mProcessors[0], *mProcessors[1], DISTANCE_GEODESIC, std::vector<double>(), mShapeMatcher.m_randPairs, 0) << endl; },
		[&](){ std::cout << "Error biharmonic1: " << ShapeMatcher::evaluateDistance(*mProcessors[0], *mProcessors[1], DISTANCE_BIHARMONIC, std::vector<double>(), mShapeMatcher.m_randPairs, 0) << endl; },
		[&](){ std::cout << "Error biharmonic2: " << ShapeMatcher::evaluateDistance(*mProcessors[0], *mProcessors[1], DISTANCE_BIHARMONIC, std::vector<double>(), mShapeMatcher.m_randPairs, 500) << endl; }
//		[&](){cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 30.), shapeMatcher.m_randPairs, 0) << endl;},
//		cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 90.), shapeMatcher.m_randPairs, 0) << endl;
//		cout << "Error geodesic: " << DiffusionShapeMatcher::evaluateDistance(&vMP[0], &vMP[1], DISTANCE_HK, std::vector<double>(1, 270.), shapeMatcher.m_randPairs, 0) << endl;
	);
//	cout << "Error geodesic1: " << DiffusionShapeMatcher::evaluateDistance(vMP[0], vMP[1], DISTANCE_GEODESIC, std::vector<double>(), shapeMatcher.m_randPairs, 500) << endl;
	timer.stopTimer();
	std::cout << "Evaluate Dist time (ppl): " << timer.getElapsedTime() << endl;

}

void QZGeometryWindow::decomposeSingleLaplacian( int obj, int nEigVec, LaplacianType laplacianType /*= CotFormula*/ )
{
	DifferentialMeshProcessor& mp = *mProcessors[obj];
	const CMesh& mesh = *mMeshes[obj];
	const int vertCount = mesh.vertCount();
	if (nEigVec >= vertCount) nEigVec = vertCount - 1;

	if (!mp.getMHB(laplacianType).empty()) return;
	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "cache/" + mp.getMesh_const()->getMeshName() + ".mhb." + s_idx;
	
	if (gSettings.LOAD_MHB_CACHE && fileExist(pathMHB) && laplacianType != Anisotropic1 )	// MHB cache available for the current mesh
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

	std::cout << "Min EigVal: " << mp.getMHB(laplacianType).getEigVals().front() 
			  << "; Max EigVal: " << mp.getMHB(laplacianType).getEigVals().back() << std::endl;
}

void QZGeometryWindow::saveSignature()
{
	if (!mMeshes[0]->hasAttr(StrAttrOriginalSignature)) {
		qout.output("No signature available", OUT_MSGBOX);
		return;
	}

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Signature to File"),
		"./output/signature.txt",
		tr("Text Files (*.txt *.dat)"));
	const std::vector<double>& vSig = mMeshes[0]->getAttrValue<std::vector<double>,AR_VERTEX>(StrAttrOriginalSignature);

	vector2file<double>(fileName.toStdString(), vSig);
}

void QZGeometryWindow::computeLaplacian( int lapType )
{
	LaplacianType laplacianType = (LaplacianType)lapType;
	Concurrency::parallel_for(0, mMeshCount, [&](int obj) {
		mProcessors[obj]->constructLaplacian(laplacianType);
	});

	int totalToDecompose = 0;
	int nEigVec = gSettings.DEFAULT_EIGEN_SIZE;

	for (int obj = 0; obj < mMeshCount; ++obj) {
		if (!mProcessors[obj]->hasLaplacian(laplacianType))
			throw std::logic_error("Laplacian type not valid!");
		if (laplacianRequireDecompose(obj, nEigVec, laplacianType)) 
			++totalToDecompose;
	}
	std::cout << totalToDecompose << " mesh Laplacians require explicit decomposition" << std::endl;

	if (totalToDecompose <= 1 && mMeshCount > 1) {        
		Concurrency::parallel_for(0, mMeshCount, [&](int obj){			
			decomposeSingleLaplacian(obj, nEigVec, laplacianType);    
		});
	} else {    // if both need explicit decomposition, then must run decomposition in sequence in Matlab
		for(int obj = 0; obj < mMeshCount; ++obj) {
			decomposeSingleLaplacian(obj, nEigVec, laplacianType);    
		}
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

bool QZGeometryWindow::laplacianRequireDecompose( int obj, int nEigVec, LaplacianType laplacianType ) const
{
	const DifferentialMeshProcessor& mp = *mProcessors[obj];
	const CMesh& mesh = *mMeshes[obj];
	
	if (!mp.getMHB(laplacianType).empty()) return false; // already decomposed     
	if (!gSettings.LOAD_MHB_CACHE) return true;    

	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = "cache/" + mp.getMesh_const()->getMeshName() + ".mhb." + s_idx;
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
	const MeshLaplacian &lap1 = mProcessors[0]->getMeshLaplacian(CotFormula);
	const MeshLaplacian &lap2 = mProcessors[1]->getMeshLaplacian(CotFormula);
	const ZGeom::EigenSystem& mhb1 = mProcessors[0]->getMHB(CotFormula);
	const ZGeom::EigenSystem& mhb2 = mProcessors[1]->getMHB(CotFormula);
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
	for (int obj = 0; obj < mMeshCount; ++obj) {
		double areaSum(0);
		for (int i = 0; i < mMeshes[obj]->faceCount(); ++i) {
			areaSum += mMeshes[obj]->calFaceArea(i);
		}
		double weightSum(0);
		const MeshLaplacian& laplacian = mProcessors[obj]->getMeshLaplacian(CotFormula);
		std::vector<double> vAreas;
		laplacian.getW().getDiagonal(vAreas);
		for (int i = 0; i < mMeshes[obj]->vertCount(); ++i) {
			weightSum += vAreas[i];
		}
		std::cout << "Vert count: " << mMeshes[obj]->vertCount() << std::endl;
		std::cout << "Total surface area: " << areaSum << std::endl;
		std::cout << "Total vert weight: " << weightSum << std::endl;
	}
}

void QZGeometryWindow::revert()
{
	mShapeEditor.revertCoordinates();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::addColorSignature( int obj, const std::vector<double>& vVals, const std::string& sigName )
{
	auto iResult = minmax_element(vVals.begin(), vVals.end());
	double sMin = *(iResult.first); 
	double sMax = *(iResult.second);
	std::cout << "-- Signature Min = " << sMin << ", Signature Max = " << sMax << std::endl;

	std::vector<double>& vSig = mMeshes[obj]->addAttrVertScalars(StrAttrOriginalSignature).attrValue();
	vSig = vVals;

	mMeshes[obj]->addColorAttr(sigName);
	std::vector<ZGeom::Colorf>& vColors = mMeshes[obj]->getVertColors(sigName);
	signaturesToColors(vVals, vColors, mSignatureMode);	

	/* update color legend */
	ui.glMeshWidget->mLegendColors.resize(256);
	std::vector<double> vLegendVals(256);
	for (int i = 0; i <= 255; ++i) vLegendVals[i] = sMin + (double)i/255.*(sMax-sMin);
	signaturesToColors(vLegendVals, ui.glMeshWidget->mLegendColors, mSignatureMode);
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

void QZGeometryWindow::reconstructSGW()
{
	mShapeEditor.reconstructSpectralWavelet();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::computeGeodesics()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		DifferentialMeshProcessor& mp = *mProcessors[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();

		std::vector<double> values(meshSize);
		for (int vIdx = 0; vIdx < meshSize; ++vIdx) {
			values[vIdx] = mMeshes[obj]->calGeodesic(refPoint, vIdx);
		}

		addColorSignature(obj, values, StrAttrColorGeodesics);
	}

	displaySignature(StrAttrColorGeodesics.c_str());
	updateDisplaySignatureMenu();
	mLastOperation = Compute_Geodesics;
}

void QZGeometryWindow::computeHeatTransfer()
{
	double tMultiplier = 1.0;
	for (int obj = 0; obj < mMeshCount; ++obj) {
		DifferentialMeshProcessor *mp = mProcessors[obj];
		const int vertCount = mMeshes[obj]->vertCount();
		const int vSrc = mp->getRefPointIndex();
		std::vector<double> vHeat;

		mp->calHeat(vSrc, tMultiplier, vHeat);
		addColorSignature(obj, vHeat, StrAttrColorHeat);
	}

	displaySignature(StrAttrColorHeat.c_str());
	updateDisplaySignatureMenu();
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
	double phi = 0.1;
	mShapeEditor.addNoise(phi);

	std::cout << "Add Gauss noise with phi=" << phi << std::endl;
	ui.glMeshWidget->update();
}

void QZGeometryWindow::updateSignature( SignatureMode smode )
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		const std::string& currentSig = mRenderManagers[obj]->mActiveColorSignatureName;
		if (!mMeshes[obj]->hasAttr(currentSig) || !mMeshes[obj]->hasAttr(StrAttrOriginalSignature)) continue;

		const std::vector<double>& vSig = mMeshes[obj]->getVertScalars(StrAttrOriginalSignature);
		std::vector<ZGeom::Colorf>& vColors = mMeshes[obj]->getVertColors(currentSig);
		signaturesToColors(vSig, vColors, smode);

		if (obj == 0) {
			auto iResult = minmax_element(vSig.begin(), vSig.end());
			double sMin = *(iResult.first); 
			double sMax = *(iResult.second);
			ui.glMeshWidget->mLegendColors.resize(256);
			std::vector<double> vLegendVals(256);
			for (int i = 0; i <= 255; ++i) vLegendVals[i] = sMin + (double)i/255.*(sMax-sMin);
			signaturesToColors(vLegendVals, ui.glMeshWidget->mLegendColors, mSignatureMode);
		}
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
		mProcessors[obj]->clearAllHandles();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::nextCoordinate()
{
	mShapeEditor.nextCoordinates();
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
		DifferentialMeshProcessor& mp = *mProcessors[i];
		std::vector<double> eigVec = mShapeEditor.mEditBasis[idx].toStdVector();
		addColorSignature(i, eigVec, StrAttrColorWaveletBasis);
	}

	displaySignature(StrAttrColorWaveletBasis.c_str());
	mLastOperation = Compute_Edit_Basis;
	qout.output("Show basis #" + Int2String(select_basis), OUT_STATUS);
	updateDisplaySignatureMenu();

	if (mSignatureMode == SM_BandCurved) {
		ui.sliderSigMin->triggerAction(QAbstractSlider::SliderToMinimum);
		ui.sliderSigMax->triggerAction(QAbstractSlider::SliderToMaximum);
		updateSignatureMin(ui.sliderSigMin->minimum());
		updateSignatureMax(ui.sliderSigMax->maximum());
	}
}

void QZGeometryWindow::setSignatureMode( const QString& sigModeName )
{
	for (int i = 0; i < SM_CountSigModes; ++i)
		if (sigModeName == StrSignatureModeNames[i].c_str()) mSignatureMode = (SignatureMode)i;
	
	if (mSignatureMode == SM_BandCurved) {
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
	if (mSignatureMode != SM_BandCurved) return;
	if (!mMeshes[0]->hasAttr(StrAttrOriginalSignature)) return;
	if (sMin >= ui.sliderSigMax->value()) return;

	std::vector<double>& vSig = mMeshes[0]->getAttrValue<std::vector<double>,AR_VERTEX>(StrAttrOriginalSignature);
	auto mmp = std::minmax_element(vSig.begin(), vSig.end());
	double vMin = *mmp.first, vMax = *mmp.second;

	double newVal = (double)sMin / (double)ui.sliderSigMin->maximum() * (vMax - vMin) + vMin;
	ui.labelSigMin->setText("Min: " + QString::number(newVal));
	
	updateSignature(SM_BandCurved);
}

void QZGeometryWindow::updateSignatureMax( int sMax )
{
	if (mSignatureMode != SM_BandCurved) return;
	if (!mMeshes[0]->hasAttr(StrAttrOriginalSignature)) return;
	if (sMax <= ui.sliderSigMin->value()) return;

	std::vector<double>& vSig = mMeshes[0]->getAttrValue<std::vector<double>,AR_VERTEX>(StrAttrOriginalSignature);
	auto mmp = std::minmax_element(vSig.begin(), vSig.end());
	double vMin = *mmp.first, vMax = *mmp.second;

	double newVal = (double)sMax / (double)ui.sliderSigMax->maximum() * (vMax - vMin) + vMin;
	ui.labelSigMax->setText("Max: " + QString::number(newVal));

	updateSignature(SM_BandCurved);
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
	QImage img = ui.glMeshWidget->grabFrameBuffer();
	QString filename = "output/screenshots/" + QDateTime::currentDateTime().toString("MM-dd-yyyy_hh.mm.ss") + ".png";
	
	if (img.save(filename))
		qout.output("Screenshot saved to " +  filename + "\n", OUT_CONSOLE);
}

void QZGeometryWindow::captureGLAs()
{
	QImage img = ui.glMeshWidget->grabFrameBuffer();
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

void QZGeometryWindow::nextBasisScale()
{
	if (mShapeEditor.mTotalScales <= 0) return;
	mCurrentBasisScale = (mCurrentBasisScale + 1) % mShapeEditor.mTotalScales;
	std::cout << "Current basis scale: " << mCurrentBasisScale << std::endl;

	computeDictAtom();
}

void QZGeometryWindow::computeDictAtom()
{
	DifferentialMeshProcessor& mp = *mProcessors[0];
	const int meshSize = mp.getMesh()->vertCount();
	const int refPoint = mp.getRefPointIndex();

	std::vector<double> values(meshSize);
	values = mp.getWaveletMat().getRowVec(meshSize * mCurrentBasisScale + refPoint).toStdVector();

	addColorSignature(0, values, StrAttrColorDictAtom);	

	displaySignature(StrAttrColorDictAtom.c_str());
	updateDisplaySignatureMenu();
	mLastOperation = Compute_Dict_Atom;
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

	const int vertCount = mMeshes[0]->vertCount();
	std::vector<double> vDiff;
	vDiff.resize(vertCount);

	const MeshCoordinates &oldCoord = mShapeEditor.getOldMeshCoord(),
		                  &newCoord = mShapeEditor.getApproximateCoordinate(selectedApprox, coordIdx);

	for (int i = 0; i < vertCount; ++i) {
		vDiff[i] = (oldCoord[i] - newCoord[i]).length();
	}

	auto iResult = std::minmax_element(vDiff.begin(), vDiff.end());
	double sMin = *(iResult.first); 
	double sMax = *(iResult.second);
	std::cout << "-- Signature Min = " << sMin << ", Signature Max = " << sMax << std::endl;

	mMeshes[0]->addColorAttr(StrAttrColorPosDiff);
	std::vector<ZGeom::Colorf>& vColors = mMeshes[0]->getVertColors(StrAttrColorPosDiff);
	for (int i = 0; i < vertCount; ++i) {
		if (vDiff[i] <= mDiffMax) vColors[i].falseColor(vDiff[i]/mDiffMax, 1.f, mColorMapType);
		else vColors[i].falseColor(1.0f, 1.f, mColorMapType);
	}
	
	/* update color legend */
	ui.glMeshWidget->mLegendColors.resize(256);
	std::vector<double> vLegendVals(256);
	for (int i = 0; i <= 255; ++i) vLegendVals[i] = float(i)/255.0;
	signaturesToColors(vLegendVals, ui.glMeshWidget->mLegendColors, mSignatureMode);

	displaySignature(StrAttrColorPosDiff.c_str());
	updateDisplaySignatureMenu();
}

bool QZGeometryWindow::isMeshSelected( int obj )
{
	if (mObjInFocus == -1) return true;
	else return mObjInFocus == obj;
}

void QZGeometryWindow::listMeshAttributes()
{
	std::vector<std::string> attrList = mMeshes[0]->getAttrNamesList();
	std::ostringstream ostr;
	ostr << "List all attributes:";
	for (size_t i = 0; i < attrList.size(); ++i) 
		ostr << "\n " << i << ": " << attrList[i];
	qout.output(ostr.str(), OUT_CONSOLE);
}

void QZGeometryWindow::computeVertNormals()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		CMesh* mesh = mMeshes[obj];
		int vertCount = mesh->vertCount();
		const std::vector<Vector3D>& vNormals = mesh->getVertNormals();
		assert(vNormals.size() == vertCount);
		
		MeshVectorList mvl;
		for (int i = 0; i < vertCount; ++i)	{
			const Vector3D& vi = mesh->getVertexPosition(i);
			mvl.push_back(std::make_pair(vi, vNormals[i]));
		}
		
		mesh->addAttr<MeshVectorList,AR_UNIFORM>(mvl, StrAttrVecVertNormal, AT_VEC3);
		mRenderManagers[obj]->mActiveVectorName = StrAttrVecVertNormal;
	}

	if (!ui.glMeshWidget->m_bShowVectors) toggleShowVectors(true);
	else ui.glMeshWidget->update();
}

void QZGeometryWindow::computeFaceNormals()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
		CMesh* mesh = mMeshes[obj];
		int faceCount = mesh->faceCount();
		const std::vector<Vector3D>& fNormals = mesh->getFaceNormals();
		assert(fNormals.size() == faceCount);

		MeshVectorList mvl;

		for (int fIdx = 0; fIdx < faceCount; ++fIdx)	{
			Vector3D vc = mesh->getFace(fIdx)->calBarycenter();
			mvl.push_back(std::make_pair(vc, fNormals[fIdx]));
		}

		mesh->addAttr<MeshVectorList,AR_UNIFORM>(mvl, StrAttrVecFaceNormal, AT_VEC3);
		mRenderManagers[obj]->mActiveVectorName = StrAttrVecFaceNormal;
	}

	if (!ui.glMeshWidget->m_bShowVectors) toggleShowVectors(true);
	else ui.glMeshWidget->update();
}
