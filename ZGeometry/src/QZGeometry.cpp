#include "QZGeometry.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
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
#include <QColorDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QTime>
#include <QImage>
#include <ZGeom/geometry_processing.h>
#include <ZGeom/util.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/DenseMatrix.h>
#include "heat_diffusion.h"
#include "hole_fairing.h"

using namespace concurrency;

using std::vector;
using ZGeom::Colorf;
using ZGeom::logic_assert;
using ZGeom::runtime_assert;
using ZGeom::MatlabEngineWrapper;
using ZGeom::ColorSignature;
using ZGeom::MeshRegion;

vector<double> calHoleVertDistToHoleFace(const CMesh& mesh1, const ZGeom::MeshRegion& hole1, const CMesh& mesh2, const ZGeom::MeshRegion& hole2)
{
    using namespace ZGeom;
    int nVert1 = (int)hole1.vert_inside.size();
    int nFace2 = (int)hole2.face_inside.size();

    vector<vector<Vec3d>> holeFace2Tri(nFace2);
    for (int k = 0; k < nFace2; ++k)
        holeFace2Tri[k] = mesh2.getFace(hole2.face_inside[k])->getAllVertPos();

    vector<double> result(nVert1);
    concurrency::parallel_for(0, nVert1, [&](int k)
    {
        const Vec3d &vPos = mesh1.vertPos(hole1.vert_inside[k]);
        double minDistVi = 1e15;
        for (const auto& tri : holeFace2Tri)
            minDistVi = std::min(minDistVi, distPointTriangle(vPos, tri).distance);
        result[k] = minDistVi;
    });

    return result;
}

double rmseVerts(const std::vector<int>& verts, const std::vector<double>& all_vert_areas, const std::vector<double>& all_vert_error)
{
    double area_sum(0);
    double s2_error_sum(0);
    for (int vi : verts) {
        s2_error_sum += all_vert_areas[vi] * all_vert_error[vi] * all_vert_error[vi];
        area_sum += all_vert_areas[vi];
    }

    return std::sqrt(s2_error_sum / area_sum);
}

bool getfairHoleL1LsParameters(QWidget* parent, ParaL1LsInpainting& para)
{
    using namespace std;

    QDialog dialog(parent);

    // Use a layout allowing to have a label next to each field
    QFormLayout form(&dialog);

    unique_ptr<QSpinBox> ring_input = make_unique<QSpinBox>(&dialog);
    ring_input->setRange(0, 100);
    ring_input->setValue(para.fitting_ring);
    
    unique_ptr<QSpinBox> eigen_input = make_unique<QSpinBox>(&dialog);
    eigen_input->setRange(-1, 10000);
    eigen_input->setValue(para.eigen_count);

    unique_ptr<QDoubleSpinBox> tol_input = make_unique<QDoubleSpinBox>(&dialog);
    tol_input->setDecimals(4);
    tol_input->setRange(1e-4, 1e-2);
    tol_input->setSingleStep(1e-4);
    tol_input->setValue(para.tol);

    unique_ptr<QDoubleSpinBox> lambda_input = make_unique<QDoubleSpinBox>(&dialog);
    lambda_input->setDecimals(4);
    lambda_input->setRange(1e-4, 0.9999);
    lambda_input->setSingleStep(1e-3);
    lambda_input->setValue(para.lambda);

    form.addRow("Surrounding ring", ring_input.get());
    form.addRow("#Eigenfunctions", eigen_input.get());
    form.addRow("tol", tol_input.get());
    form.addRow("lambda", lambda_input.get());

    // Add some standard buttons (Cancel/Ok) at the bottom of the dialog
    QDialogButtonBox buttonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
        Qt::Horizontal, &dialog);
    form.addRow(&buttonBox);
    QObject::connect(&buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
    QObject::connect(&buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

    // Show the dialog as modal
    if (dialog.exec() == QDialog::Accepted) {
        para.fitting_ring = ring_input->value();
        para.eigen_count = eigen_input->value();
        para.tol = tol_input->value();
        para.lambda = lambda_input->value();
        return true;
    }
    else return false;
}

QZGeometryWindow::QZGeometryWindow(QWidget *parent,  Qt::WindowFlags flags) : QMainWindow(parent, flags)
{
	mMeshCount				= 0;
	mObjInFocus				= -1;
	mCommonParameter		= gSettings.PARAMETER_SLIDER_CENTER;
	mLastOperation			= None;
	mDeformType				= DEFORM_Simple;
	active_lap_type			= Umbrella;
	mDiffMax				= 2.0;
	mCurrentBasisScale		= 0;

	/* setup ui and connections */
	ui.setupUi(this);
	this->makeConnections();
	
	// comboBoxLaplacian
	ui.comboBoxLaplacian->clear();
    for (const std::string& lap_type_str : StrLaplacianTypes) {
        ui.comboBoxLaplacian->addItem(lap_type_str.c_str());
    }
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
	for (QAction* a : m_actionDisplayFeatures) delete a;
    for (QAction* a : m_actionDisplayLines) delete a;

	delete m_signatureSignalMapper;
	delete m_featureSignalMapper;
    delete m_linesSignalMapper;
}

void QZGeometryWindow::makeConnections()
{
	QObject::connect(ui.actionAboutQt, SIGNAL(triggered()), this, SIGNAL(displayQtVersion()));

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
    QObject::connect(ui.actionListAttributes, SIGNAL(triggered()), this, SLOT(listMeshAttributes()));
    QObject::connect(ui.actionRevertCoordinate, SIGNAL(triggered()), this, SLOT(revertCoord()));
    QObject::connect(ui.actionNextCoordinate, SIGNAL(triggered()), this, SLOT(switchToNextCoordinate()));
    QObject::connect(ui.actionSwitchMesh, SIGNAL(triggered()), this, SLOT(switchToNextMesh()));
    QObject::connect(ui.actionRemoveCurrentMesh, SIGNAL(triggered()), this, SLOT(removeCurrentMesh()));
    QObject::connect(ui.actionDeleteCurrentCoordinate, SIGNAL(triggered()), this, SLOT(deleteCurrentCoordinate()));

	////  Compute  ////
    QObject::connect(ui.actionComputeGraphLaplacian, SIGNAL(triggered()), this, SLOT(computeGraphLaplacian()));
    QObject::connect(ui.actionComputeGeoLaplacian, SIGNAL(triggered()), this, SLOT(computeGeoLaplacian()));
	QObject::connect(ui.actionEigenfunction, SIGNAL(triggered()), this, SLOT(computeEigenfunction()));
    QObject::connect(ui.actionComputeCurvatures, SIGNAL(triggered()), this, SLOT(computeCurvatures()));
    QObject::connect(ui.actionComputeShapeIndex, SIGNAL(triggered()), this, SLOT(computeShapeIndex()));
    QObject::connect(ui.actionComputeBiharmonicDistance, SIGNAL(triggered()), this, SLOT(computeBiharmonicDistField()));
    QObject::connect(ui.actionComputeSGW, SIGNAL(triggered()), this, SLOT(computeSGW()));
    QObject::connect(ui.actionComputeHK, SIGNAL(triggered()), this, SLOT(computeHK()));
	QObject::connect(ui.actionComputeHKS, SIGNAL(triggered()), this, SLOT(computeHKS()));
	QObject::connect(ui.actionComputeHKSFeatures, SIGNAL(triggered()), this, SLOT(computeHksFeatures()));	
	QObject::connect(ui.actionComputeGeodesics, SIGNAL(triggered()), this, SLOT(computeGeodesics()));
	QObject::connect(ui.actionComputeHeatTransfer, SIGNAL(triggered()), this, SLOT(computeHeatTransfer()));
	QObject::connect(ui.actionComputeVertNormals, SIGNAL(triggered()), this, SLOT(computeVertNormals()));
	QObject::connect(ui.actionComputeFaceNormals, SIGNAL(triggered()), this, SLOT(computeFaceNormals()));
    QObject::connect(ui.actionComputeHoleNeighbors, SIGNAL(triggered()), this, SLOT(computeHoleNeighbors()));

	////  Edit  ////
	QObject::connect(ui.actionClearHandles, SIGNAL(triggered()), this, SLOT(clearHandles()));
	QObject::connect(ui.actionClone, SIGNAL(triggered()), this, SLOT(clone()));
	QObject::connect(ui.actionReconstructMHB, SIGNAL(triggered()), this, SLOT(reconstructMHB()));
	QObject::connect(ui.actionDeformSimple, SIGNAL(triggered()), this, SLOT(deformSimple()));
	QObject::connect(ui.actionDeformLaplace, SIGNAL(triggered()), this, SLOT(deformLaplace()));
	QObject::connect(ui.actionDeformLaplace2, SIGNAL(triggered()), this, SLOT(deformLaplace2()));
	QObject::connect(ui.actionDeformBiLaplace, SIGNAL(triggered()), this, SLOT(deformBiLaplace()));
	QObject::connect(ui.actionDeformMixedLaplace, SIGNAL(triggered()), this, SLOT(deformMixedLaplace()));
	QObject::connect(ui.actionDiffusionFlow, SIGNAL(triggered()), this, SLOT(diffusionFlow()));
    QObject::connect(ui.actionFillHoles, SIGNAL(triggered()), this, SLOT(fillHoles()));

    /* inpainting related */
    QObject::connect(ui.actionGenerateHoles, SIGNAL(triggered()), this, SLOT(generateHoles()));
    QObject::connect(ui.actionGenerateRingHoles, SIGNAL(triggered()), this, SLOT(generateRingHoles()));
    QObject::connect(ui.actionGenerateBandHole, SIGNAL(triggered()), this, SLOT(generateBandHole()));
    QObject::connect(ui.actionAutoGenHoles, SIGNAL(triggered()), this, SLOT(autoGenerateHoles()));
    QObject::connect(ui.actionIgnoreBoundary, SIGNAL(triggered()), this, SLOT(ignoreOuterBoundary()));
    QObject::connect(ui.actionDegradeHoles, SIGNAL(triggered()), this, SLOT(degradeHoles()));
    QObject::connect(ui.actionCutHoles, SIGNAL(triggered()), this, SLOT(cutHoles()));
    QObject::connect(ui.actionCutToSelected, SIGNAL(triggered()), this, SLOT(cutToSelected()));
    QObject::connect(ui.actionTriangulateHoles, SIGNAL(triggered()), this, SLOT(triangulateHoles()));
    QObject::connect(ui.actionRefineHoles, SIGNAL(triggered()), this, SLOT(refineHoles()));
    QObject::connect(ui.actionRefineHoles2, SIGNAL(triggered()), this, SLOT(refineHoles2()));
    QObject::connect(ui.actionRefineHoleByVertNum, SIGNAL(triggered()), this, SLOT(refineHolesByVertNum()));
    QObject::connect(ui.actionCopyMeshWithHoles, SIGNAL(triggered()), this, SLOT(copyMeshWithHoles()));
    QObject::connect(ui.actionEvaluateInpainting, SIGNAL(triggered()), this, SLOT(evaluateCurrentInpainting()));
    QObject::connect(ui.actionEvaluateInpaintingVert, SIGNAL(triggered()), this, SLOT(evaluateInpainting2()));
    QObject::connect(ui.actionHoleFairingLeastSquares, SIGNAL(triggered()), this, SLOT(fairHoleLeastSquares()));
    QObject::connect(ui.actionHoleFairingL1LS, SIGNAL(triggered()), this, SLOT(fairHoleL1LS()));
    QObject::connect(ui.actionHoleSmoothingDLRS, SIGNAL(triggered()), this, SLOT(smoothingHoleDLRS()));

    /* retrieval related */
    QObject::connect(ui.actionRegionByDistField, SIGNAL(triggered()), this, SLOT(regionByDistanceField()));

    //// Experiment ////
    QObject::connect(ui.actionExperiment1, SIGNAL(triggered()), this, SLOT(doExperiment1()));

	////  Display  ////
	QObject::connect(ui.actionDisplayMesh, SIGNAL(triggered()), this, SLOT(setDisplayMesh()));
    QObject::connect(ui.actionChangeShadeMode, SIGNAL(triggered()), ui.glMeshWidget, SLOT(changeShadeMode()));
    QObject::connect(ui.actionCurveSignature, SIGNAL(triggered()), this, SLOT(curveSignature()));
    QObject::connect(ui.actionSetColor, SIGNAL(triggered()), this, SLOT(setColor()));
    QObject::connect(ui.actionSetLegend, SIGNAL(triggered()), this, SLOT(setLegend()));

	QObject::connect(ui.actionDisplayWireframe, SIGNAL(triggered()), this, SLOT(setDisplayWireframe()));
	QObject::connect(ui.actionDisplayPointCloud, SIGNAL(triggered()), this, SLOT(setDisplayPointCloud()));
	QObject::connect(ui.actionDisplayNeighbors, SIGNAL(triggered()), this, SLOT(displayNeighborVertices()));
	QObject::connect(ui.actionShowFeatures, SIGNAL(triggered(bool)), this, SLOT(toggleShowFeatures(bool)));
	QObject::connect(ui.actionShowRefPoint, SIGNAL(triggered(bool)), this, SLOT(toggleShowRefPoint(bool)));
	QObject::connect(ui.actionShowSignature, SIGNAL(triggered(bool)), this, SLOT(toggleShowSignature(bool)));
	QObject::connect(ui.actionShowWireframeOverlay, SIGNAL(triggered(bool)), this, SLOT(toggleShowWireframeOverlay(bool)));
    QObject::connect(ui.actionShowHoles, SIGNAL(triggered(bool)), this, SLOT(toggleShowHoles(bool)));
    QObject::connect(ui.actionShowHoleHollow, SIGNAL(triggered(bool)), this, SLOT(toggleShowHoleHollow(bool)));
    QObject::connect(ui.actionShowHoleBoundary, SIGNAL(triggered(bool)), this, SLOT(toggleShowHoleBoundary(bool)));
    QObject::connect(ui.actionShowSurrounding, SIGNAL(triggered(bool)), this, SLOT(toggleShowSurrounding(bool)));
    QObject::connect(ui.actionShowHoleErrors, SIGNAL(triggered(bool)), this, SLOT(toggleShowHoleErrors(bool)));
    QObject::connect(ui.actionShowBbox, SIGNAL(triggered(bool)), this, SLOT(toggleShowBoundingBox(bool)));
	QObject::connect(ui.actionShowColorLegend, SIGNAL(triggered(bool)), this, SLOT(toggleShowColorLegend(bool)));
	QObject::connect(ui.actionShowVectors, SIGNAL(triggered(bool)), this, SLOT(toggleShowLines(bool)));
	
    QObject::connect(ui.actionDrawMatching, SIGNAL(triggered(bool)), this, SLOT(toggleDrawMatching(bool)));
	QObject::connect(ui.actionShowMatchingLines, SIGNAL(triggered(bool)), this, SLOT(toggleShowMatchingLines(bool)));
	QObject::connect(ui.actionDrawRegistration, SIGNAL(triggered(bool)), this, SLOT(toggleDrawRegistration(bool)));	
    
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
		gSettings.ACTIVE_COLOR_MAP_TYPE = ZGeom::ColorMapType((gSettings.ACTIVE_COLOR_MAP_TYPE + 1) % ZGeom::COLOR_MAP_COUNT);
        updateSignature();
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
		break;

	case Qt::Key_Y:
		break;

	case Qt::Key_Z:
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

void QZGeometryWindow::updateUI()
{
    updateMenuDisplaySignature();
    updateMenuDisplayFeatures();
    updateMenuDisplayLines();
    displaySignature(CMesh::StrAttrColorSigDefault.c_str());
}

bool QZGeometryWindow::initialize(const std::string& mesh_list_name)
{
	qout.output("******** Welcome ********", OUT_CONSOLE);
	qout.outputDateTime(OUT_CONSOLE);
	qout.output('*', 24, OUT_CONSOLE);

	switch (g_task) {
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
    if (vMeshFiles.size() < mMeshCount) {
        throw std::runtime_error("Not enough meshes in mesh list!");
    }

    bool scale_to_unit = (gSettings.INPUT_SCALE_TO_UNIT == 1 ? true : false);

	allocateStorage(mMeshCount);
	parallel_for(0, mMeshCount, [&](int obj) {
        loadMesh(vMeshFiles[obj], obj, scale_to_unit); 
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

void QZGeometryWindow::loadMesh(std::string mesh_filename, int obj, bool scale_to_unit /*= false*/)
{
    assert(obj < mMeshCount);
    CMesh *newMesh = new CMesh();
    newMesh->load(mesh_filename);
    newMesh->moveToOrigin();
    if (scale_to_unit) {
        double unit = 1.0;
        newMesh->scaleToUnitBox(unit);
    }    
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

void QZGeometryWindow::addMesh()
{
	QStringList filenames =  QFileDialog::getOpenFileNames(this, 
            "Select one or more mesh files to open", 
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

void QZGeometryWindow::removeCurrentMesh()
{
    mMeshHelper[0].removeCurrentMesh();
    mShapeEditor.init(mMeshHelper[0]);

    updateUI();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToNextMesh()
{
    mMeshHelper[0].nextMesh();
    mShapeEditor.init(mMeshHelper[0]);
    qout.outputStatus(QString("Switch to: ") + getMesh(0)->getMeshDescription().c_str());

    updateUI();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToPreviousMesh()
{
    mMeshHelper[0].prevMesh();
    mShapeEditor.init(mMeshHelper[0]);
    qout.outputStatus(QString("Switch to: ") + getMesh(0)->getMeshDescription().c_str());

    updateUI();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToNextCoordinate()
{
    QString coord_name = QString(mMeshHelper[0].getMesh()->switchNextCoordinate().c_str());
    qout.outputStatus("Coord: " + coord_name);
	
    ui.glMeshWidget->update();
}

void QZGeometryWindow::switchToPrevCoordinate()
{
    QString coord_name = QString(mMeshHelper[0].getMesh()->switchPrevCoordinate().c_str());
    qout.outputStatus("Coord: " + coord_name);
    
    ui.glMeshWidget->update();
}

void QZGeometryWindow::revertCoord()
{
    QString coord_name = QString(mMeshHelper[0].getMesh()->revertCoordinate().c_str());
    qout.outputStatus("Coord: " + coord_name);

	ui.glMeshWidget->update();
}

void QZGeometryWindow::deleteCurrentCoordinate()
{
    QString coord_name = QString(mMeshHelper[0].getMesh()->deleteCoordinate().c_str());
    qout.outputStatus("Coord: " + coord_name);

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
            || mLastOperation == Compute_HKS_Feature
		    || mLastOperation == Compute_MHWS || mLastOperation == Compute_MHW
		    || mLastOperation == Compute_SGWS || mLastOperation == Compute_SGW)
	{
		double time_scale;
		if (mCommonParameter <= sliderCenter) 
			time_scale = std::exp(std::log(gSettings.DEFAULT_HK_TIMESCALE / gSettings.MIN_HK_TIMESCALE) * ((double)mCommonParameter / (double)sliderCenter) + std::log(gSettings.MIN_HK_TIMESCALE));
		else 
			time_scale = std::exp(std::log(gSettings.MAX_HK_TIMESCALE / gSettings.DEFAULT_HK_TIMESCALE) * (double(mCommonParameter - sliderCenter) / sliderCenter) + std::log(gSettings.DEFAULT_HK_TIMESCALE)); 
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
	
    for (int obj = 0; obj < mMeshCount; ++obj) {
        if (!bToShow) {
            getMesh(obj)->colorize(CMesh::StrAttrColorSigDefault);
        }
        else {
            getMesh(obj)->colorize(mRenderManagers[obj].mActiveColorSignatureName);
        }
    }

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

void QZGeometryWindow::toggleShowHoleBoundary(bool show)
{
    bool bToShow = !ui.glMeshWidget->m_bShowHoleBoundary;
    ui.glMeshWidget->m_bShowHoleBoundary = bToShow;
    ui.actionShowHoles->setChecked(bToShow);

    ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowHoleHollow(bool show /*= false*/)
{
    bool bToShow = !ui.glMeshWidget->m_bShowHoleHollow;
    ui.glMeshWidget->m_bShowHoleHollow = bToShow;
    ui.actionShowHoles->setChecked(bToShow);

    ui.glMeshWidget->update();
}

void QZGeometryWindow::toggleShowHoleErrors(bool show /*= false*/)
{
    bool bToShow = !ui.glMeshWidget->m_bShowHoleError;
    ui.glMeshWidget->m_bShowHoleError = bToShow;
    ui.actionShowHoleErrors->setChecked(bToShow);

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
    using ZGeom::VertCurvature;
	for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh& mesh = *getMesh(obj);
        if (!mesh.hasAttr(ZGeom::StrAttrVertAllCurvatures)) {
            ZGeom::computeMeshCurvatures(mesh, true);
        }
        const vector<double>& vCM = ZGeom::getMeshCurvatures(*getMesh(obj), VertCurvature::MEAN);
        const vector<double>& vCG = ZGeom::getMeshCurvatures(*getMesh(obj), VertCurvature::GAUSS);
        const vector<double>& vCP1 = ZGeom::getMeshCurvatures(*getMesh(obj), VertCurvature::PRINCIPAL_1);
        const vector<double>& vCP2 = ZGeom::getMeshCurvatures(*getMesh(obj), VertCurvature::PRINCIPAL_2);

        ColorSignature colorCM(vCM, gSettings.ACTIVE_COLOR_MAP_TYPE);
        ColorSignature colorCG(vCG, gSettings.ACTIVE_COLOR_MAP_TYPE);
        ColorSignature colorCP1(vCP1, gSettings.ACTIVE_COLOR_MAP_TYPE);
        ColorSignature colorCP2(vCP2, gSettings.ACTIVE_COLOR_MAP_TYPE);

        getMesh(obj)->addColorSigAttr("color_mean_curvature", colorCM);
        getMesh(obj)->addColorSigAttr("color_gauss_curvature", colorCG);
        getMesh(obj)->addColorSigAttr("color_principal_curvature_1", colorCP1);
        getMesh(obj)->addColorSigAttr("color_principal_curvature_2", colorCP2);

        auto mm1 = std::minmax_element(vCM.begin(), vCM.end());
        auto mm2 = std::minmax_element(vCG.begin(), vCG.end());
        qout.output(QString("- Mean curvature -  min: %1, max: %2").arg(*mm1.first).arg(*mm1.second), OUT_TERMINAL);
        qout.output(QString("- Gauss curvature -  min: %1, max: %2").arg(QString::number(*mm2.first), QString::number(*mm2.second)), OUT_TERMINAL);
	}

	updateMenuDisplaySignature();
    displaySignature("color_mean_curvature");
    toggleShowColorLegend(true);
	qout.output("Mean curvature visualized");	
}

void QZGeometryWindow::computeShapeIndex()
{
    for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh& mesh = *getMesh(obj);
        ZGeom::computeShapeIndex(mesh);
        const vector<double>& vec_shape_index = mesh.getAttrValue<vector<double>>(ZGeom::StrAttrVertShapeIndex);
        ColorSignature colorSI(vec_shape_index, gSettings.ACTIVE_COLOR_MAP_TYPE);
        mesh.addColorSigAttr("color_shape_index", colorSI);
    }

    updateMenuDisplaySignature();
    displaySignature("color_shape_index");
    toggleShowColorLegend(true);
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
}

void QZGeometryWindow::reconstructMHB()
{
	int sliderCenter = ui.horizontalSliderParamter->maximum() / 2;
	double ratio = std::min((double)mCommonParameter/sliderCenter, 1.0);
	int nEig = mMeshHelper[0].getEigenSystem(CotFormula).eigVecCount() * ratio;
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
    if (hole == nullptr && mesh.hasAttr(ZGeom::StrAttrManualHoleRegions)) {
        auto &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrManualHoleRegions);
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

    vector<int> hole_verts(hole->vert_inside);
    for (int vi : hole->vert_on_boundary) hole_verts.push_back(vi);
    vector<int> surrounding_verts = ZGeom::vertSurroundingVerts(mesh, hole_verts, ring-1);
    for (int vi : hole->vert_on_boundary) surrounding_verts.push_back(vi);

    mesh.addAttrMeshFeatures(MeshFeatureList(surrounding_verts, ZGeom::ColorMagenta), "hole_ring_neighbor_verts");
    updateMenuDisplayFeatures();

    vector<int> surrounding_faces = ZGeom::getFaceEncompassedByVerts(mesh, surrounding_verts);
    mesh.addAttr<vector<int>>(surrounding_faces, StrAttrHoleSurroundingFaces, AR_UNIFORM, AT_VEC_INT);
    
    ui.glMeshWidget->update();
}

void QZGeometryWindow::computeEigenfunction()
{
	LaplacianType lap_type = active_lap_type;
    if (!mMeshHelper[0].hasEigenSystem(lap_type)) {
        QMessageBox::warning(this, "No Laplacian", "Compute Laplacian first!");
        return;
    }

    ZGeom::EigenSystem& es = mMeshHelper[0].getEigenSystem(lap_type);
    bool ok;
    int selected_eig = QInputDialog::getInt(this, "Select eigenfunction",
        "i-th eigenfunction", 1, 1, es.eigVecSize() - 1, 1, &ok);
    if (!ok) return;

	for (int obj = 0; obj < mMeshCount; ++obj) {
		MeshHelper& mp = mMeshHelper[obj];
		vector<double> vec_eigfunc = mp.getEigenSystem(lap_type).getEigVec(selected_eig).toStdVector();
        getMesh(obj)->addColorSigAttr(StrAttrColorEigenFunction, 
            ColorSignature(vec_eigfunc, gSettings.ACTIVE_COLOR_MAP_TYPE));
	}

	displaySignature(StrAttrColorEigenFunction.c_str());
	mLastOperation = Compute_Eig_Func;
    qout.output("Show eigenfunction #" + Int2String(selected_eig), OUT_CONSOLE);
	updateMenuDisplaySignature();
}

void QZGeometryWindow::computeBiharmonicDistField()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh* cur_mesh = getMesh(obj);
		MeshHelper& mp = mMeshHelper[obj];
		const int vertCount = getMesh(obj)->vertCount();
		const int refPoint = mp.getRefPointIndex();    
        ZGeom::EigenSystem& es = mp.getEigenSystem(CotFormula);
        vector<double> vVals = calAllBiharmonicDist(es, refPoint);

        addColorSignature(obj, vVals, StrAttrColorBiharmonicField);
	}

	displaySignature(StrAttrColorBiharmonicField.c_str());
	updateMenuDisplaySignature();
	mLastOperation = Compute_Biharmonic;
}

void QZGeometryWindow::computeHK()
{
	double time_scale = parameterFromSlider(gSettings.DEFAULT_HK_TIMESCALE, gSettings.MIN_HK_TIMESCALE, gSettings.MAX_HK_TIMESCALE);

	for (int obj = 0; obj < mMeshCount; ++obj) {
		MeshHelper& mp = mMeshHelper[obj];
		const int meshSize = mp.getMesh()->vertCount();
		const int refPoint = mp.getRefPointIndex();
		const ZGeom::EigenSystem &mhb = mp.getEigenSystem(active_lap_type);
        std::vector<double> values = ZGeom::calAllHeatKernel(mhb, refPoint, time_scale);
		addColorSignature(obj, values, StrAttrColorHK);
	}

	displaySignature(StrAttrColorHK.c_str());
	qout.output(QString().sprintf("HK with timescale: %f", time_scale));
	updateMenuDisplaySignature();
	mLastOperation = Compute_HK;
}

void QZGeometryWindow::computeHKS()
{
    mLastOperation = Compute_HKS_Feature;

	for (int i = 0; i < mMeshCount; ++i) {
		MeshHelper& mp = mMeshHelper[i];
		const int meshSize = mp.getMesh()->vertCount();
		const ZGeom::EigenSystem& mhb = mp.getEigenSystem(active_lap_type);

        double t_min = 4 * std::log(10.0) / mhb.getEigVal(299), t_max = 4 * std::log(10.0) / mhb.getEigVal(1);
        gSettings.MIN_HK_TIMESCALE = t_min; 
        gSettings.MAX_HK_TIMESCALE = t_max;
        gSettings.DEFAULT_HK_TIMESCALE = std::sqrt(t_min * t_max);
        double time_scale = parameterFromSlider(gSettings.DEFAULT_HK_TIMESCALE, t_min, t_max);
        
        std::vector<double> values = ZGeom::calHeatKernelSignature(mhb, time_scale);
		addColorSignature(i, values, StrAttrColorHKS);
        qout.output(QString().sprintf("HKS with time-scale: %f", time_scale));
	}

    displaySignature(StrAttrColorHKS.c_str());
	updateMenuDisplaySignature();
}

void QZGeometryWindow::computeHksFeatures()
{
	for (int i = 0; i < mMeshCount; ++i) {
        CMesh& mesh = *getMesh(i);
        MeshHelper& mp = mMeshHelper[i];
        const int meshSize = mp.getMesh()->vertCount();
        const ZGeom::EigenSystem& mhb = mp.getEigenSystem(active_lap_type);

        double t_min = 4 * std::log(10.0) / mhb.getAllEigVals().back(), t_max = 4 * std::log(10.0) / mhb.getEigVal(1);
        gSettings.MIN_HK_TIMESCALE = t_min;
        gSettings.MAX_HK_TIMESCALE = t_max;
        gSettings.DEFAULT_HK_TIMESCALE = std::sqrt(t_min * t_max);
        double time_scale = parameterFromSlider(gSettings.DEFAULT_HK_TIMESCALE, t_min, t_max);
        
        std::vector<double> values = ZGeom::calHeatKernelSignature(mhb, time_scale);
        vector<int> features = ZGeom::extractMeshExtrema(mesh, values, 2);
        mesh.addAttrMeshFeatures(features, "HKS features");
        updateMenuDisplayFeatures();
        
        if (mRenderManagers[i].mActivePointFeatures.find("HKS features") == mRenderManagers[i].mActivePointFeatures.end())
            displayFeature("HKS features");
	}
	
	if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures();
}

void QZGeometryWindow::computeSGW()
{
    for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh* cur_mesh = getMesh(obj);
        MeshHelper& mp = mMeshHelper[obj];
        const int vertCount = getMesh(obj)->vertCount();
        const int source_point = mp.getRefPointIndex();
        ZGeom::EigenSystem& es = mp.getEigenSystem(CotFormula);

        vector<double> timescales = computeSgwScales(es, 3);
        for (size_t s = 0; s < timescales.size(); ++s) {
            vector<double> vals = calAllSgwWavelet(es, timescales[s], source_point);
            cur_mesh->addColorSigAttr(StrAttrColorSGW[s], ZGeom::ColorSignature(vals, gSettings.ACTIVE_COLOR_MAP_TYPE, true));
        }
    }

    displaySignature(StrAttrColorSGW[0].c_str());
    updateMenuDisplaySignature();
}


void QZGeometryWindow::repeatOperation()
{
	switch(mLastOperation)
	{
	case Compute_Eig_Func:
		computeEigenfunction(); break;
	
	case Compute_HKS:
		computeHKS(); break;

	case Compute_HK:
		computeHK(); break;

	case Compute_Biharmonic:
		computeBiharmonicDistField(); break;

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
            getMesh(i)->colorize(sigName.toStdString());
        }        
	}
    
    for (QAction *qa : m_actionDisplaySignatures) {
        if (qa->text() == sigName) qa->setChecked(true);
        else qa->setChecked(false);        
    }

    if (ui.glMeshWidget->m_bShowSignature == false) {
        toggleShowSignature(true);
    }
	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayFeature(QString selected_feature_name)
{
	for (int obj = 0; obj < mMeshCount; ++obj) 
    {
		if (!isMeshSelected(obj)) continue;
        CMesh& mesh = *getMesh(obj);

        std::string new_feature_name = selected_feature_name.toStdString();
        auto &active_feature_names = mRenderManagers[obj].mActivePointFeatures;
        QAction *active_action = nullptr;
        for (QAction *qa : m_actionDisplayFeatures) {
            if (qa->text() == selected_feature_name) {
                active_action = qa;
                break;
            }
        }
        logic_assert(active_action != nullptr);

        if (!setHas(active_feature_names, new_feature_name)) {
            active_feature_names.insert(new_feature_name);
            active_action->setChecked(true);
            if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures(true);
        }
        else {
            active_feature_names.erase(new_feature_name);
            active_action->setChecked(false);            
        }
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::displayLine(QString line_feature_name)
{
    for (int obj = 0; obj < mMeshCount; ++obj) {
        if (!isMeshSelected(obj)) continue;
        auto &active_line_features = mRenderManagers[obj].mActiveLineNames;
        std::string new_line_feature = line_feature_name.toStdString();
        QAction *active_action = nullptr;
        for (QAction *qa : m_actionDisplayLines) {
            if (qa->text() == line_feature_name) {
                active_action = qa;
                break;
            }
        }
        logic_assert(active_action != nullptr);

        if (!setHas(active_line_features, new_line_feature)) {
            active_line_features.insert(new_line_feature);
            active_action->setChecked(true);
            if (!ui.glMeshWidget->m_bShowLines) toggleShowLines(true);
        }
        else {
            active_line_features.erase(new_line_feature);
            active_action->setChecked(false);
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

	if (!ui.glMeshWidget->m_bShowFeatures) toggleShowFeatures();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::matchFeatures()
{
	using namespace std;

	int force_matching = g_configMgr.getConfigValueInt("FORCE_MATCHING");
	string string_override = "";
	bool alreadyMached = false;
	
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

void QZGeometryWindow::computeGraphLaplacian()
{
    computeLaplacian(Umbrella);
}

void QZGeometryWindow::computeGeoLaplacian()
{
    computeLaplacian(CotFormula);
}

void QZGeometryWindow::computeLaplacian(LaplacianType lap_type)
{
	parallel_for(0, mMeshCount, [&](int obj) {
        if (!mMeshHelper[obj].hasLaplacian(lap_type))
		    mMeshHelper[obj].constructLaplacian(lap_type);
	});

	for (int obj = 0; obj < mMeshCount; ++obj) {
        int num_eig = getMesh(obj)->vertCount() - 2;
        num_eig = QInputDialog::getInt(this, "Select eig_num",
            QString().sprintf("Num of eigen (max: %d)", num_eig),
            num_eig, 1, num_eig, 1);

        MeshHelper& mp = mMeshHelper[obj];
        mp.computeLaplacian(num_eig, lap_type);
        auto& all_eigvals = mp.getEigenSystem(lap_type).getAllEigVals();
        std::cout << "Min EigVal: " << all_eigvals.front()
                  << "; Max EigVal: " << all_eigvals.back() 
                  << std::endl;
    }

    active_lap_type = lap_type;
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
	MeshHelper& mp = mMeshHelper[obj];
    CMesh& mesh = *getMesh(obj);
	
	if (!mp.getEigenSystem(laplacianType).empty()) return false; // already decomposed     
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
	const ZGeom::EigenSystem& mhb1 = mMeshHelper[0].getEigenSystem(CotFormula);
	const ZGeom::EigenSystem& mhb2 = mMeshHelper[1].getEigenSystem(CotFormula);
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
	
// 	funcMap1.print("output/funcmap1.txt");
// 	funcMap2.print("output/funcmap2.txt");
// 
// 	lap1.getLS().print("output/Ls1.txt");
// 	lap1.getW().print("output/W1.txt");
// 	lap2.getLS().print("output/Ls2.txt");
// 	lap2.getW().print("output/W2.txt");
// 
// 	mhb1.printEigVals("output/eigvals1.txt");
// 	mhb2.printEigVals("output/eigvals2.txt");
}

void QZGeometryWindow::addColorSignature( int obj, const std::vector<double>& vVals, const std::string& sigName )
{
	auto iResult = minmax_element(vVals.begin(), vVals.end());
	double sMin = *(iResult.first); 
	double sMax = *(iResult.second);
	std::cout << "-- Signature Min: " << sMin << ", Signature Max: " << sMax << std::endl;

    std::vector<double>& vSig = getMesh(obj)->addAttrVertScalars(StrAttrOriginalSignature).attrValue();
	vSig = vVals;
    getMesh(obj)->addColorSigAttr(sigName, ColorSignature(vVals, gSettings.ACTIVE_COLOR_MAP_TYPE, true));
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

void QZGeometryWindow::updateSignature()
{
	for (int obj = 0; obj < mMeshCount; ++obj) {
        CMesh& cur_mesh = *getMesh(obj);
        std::string cur_sig_name = mRenderManagers[0].mActiveColorSignatureName;
        if (cur_sig_name.empty() || !cur_mesh.hasAttr(cur_sig_name)) continue;
        ColorSignature& cur_sig = cur_mesh.getColorSignature(cur_sig_name);
        cur_sig.changeColorMap(gSettings.ACTIVE_COLOR_MAP_TYPE);
        cur_mesh.colorize(cur_sig_name);
	}

	ui.glMeshWidget->update();
}

void QZGeometryWindow::clearHandles()
{
	for (int obj = 0; obj < mMeshCount; ++obj) 
		mMeshHelper[obj].clearAllHandles();
	ui.glMeshWidget->update();
}

void QZGeometryWindow::deformLaplace2()
{
	mShapeEditor.deformLaplacian_v2();
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

void QZGeometryWindow::updateSignatureMin( int sMin )
{
    if (!getMesh(0)->hasAttr(StrAttrOriginalSignature)) return;
	if (sMin >= ui.sliderSigMax->value()) return;

    std::vector<double>& vSig = getMesh(0)->getAttrValue<std::vector<double>>(StrAttrOriginalSignature);
	auto mmp = std::minmax_element(vSig.begin(), vSig.end());
	double vMin = *mmp.first, vMax = *mmp.second;

	double newVal = (double)sMin / (double)ui.sliderSigMin->maximum() * (vMax - vMin) + vMin;
	ui.labelSigMin->setText("Min: " + QString::number(newVal));
	
    updateSignature();
}

void QZGeometryWindow::updateSignatureMax( int sMax )
{
    if (!getMesh(0)->hasAttr(StrAttrOriginalSignature)) return;
	if (sMax <= ui.sliderSigMin->value()) return;

    std::vector<double>& vSig = getMesh(0)->getAttrValue<std::vector<double>>(StrAttrOriginalSignature);
	auto mmp = std::minmax_element(vSig.begin(), vSig.end());
	double vMin = *mmp.first, vMax = *mmp.second;

	double newVal = (double)sMax / (double)ui.sliderSigMax->maximum() * (vMax - vMin) + vMin;
	ui.labelSigMax->setText("Max: " + QString::number(newVal));

    updateSignature();
}

void QZGeometryWindow::setLaplacianType( const QString& laplacianTypeName )
{
	qout.output("Select " + laplacianTypeName, OUT_STATUS);

	if (laplacianTypeName == "Umbrella") active_lap_type = Umbrella;
	else if (laplacianTypeName == "CotFormula") active_lap_type = CotFormula;
	else if (laplacianTypeName == "Anisotropic1") active_lap_type = Anisotropic1;
	else if (laplacianTypeName == "Anisotropic2") active_lap_type = Anisotropic2;
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

void QZGeometryWindow::setColor()
{
    vector<QString> property_names = { "hole", "wireframe", "boundary", "neighbor" };
    QStringList items;
    for (auto& str : property_names)
        items << str;

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select element"),
        tr("Element:"), items, 0, false, &ok);
    if (!ok) return;

    QColor *default_color(nullptr);
    if (item == "hole") default_color = &ui.glMeshWidget->m_regionColor;
    else if (item == "wireframe") default_color = &ui.glMeshWidget->m_wireframeColor;
    else if (item == "boundary") default_color = &ui.glMeshWidget->m_boundaryColor;
    else return;

    QColor new_color = QColorDialog::getColor(*default_color, this, "select new color");
    if (!new_color.isValid()) return;
    *default_color = new_color;

    ui.glMeshWidget->update();
}

void QZGeometryWindow::setLegend()
{
    vector<QString> property_names = { "width", "height" };
    QStringList items;
    for (auto& str : property_names)
        items << str;

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select element"),
        tr("Element:"), items, 0, false, &ok);
    if (!ok) return;

    int *dimension = nullptr;
    if (item == "width") dimension = &ui.glMeshWidget->m_colorBarWidth;
    else dimension = &ui.glMeshWidget->m_colorBarHeight;

    int newDim = QInputDialog::getInt(this, "Set color bar size", "size", *dimension, 1);
    *dimension = newDim;

    ui.glMeshWidget->update();
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

void QZGeometryWindow::generateHoles()
{
    int refIdx = mMeshHelper[0].getRefPointIndex();
    int holeVertCount = 30;

    bool ok;
    int i = QInputDialog::getInt(this, tr("Input hole size"),
        tr("Hole size:"), 25, 1, 10000, 1, &ok);
    if (ok) holeVertCount = i; 
    else return;
    
    MeshRegion generated_holes = ZGeom::generateRandomMeshRegion(*getMesh(0), vector<int>{refIdx}, holeVertCount);
    getMesh(0)->addAttr<vector<MeshRegion>>(vector < MeshRegion > {generated_holes}, ZGeom::StrAttrManualHoleRegions, AR_UNIFORM);    
    std::cout << "-- #inside_vert: " << generated_holes.vert_inside.size() << std::endl;

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
    int hole_count = 1;

    double missing_ratio = 0.2;
    missing_ratio = QInputDialog::getDouble(this, tr("Missing vertex ratio"),
        tr("missing_ratio:"), missing_ratio, 0.001, 0.75, 3, &ok);
    if (!ok) return;
    int holeVertCount = std::round(missing_ratio * N);

    vector<int> seedVerts;

    auto &handles = mMeshHelper[0].getHandles();
    if (!handles.empty()) {
        for (auto& elem : handles) seedVerts.push_back(elem.first);
        hole_count = handles.size();
    }
    else {
        hole_count = QInputDialog::getInt(this, tr("Input number of holes"),
            tr("#Holes"), hole_count, 0, holeVertCount, 1, &ok);
        if (!ok) return;
        if (hole_count <= 0) hole_count = holeVertCount;

        std::vector<int> all_verts(N);
        for (int i = 0; i < N; ++i) all_verts[i] = i;
        std::srand(std::time(0));
        random_unique<vector<int>::iterator>(all_verts.begin(), all_verts.end(), hole_count);
        seedVerts = vector < int > {all_verts.begin(), all_verts.begin() + hole_count};
    }

    vector<MeshRegion> generated_holes;
    for (int vi : seedVerts) {
        generated_holes.push_back(ZGeom::generateRandomMeshRegion(*getMesh(0), vector<int>{vi}, holeVertCount / hole_count));
    }

    if (hole_count < holeVertCount) {
        mergeMeshRegions(*getMesh(0), generated_holes);
    }

    getMesh(0)->addAttr<vector<MeshRegion>>(generated_holes, ZGeom::StrAttrManualHoleRegions, AR_UNIFORM);
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(getMeshRegionsInsideVerts(generated_holes), ZGeom::ColorGreen), "hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(getMeshRegionsBoundaryVerts(generated_holes), ZGeom::ColorRed), "hole_boundary_verts");
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

    MeshRegion generated_holes = ZGeom::generateRingMeshRegion(*getMesh(0), refIdx, ring);
    getMesh(0)->addAttr<vector<MeshRegion>>(vector < MeshRegion > {generated_holes}, ZGeom::StrAttrManualHoleRegions, AR_UNIFORM);
    std::cout << "#inside_vert: " << generated_holes.vert_inside.size() << std::endl;

    MeshRegion &hole = generated_holes;
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_inside, ZGeom::ColorGreen), "mesh_hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_on_boundary, ZGeom::ColorRed), "mesh_hole_boundary_verts");
    updateMenuDisplayFeatures();

    ui.glMeshWidget->update();
}

void QZGeometryWindow::generateBandHole()
{
    CMesh& mesh = *getMesh(0);
    int refIdx = mMeshHelper[0].getRefPointIndex();
    int inner_ring = 0;    

    bool ok;
    inner_ring = QInputDialog::getInt(this, tr("Input inner hole rings"),
        tr("inner ring"), inner_ring, 0, 100, 1, &ok);
    if (!ok) return;

    int outer_ring = inner_ring + 1;
    outer_ring = QInputDialog::getInt(this, tr("Input outer hole rings"),
        tr("outer ring"), outer_ring, inner_ring + 1, 100, 1, &ok);
    if (!ok) return;

    std::set<int> inner_verts{ refIdx }, outer_verts;

    if (inner_ring >= 1) {
        vector<int> ring_verts = ZGeom::vertSurroundingVerts(mesh, vector < int > {refIdx}, inner_ring);
        for (int vi : ring_verts) inner_verts.insert(vi);
    }

    vector<int> band_verts = ZGeom::vertSurroundingVerts(mesh, vector < int > {inner_verts.begin(), inner_verts.end()}, outer_ring - inner_ring);
    

    vector<MeshRegion> generated_holes{ ZGeom::meshRegionFromInsideVerts(mesh, band_verts) };

    getMesh(0)->addAttr<vector<MeshRegion>>(generated_holes, ZGeom::StrAttrManualHoleRegions, AR_UNIFORM);
    std::cout << "#inside_vert: " << generated_holes[0].vert_inside.size() << std::endl;

    MeshRegion &hole = generated_holes[0];
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_inside, ZGeom::ColorGreen), "mesh_hole_vertex");
    getMesh(0)->addAttrMeshFeatures(MeshFeatureList(hole.vert_on_boundary, ZGeom::ColorRed), "mesh_hole_boundary_verts");
    updateMenuDisplayFeatures();

    ui.glMeshWidget->update();
}

void QZGeometryWindow::cutHoles()
{
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    if (!original_mesh->hasAttr(ZGeom::StrAttrManualHoleRegions)) {
        std::cout << "No faces selected to cut" << std::endl;
        return;
    }
    vector<MeshRegion>& generated_holes = original_mesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrManualHoleRegions);
    vector<int> all_inside_faces = ZGeom::getMeshRegionsFaces(generated_holes);

    std::unique_ptr<CMesh> newMesh = std::move(ZGeom::cutFromMesh(*original_mesh, all_inside_faces));
        
    newMesh->initNamedCoordinates();
    mMeshHelper[0].addMesh(std::move(newMesh), "selected_partial_mesh");

    mShapeEditor.init(mMeshHelper[0]);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::cutToSelected()
{
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    if (!original_mesh->hasAttr(ZGeom::StrAttrManualHoleRegions)) {
        std::cout << "No faces selected to cut" << std::endl;
        return;
    }
    vector<MeshRegion>& generated_holes = original_mesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrManualHoleRegions);
    vector<int> all_inside_faces = ZGeom::getMeshRegionsFaces(generated_holes);

    std::unique_ptr<CMesh> newMesh = std::move(ZGeom::cutMeshTo(*mMeshHelper[0].getMesh(), all_inside_faces));
    
    newMesh->initNamedCoordinates();
    mMeshHelper[0].addMesh(std::move(newMesh), "hole_cut_mesh");

    mShapeEditor.init(mMeshHelper[0]);
    ui.glMeshWidget->update();
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
    const auto& vHoles = ZGeom::getMeshBoundaryLoops(*oldMesh);
    if (!vHoles.empty()) 
    {
        if (vHoles[0].face_inside.empty()) {
            std::cout << "Hole need triangulation first!" << std::endl;
            return;
        }

        std::unique_ptr<CMesh> newMesh = std::make_unique<CMesh>(*oldMesh);
        double lambda = 0.7;
        bool ok;
        lambda = QInputDialog::getDouble(this, tr("Input refine coefficient"),
            tr("lambda:"), lambda, 0.1, 2, 2, &ok);
        if (!ok) return;

        ZGeom::refineMeshHoles3(*newMesh.get(), lambda);
        mMeshHelper[0].addMesh(std::move(newMesh), "hole_refined_mesh");
        std::cout << "Mesh hole refined!" << std::endl;

        evaluateCurrentInpainting();
        ui.glMeshWidget->update();
    }
    else std::cout << "No holes found! Do nothing." << std::endl;
}

void QZGeometryWindow::refineHoles2()
{
    CMesh* oldMesh = mMeshHelper[0].getMesh();
    const auto& vHoles = ZGeom::getMeshBoundaryLoops(*oldMesh);
    if (!vHoles.empty())
    {
        if (vHoles[0].face_inside.empty()) {
            std::cout << "Hole need triangulation first!" << std::endl;
            return;
        }

        std::unique_ptr<CMesh> newMesh = std::make_unique<CMesh>(*oldMesh);
        double lambda = 0.7;
        bool ok;
        lambda = QInputDialog::getDouble(this, tr("Input refine coefficient"),
            tr("lambda:"), lambda, 0.1, 2, 2, &ok);
        if (!ok) return;

        ZGeom::refineMeshHoles2(*newMesh.get(), lambda);
        mMeshHelper[0].addMesh(std::move(newMesh), "hole_refined_mesh");
        std::cout << "Mesh hole refined!" << std::endl;

        evaluateCurrentInpainting();
        ui.glMeshWidget->update();
    }
    else std::cout << "No holes found! Do nothing." << std::endl;
}

void QZGeometryWindow::refineHolesByVertNum()
{
    CMesh* oldMesh = mMeshHelper[0].getMesh();
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();

    vector<MeshRegion> vHoles = ZGeom::getMeshBoundaryLoops(*oldMesh);
    if (vHoles.empty()) {
        std::cout << "No holes found! Do nothing." << std::endl;
        return;
    }

    if (vHoles[0].face_inside.empty()) {
        std::cout << "Hole need triangulation first!" << std::endl;
        return;
    }

    if (!original_mesh->hasAttr(ZGeom::StrAttrManualHoleRegions)) return;
    int original_inside_vert_count = original_mesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrManualHoleRegions)[0].vert_inside.size();

    int new_vert_count = original_inside_vert_count;
    bool ok(false);
    new_vert_count = QInputDialog::getInt(this, tr("Input number of new vertices to add"),
        tr("#New Vertices"), new_vert_count, 0, 10000, 1, &ok);
    if (!ok) return;

    std::unique_ptr<CMesh> newMesh = std::make_unique<CMesh>(*oldMesh);
    
    if (!ZGeom::refineMeshHoleByNum(*newMesh.get(), vHoles[0], new_vert_count)) {
        std::cout << "Fail to refine mesh hole!" << std::endl;
        return;
    }
    
    mMeshHelper[0].addMesh(std::move(newMesh), "hole_refined_mesh");
    std::cout << "Mesh hole refined!" << std::endl;

    evaluateCurrentInpainting();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::copyMeshWithHoles()
{
    CMesh* oldMesh = mMeshHelper[0].getMesh();
    if (oldMesh->hasAttr(ZGeom::StrAttrManualHoleRegions)) 
    {
        auto& manualHoles = oldMesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrManualHoleRegions);
        if (!manualHoles.empty()) {
            std::unique_ptr<CMesh> newMesh(new CMesh(*oldMesh));
            newMesh->addAttr<vector<MeshRegion>>(manualHoles, ZGeom::StrAttrMeshHoleRegions, AR_UNIFORM);
            newMesh->removeAttr(ZGeom::StrAttrManualHoleRegions);
            mMeshHelper[0].addMesh(std::move(newMesh), "hole_copied_mesh");
            std::cout << "Manual hole copied to new mesh!" << std::endl;

            ui.glMeshWidget->update();
        }
    }
    else std::cout << "No manual hole found in the current mesh" << std::endl;
}

void QZGeometryWindow::evaluateCurrentInpainting()
{
    using ZGeom::Vec3d;

    /* compare inpainting result with the original generated mesh */
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    CMesh* cur_mesh = mMeshHelper[0].getMesh();
    if (!original_mesh->hasAttr(ZGeom::StrAttrManualHoleRegions) || !cur_mesh->hasAttr(ZGeom::StrAttrMeshHoleRegions)) return;
    
    vector<MeshRegion>& manualHoles = original_mesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrManualHoleRegions);
    vector<MeshRegion>& refinedHoles = cur_mesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (refinedHoles.size() == manualHoles.size() && refinedHoles.size() > 0) 
    {
        vector<double> vertArea1 = ZGeom::getMeshVertMixedAreas(*original_mesh);
        vector<double> vertArea2 = ZGeom::computeMeshVertArea(*cur_mesh);

        vector<int> vFaces1 = getMeshRegionsFaces(manualHoles);
        vector<int> vVerts1 = getMeshRegionsInsideVerts(manualHoles); 
        vector<int> vFaces2 = getMeshRegionsFaces(refinedHoles);
        vector<int> vVerts2 = getMeshRegionsInsideVerts(refinedHoles);

        vector<vector<Vec3d>> meshTri1(vFaces1.size()), meshTri2(vFaces2.size());
        for (int i = 0; i < (int)vFaces1.size(); ++i)
            meshTri1[i] = original_mesh->getFace(vFaces1[i])->getAllVertPos();
        for (int i = 0; i < (int)vFaces2.size(); ++i)
            meshTri2[i] = cur_mesh->getFace(vFaces2[i])->getAllVertPos();

        vector<double> vertDist1(original_mesh->vertCount(), 0);
        vector<double> vertDist2(cur_mesh->vertCount(), 0);

        concurrency::parallel_for_each(vVerts1.begin(), vVerts1.end(), [&](int vi) {
            const ZGeom::Vec3d &vPos = original_mesh->vertPos(vi);
            double minDistVi = 1e15;
            for (const vector<Vec3d>& tri : meshTri2)
                minDistVi = std::min(minDistVi, distPointTriangle(vPos, tri).distance);
            vertDist1[vi] = minDistVi;
        });

        concurrency::parallel_for_each(vVerts2.begin(), vVerts2.end(), [&](int vi) {
            const Vec3d &vPos = cur_mesh->vertPos(vi);
            double minDistVi = 1e15;
            for (const vector<Vec3d>& tri : meshTri1)
                minDistVi = std::min(minDistVi, distPointTriangle(vPos, tri).distance);
            vertDist2[vi] = minDistVi;
        });

        double rmse1 = rmseVerts(vVerts1, vertArea1, vertDist1);
        double rmse2 = rmseVerts(vVerts2, vertArea2, vertDist2);
        std::cout << "-- Inpainting Error new to original: " << rmse2 << "\n";
        std::cout << "-- Inpainting Error original to new: " << rmse1 << std::endl;

        cur_mesh->addColorSigAttr(StrAttrColorInpaintError, ZGeom::ColorSignature(vertDist2, ZGeom::CM_PARULA, true));
        cur_mesh->getColorSignature(StrAttrColorInpaintError).curve(0, inpainting_error_curving_max);

        updateMenuDisplaySignature();
        ui.glMeshWidget->update();
    }    
}

void QZGeometryWindow::evaluateInpainting2()
{
    using namespace ZGeom;
    CMesh* original_mesh = mMeshHelper[0].getOriginalMesh();
    CMesh* cur_mesh = mMeshHelper[0].getMesh();
    
    if (!cur_mesh->hasAttr(StrAttrMeshHoleRegions)) return;
    if (original_mesh->vertCount() != cur_mesh->vertCount()) return;

    vector<MeshRegion>& refinedHoles = cur_mesh->getAttrValue<vector<MeshRegion>>(StrAttrMeshHoleRegions);
    vector<int> inside_verts = getMeshRegionsInsideVerts(refinedHoles);
    
    vector<double> vErrors(cur_mesh->vertCount(), 0);
    double errorSum(0);
    for (int vi : inside_verts) {
        vErrors[vi] = (cur_mesh->vertPos(vi) - original_mesh->vertPos(vi)).length();
        errorSum += ZGeom::sqr(vErrors[vi]);        
    }

    cur_mesh->addColorSigAttr(StrAttrColorInpaintError, ZGeom::ColorSignature(vErrors, ZGeom::CM_PARULA, true));
    cur_mesh->getColorSignature(StrAttrColorInpaintError).curve(0, inpainting_error_curving_max);

    updateMenuDisplaySignature();
    ui.glMeshWidget->update();

    double rmse = std::sqrt(errorSum / (double)inside_verts.size());
    std::cout << "-- Vert RMSE: " << rmse << "\n";
}

void QZGeometryWindow::curveSignature()
{
    CMesh* cur_mesh = getMesh(0);
    std::string cur_sig_name = mRenderManagers[0].mActiveColorSignatureName;
    if (cur_sig_name.empty() || !cur_mesh->hasAttr(cur_sig_name)) return;
    ColorSignature& cur_sig = cur_mesh->getColorSignature(cur_sig_name);
    if (!cur_sig.hasOriginalValues()) return;

    double lower_bound = cur_sig.getCurveMin();
    double upper_bound = cur_sig.getCurveMax();

    bool ok;
    lower_bound = QInputDialog::getDouble(this, tr("Curve min"), 
            tr("Curve min"), lower_bound, 
            (double)INT_MIN, (double)INT_MAX, 4, &ok);
    if (!ok) return;
    upper_bound = QInputDialog::getDouble(this, tr("Curve max"),
            tr("Curve max"), upper_bound, 
            lower_bound, (double)INT_MAX, 4, &ok);
    if (!ok) return;

    cur_sig.curve(lower_bound, upper_bound);
    cur_mesh->colorize(cur_sig_name);
    ui.glMeshWidget->update();
}

void QZGeometryWindow::fairHoleLeastSquares()
{
    CMesh& mesh = *mMeshHelper[0].getMesh();
    vector<MeshRegion> &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) {
        std::cout << "No holes found!" << std::endl;
        return;
    }

    /* input parameters */
    int anchor_ring = 3;
    double anchor_weight = 1.0;
    bool ok;

    anchor_ring = QInputDialog::getInt(this, tr("Input surrounding ring"),
        tr("Ring:"), anchor_ring, 0, 50, 1, &ok);
    if (!ok) return;  
    anchor_weight = QInputDialog::getDouble(this, tr("Input anchor weight"),
        tr("Weight:"), anchor_weight, 0.1, 10000, 2, &ok);
    if (!ok) return;

    MeshCoordinates coord_ls = least_square_hole_inpainting(mesh, vHoles, anchor_ring, anchor_weight);
    mesh.addNamedCoordinate(coord_ls, "ls_hole_fairing");
    
    evaluateCurrentInpainting();
    //evaluateInpainting2();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::fairHoleThinPlateEnergy()
{
    CMesh& mesh = *mMeshHelper[0].getMesh();
    vector<MeshRegion> &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) {
        std::cout << "No holes found!" << std::endl;
        return;
    }
    const ZGeom::MeshRegion& hole_region = vHoles[0];

    int nIter = 1;
    double eps = 1e-3;

    MeshCoordinates coord_ls = thin_plate_energy_hole_inpainting(mesh, hole_region, nIter, eps);
    mesh.addNamedCoordinate(coord_ls, "thin_plate_hole_fairing");

    evaluateCurrentInpainting();
    ui.glMeshWidget->update();
}

void QZGeometryWindow::fairHoleL1LS()
{
    CMesh& mesh = *mMeshHelper[0].getMesh();
    vector<MeshRegion> &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) {
        std::cout << "No holes found!" << std::endl;
        return;
    }

    /* input parameters */
    ParaL1LsInpainting para;
    para.fitting_ring = 0;
    para.eigen_count = -1;
    para.tol = 1e-3;
    para.lambda = 1e-3;

    if (!getfairHoleL1LsParameters(this, para)) return;

    MeshCoordinates coord_ls = l1_ls_hole_inpainting(mesh, vHoles, para);
    mesh.addNamedCoordinate(coord_ls, "l1ls_hole_inpainting");
    
    ui.glMeshWidget->update();
}

void QZGeometryWindow::doExperiment1()
{
    std::cout << "Experiment 1 - Random vertex recovery with different eigendecomposition" << std::endl;

    CMesh& mesh = *mMeshHelper[0].getMesh();
    const int total_vert_count = mesh.vertCount();
    vector<MeshRegion> &hole_regions = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    vector<int> hole_verts = getMeshRegionsInsideVerts(hole_regions);
    if (hole_verts.empty()) {
        std::cerr << "Missing vertices not specified!\nDid you forgot to copy the mesh with holes?" << std::endl;
        return;
    }

    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(&mesh);
    int eigenCount = total_vert_count - 1;

    ZGeom::EigenSystem es;
    graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, es);
    
    ZGeom::Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    MeshCoordinates coord_old = mesh.getVertCoordinates();
    ZGeom::DenseMatrixd mat_coord_old = mesh.getVertCoordinates().toDenseMatrix();
    
    double tol = 1e-3;
    double lambda = 1e-3;

    CStopWatch timer;
    for (int i = 1; i <= 10; i+=1) {
        int dict_size = int(eigenCount * (double(i) / 10.0));
        ZGeom::DenseMatrixd mat_dict = dictMHB.toDenseMatrix(dict_size);
        std::cout << "- " << i <<  ": #eigen = " << dict_size << std::endl;

        timer.startTimer();
        ZGeom::DenseMatrixd mat_coord_inpainted = matlab_inpaintL1LS(mat_coord_old, mat_dict, hole_verts, lambda, tol);
        timer.stopTimer("-- L1_Ls inpainting time: ");

        MeshCoordinates coord_inpainted;
        coord_inpainted.fromDenseMatrix(mat_coord_inpainted);
        double errorSum(0);
        for (int vi : hole_verts) {
            double error = (coord_inpainted[vi] - coord_old[vi]).length();
            errorSum += error * error;
        }

        double rmse = std::sqrt(errorSum / (double)hole_verts.size());
        std::cout << "-- RMSE: " << rmse << std::endl;

        // use the following 4 lines if only partial assignment
        MeshCoordinates result(coord_old);
        for (int vi : hole_verts) {
            result.setVertCoord(vi, coord_inpainted[vi]);
        }

        // use the following line if full assignment
        //MeshCoordinates& result = coord_inpainted;
        mesh.addNamedCoordinate(result, "l1ls_hole_inpainting_" + Int2String(i));
    }
}

void QZGeometryWindow::degradeHoles()
{
    CMesh& mesh = *getMesh(0);
    if (!mesh.hasAttr(ZGeom::StrAttrMeshHoleRegions)) return;
    vector<MeshRegion> &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) {
        std::cout << "No holes found!" << std::endl;
        return;
    }
    vector<int> selectedVerts = getMeshRegionsInsideVerts(vHoles);

    double phi = 0.02;
    bool ok;
    phi = QInputDialog::getDouble(this, tr("Add Gauss noise to mesh"),
        tr("phi"), phi, 0.001, 1, 3, &ok);
    if (!ok) return;


    MeshCoordinates noisyCoord = ZGeom::addMeshNoise(mesh, phi, selectedVerts);
    mesh.addNamedCoordinate(noisyCoord, "hole_noisy_coord");
    std::cout << "Add Gauss noise with phi=" << phi << std::endl;
    ui.glMeshWidget->update();
}

void QZGeometryWindow::smoothingHoleDLRS()
{
    CMesh& mesh = *getMesh(0);
    if (!mesh.hasAttr(ZGeom::StrAttrMeshHoleRegions)) return;
    vector<MeshRegion> &vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) {
        std::cout << "No holes found!" << std::endl;
        return;
    }
    vector<int> selectedVerts = getMeshRegionsInsideVerts(vHoles);

    double lambda = 0.8;
    int anchor_ring = 5;
    bool ok;
    anchor_ring = QInputDialog::getInt(this, tr("Input surrounding ring"),
        tr("Ring:"), anchor_ring, 0, 50, 1, &ok);
    if (!ok) return;
    lambda = QInputDialog::getDouble(this, tr("DLRS denoising coefficient"),
        tr("lambda"), lambda, 0, 100, 3, &ok);
    if (!ok) return;
    
    MeshCoordinates denoisedCoord = meshHoleDLRS(mesh, vHoles, lambda, anchor_ring);
    mesh.addNamedCoordinate(denoisedCoord, "hole_denoised");
    std::cout << "Denoising with DLRS; lambda = " << lambda << std::endl;
    ui.glMeshWidget->update();
}

void QZGeometryWindow::ignoreOuterBoundary()
{
    CMesh& mesh = *getMesh(0);
    vector<MeshRegion>& vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    auto iter = std::max_element(vHoles.begin(), vHoles.end(), [](const MeshRegion& mr1, const MeshRegion& mr2) {
        return mr1.he_on_boundary.size() < mr2.he_on_boundary.size();
    });

    std::cout << "Max boundary length: " << iter->he_on_boundary.size() << std::endl;
    vHoles.erase(iter);

    ui.glMeshWidget->update();
}

void QZGeometryWindow::regionByDistanceField()
{
    CMesh& mesh = *getMesh(0);
    computeBiharmonicDistField();
    const vector<double>& biharm_dist = getMesh(0)->getAttrValue<ColorSignature>(StrAttrColorBiharmonicField).getOriginalVals();

    double dist_max = *std::max_element(biharm_dist.begin(), biharm_dist.end());
    double threshold = 0.25 * dist_max;
    bool ok;
    threshold = QInputDialog::getDouble(this, tr("Contour distance threshold"),
        tr("Threshold:"), threshold, 0, dist_max, 3, &ok);
    if (!ok) return;

    int source_point = mMeshHelper[0].getRefPointIndex();
    MeshRegion distRegion = ZGeom::meshRegionFromDistField(mesh, biharm_dist, 
        source_point, [=](double val) { return val <= threshold; });

    mesh.addAttr<vector<MeshRegion>>(vector<MeshRegion>{distRegion}, ZGeom::StrAttrManualHoleRegions, AR_UNIFORM);
    ui.glMeshWidget->update();
}
