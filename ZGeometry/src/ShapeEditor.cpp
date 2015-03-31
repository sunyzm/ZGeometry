#include "ShapeEditor.h"
#include <iostream>
#include <random>
#include <functional>
#include <iomanip>
#include <array>
#include <queue>
#include <ppl.h>
#include <boost/lexical_cast.hpp>
#include <ZGeom/ZGeom.h>
#include <ZGeom/util.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/sparse_approximation.h>
#include <ZGeom/SparseSolver.h>
#include <ZGeom/MCA.h>
#include <ZGeom/Geodesic.h>
#include "differential_geometry.h"
#include "heat_diffusion.h"
#include "global.h"

using std::vector;
using std::array;
using std::min;
using std::max;
using ZGeom::Vec3d;
using ZGeom::VecNd;
using ZGeom::SparseMatrix;
using ZGeom::Dictionary;
using ZGeom::SparseCoding;
using ZGeom::SparseApproxMethod;
using ZGeom::combineDictionary;
using ZGeom::DenseMatrix;
using ZGeom::DenseMatrixd;
using ZGeom::Colorf;
using ZGeom::MeshRegion;
using ZGeom::WeightSet;

std::string ShapeEditor::strOriginalCoord = "original";

bool readPursuits(const std::string& pursuitFile, ZGeom::SparseCoding *vPursuits[]) 
{
	if (!fileExist(pursuitFile)) 
		return false;

	std::ifstream ifs(pursuitFile.c_str());
	int atomCount;
	ifs >> atomCount;
	
	for (int b = 0; b < 3; ++b) {
		vPursuits[b]->resize(atomCount);
		for (int i = 0; i < atomCount; ++i) {
			ZGeom::SparseCodingItem& item = (*vPursuits[b])[i];
			ifs >> item.index() >>item.coeff();
		}
	}

	return true;
}

void writePursuits(const std::string& pursuitFile, ZGeom::SparseCoding *vPursuits[]) {
	std::ofstream ofs(pursuitFile.c_str());
	int atomCount = vPursuits[0]->size();
	ofs << atomCount << std::endl;
	for (int b = 0; b < 3; ++b) {
		for (int i = 0; i < atomCount; ++i) {
			const ZGeom::SparseCodingItem& item = (*vPursuits[b])[i];
			ofs << item.index() << ' ' << item.coeff() << std::endl;
		}
	}	
}

void colorPartitions(const std::vector<int>& partIdx, const Palette& partPalette, std::vector<ZGeom::Colorf>& vColors)
{
	int nPart = partPalette.totalColors();
	int nPoints = partIdx.size();
	vColors.resize(nPoints);
	
	for (int i = 0; i < nPoints; ++i)
		vColors[i] = partPalette.getColor(partIdx[i]);
}

/* Discrete Laplacian Regularization Smoothing */
vector<VecNd> DLRS(ZGeom::MatlabEngineWrapper& engine, const SparseMatrix<double>& matL, double lambda, const vector<VecNd>& vSignals)
{
    assert(matL.rowCount() == matL.colCount());
    int nOrder = matL.rowCount();
    int nChannel = vSignals.size();
    engine.addSparseMat(matL, "matL");
    engine.addDoubleScalar(lambda, "lambda");
    ZGeom::DenseMatrixd matSignal(vSignals);
    engine.addDenseMat(matSignal, "matSignals", true);
    
    // solve (I + lambda*L^T*L) S = P
    // cf. Decoupling noise and features via weighted L1-analysis compressed sensing (SIGGRAPH 2014)
    engine.eval("matDenoised=(eye(size(matL))+lambda*matL'*matL) \\ matSignals;");
    
    DenseMatrixd matDenoised = engine.getDenseMat("matDenoised", true);
    return matDenoised.toRowVecs();
}

vector<ZGeom::Colorf> colorCoordDiff(const MeshCoordinates& coordDiff, double diffMax, ZGeom::ColorMapType colorMapType = ZGeom::CM_JET)
{
    int vertCount = coordDiff.size();
    vector<ZGeom::Colorf> vColors(vertCount);
    for (int i = 0; i < vertCount; ++i) {
        double diff = fabs(coordDiff.getVertCoordinate(i).length());
        vColors[i].falseColor(diff / diffMax, colorMapType);
    }
    return vColors;
}

void printBeginSeparator(std::string s, char c) {
	std::cout << '\n';
	for (int i = 0; i < 8; ++i) std::cout << c;
	std::cout << ' ' << s << ' ';
	for (int i = 0; i < 8; ++i) std::cout << c;
	std::cout << '\n';
}

void printEndSeparator(char c, int num) {
	std::cout << std::setfill(c) << std::setw(num) << '\n';
}


void ShapeEditor::init(MeshHelper& processor)
{
    mMeshHelper = &processor;
    mMesh = mMeshHelper->getMesh();
    mTotalScales = 0;
    mCurCoordID = 0;
    resetStoredCoordinates();
}

void ShapeEditor::resetStoredCoordinates()
{
    mStoredCoordinates.resize(1);
    mStoredCoordinates[0] = std::make_pair(strOriginalCoord, mMesh->getVertCoordinates());
    mCurCoordID = 0;
}

void ShapeEditor::revertCoordinates()
{
    changeCoordinates(0);
}


void ShapeEditor::changeCoordinates( int coordID )
{
	if (coordID < 0 || coordID >= mStoredCoordinates.size()) return;
	mCurCoordID = coordID;
	qout.output("Change coord: " + mStoredCoordinates[mCurCoordID].first, OUT_CONSOLE);

	if (mStoredCoordinates[mCurCoordID].second.empty())
		qout.output("!! Selected coordinate is empty!", OUT_CONSOLE);
	else 
		mMesh->setVertCoordinates(mStoredCoordinates[mCurCoordID].second);
}

void ShapeEditor::changeCoordinates(std::string coordName)
{
    for (int i = 0; i < (int)mStoredCoordinates.size(); ++i) {
        if (coordName == mStoredCoordinates[i].first) {
            changeCoordinates(i); break;
        }
    }
}

void ShapeEditor::setStoredCoordinates(const MeshCoordinates& newCoord, int idx, std::string coordName)
{
	if (idx >= mStoredCoordinates.size()) mStoredCoordinates.resize(idx + 1);
	mStoredCoordinates[idx] = std::make_pair(coordName, newCoord);
}

MeshCoordinates ShapeEditor::getStoredCoordinate( int idx )
{
    if (idx >= mStoredCoordinates.size()) return MeshCoordinates();
	else return mStoredCoordinates[idx].second;
}

void ShapeEditor::addCoordinate(const MeshCoordinates &newCoord, std::string coordName)
{
    for (int i = 0; i < (int)mStoredCoordinates.size(); ++i) {
        if (mStoredCoordinates[i].first == coordName) {
            mStoredCoordinates[i].second = newCoord;
            return;
        }
    }
    mStoredCoordinates.push_back(std::make_pair(coordName, newCoord));
}

void ShapeEditor::prepareAnchors( int& anchorCount, std::vector<int>& anchorIndex, std::vector<ZGeom::Vec3d>& anchorPos ) const
{
	const auto& anchors = mMeshHelper->getHandles();
	anchorCount = anchors.size();
	anchorIndex.clear();
	anchorPos.clear();
	for (auto a : anchors) {
		anchorIndex.push_back(a.first);
		anchorPos.push_back(a.second);
	}
}

MeshCoordinates ShapeEditor::getNoisyCoord(double phi)
{
	assert(phi > 0 && phi < 1);
	const int vertCount = mMesh->vertCount();
	const double avgLen = mMesh->getAvgEdgeLength();
    double bbDiag = mMesh->getBoundingBox().length() * 2;
    MeshCoordinates newCoord = mMesh->getVertCoordinates();

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0, phi);

	for (int vIdx = 0; vIdx < vertCount; ++vIdx) {
		for (int a = 0; a < 3; ++a) {
			double noise = bbDiag * distribution(generator);
			newCoord.getCoordFunc(a)[vIdx] += noise;
		}
	}

    return newCoord;
}

void ShapeEditor::addColorSignature(const std::string& colorSigName, const std::vector<ZGeom::Colorf>& vColors)
{
    mMesh->addColorSigAttr(colorSigName, ZGeom::ColorSignature(vColors));
    emit meshSignatureAdded();
}

void ShapeEditor::continuousReconstruct( int selected, int contCoordIdx )
{
	if (selected < 0 || selected >= 4 || contCoordIdx < 0 || 
		contCoordIdx >= mContReconstructCoords[selected].size()) return;
	if (mContReconstructCoords[selected].empty()) return;

	mMesh->setVertCoordinates(mContReconstructCoords[selected][contCoordIdx]);

	/* visualize the coordinate difference against the original coordinate */
	ZGeom::ColorSignature cs = mMesh->addColorSigAttr(StrAttrColorPosDiff).attrValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, cs.getColors());
    emit meshSignatureAdded();
	emit coordinateSelected(selected, contCoordIdx);

}

void ShapeEditor::fourierReconstruct( int nEig )
{
	const int vertCount = mMesh->vertCount();
    MeshCoordinates oldCoord = mMesh->getVertCoordinates();
	MeshCoordinates newCoord(vertCount);
	LaplacianType lapType = SymCot;
		
	ZGeom::SparseMatrixCSR<double, int> matW;
	mMeshHelper->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ZGeom::EigenSystem &mhb = mMeshHelper->getMHB(lapType);
	ZGeom::logic_assert(nEig <= mhb.eigVecCount(), "Insufficient eigenvectors in mhb");

	const ZGeom::VecNd &vx = oldCoord.getCoordFunc(0),
					   &vy = oldCoord.getCoordFunc(1),
					   &vz = oldCoord.getCoordFunc(2);
	ZGeom::VecNd &xCoord = newCoord.getCoordFunc(0),
				 &yCoord = newCoord.getCoordFunc(1),
				 &zCoord = newCoord.getCoordFunc(2);
	std::vector<double> xCoeff(nEig), yCoeff(nEig), zCoeff(nEig);

	for (int i = 0; i < nEig; ++i) {
		const ZGeom::VecNd& eigVec = mhb.getEigVec(i);
		xCoeff[i] = ZGeom::innerProductSym(vx, matW, eigVec);
		yCoeff[i] = ZGeom::innerProductSym(vy, matW, eigVec);
		zCoeff[i] = ZGeom::innerProductSym(vz, matW, eigVec);
		
		xCoord += xCoeff[i] * eigVec;
		yCoord += yCoeff[i] * eigVec;
		zCoord += zCoeff[i] * eigVec;
	}

	mMesh->setVertCoordinates(newCoord);
}

void ShapeEditor::deformSimple()
{
	CStopWatch timer;	
	timer.startTimer();

	int anchorCount(0);
	std::vector<int> anchorIndex;
	std::vector<ZGeom::Vec3d> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int vertCount = mMesh->vertCount();
	int hIdx = anchorIndex[0];		// only use the first handle to deform
	ZGeom::Vec3d handleTrans = anchorPos[0] - mMesh->vert(hIdx)->pos();

	std::vector<int> vFreeIdx = mMesh->getVertNeighborVerts(hIdx, 5, true);
//	vFreeIdx.resize(vertCount);
//	for (int i = 0; i < vertCount; ++i) vFreeIdx[i] = i;

	const int freeVertCount = vFreeIdx.size();
	std::vector<double> vDist2Handle(freeVertCount);
	concurrency::parallel_for (0, freeVertCount, [&](int i) {
        double dist = ZGeom::calGeodesic(*mMesh, hIdx, vFreeIdx[i]);
		vDist2Handle[i] = dist;
	});
	double distMax = *std::max_element(vDist2Handle.begin(), vDist2Handle.end());

	std::vector<ZGeom::Vec3d> vDeformedPos(freeVertCount);
	for (int i = 0; i < freeVertCount; ++i) {
		CVertex* pv = mMesh->vert(vFreeIdx[i]);
		vDeformedPos[i] = (ZGeom::Vec3d)pv->pos() + handleTrans * (1.0 - vDist2Handle[i]/distMax);
	}

    MeshCoordinates newCoord = mMesh->getVertCoordinates();
	ZGeom::VecNd& xNewCoord = newCoord.getXCoord();
	ZGeom::VecNd& yNewCoord = newCoord.getYCoord();
	ZGeom::VecNd& zNewCoord = newCoord.getZCoord();
	
	for (int i = 0; i < freeVertCount; ++i) {
		xNewCoord[vFreeIdx[i]] = vDeformedPos[i].x;
		yNewCoord[vFreeIdx[i]] = vDeformedPos[i].y;
		zNewCoord[vFreeIdx[i]] = vDeformedPos[i].z;
	}
	mMesh->setVertCoordinates(newCoord);

	timer.stopTimer("Deformation time: ");
}

void ShapeEditor::deformLaplacian()
{
	CStopWatch timer;	
	timer.startTimer();

	int anchorCount(0);
	std::vector<int> anchorIndex;
	std::vector<ZGeom::Vec3d> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();
    MeshCoordinates oldCoord = mMesh->getVertCoordinates();
	const MeshLaplacian& ml = mMeshHelper->getMeshLaplacian(SymCot);
	const ZGeom::SparseMatrix<double>& matLs = ml.getLS();
	ZGeom::SparseMatVecMultiplier mulLs(matLs, true);	
	VecNd diffCoord[3];
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(vertCount);
		mulLs.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
	}

	vector<VecNd> solveRHS(3);
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount, 0);
		solveRHS[i].copyElements(diffCoord[i], 0);
		for (int l = 0; l < anchorCount; ++l) {
			solveRHS[i][vertCount + l] = anchorWeight * anchorPos[l][i];
		}
	}
	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount, vertCount);
	matOptS.copyElements(matLs);
	for (int a = 0; a < anchorCount; ++a)
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	std::vector<ZGeom::VecNd> vls;
	solveSparseMultiColumn(g_engineWrapper, matOptS, solveRHS, vls);
	timer.stopTimer("Deformation time: ");
	
	MeshCoordinates newCoord(vertCount, vls);
	mMesh->setVertCoordinates(newCoord);
}

void ShapeEditor::deformLaplacian_v2()
{
	CStopWatch timer;	
	timer.startTimer();

	int anchorCount(0);
	std::vector<int> anchorIndex;
	std::vector<ZGeom::Vec3d> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();

	const MeshCoordinates& oldMeshCoord = this->getOldMeshCoord();

	const ZGeom::SparseMatrix<double>& matLs = mMeshHelper->getMeshLaplacian(SymCot).getLS();

	/* | A    B | | d |    | O  |
	   |        | |   | =  |    |
	   | B^T  O | | r |    | d' | */

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount, 0);
		for (int l = 0; l < anchorCount; ++l) {
			solveRHS[i][vertCount + l] = anchorPos[l][i] - oldMeshCoord[anchorIndex[l]][i];
		}
	}
	g_engineWrapper.addColVec(solveRHS[0].c_ptr(), vertCount + anchorCount, "dcx");
	g_engineWrapper.addColVec(solveRHS[1].c_ptr(), vertCount + anchorCount, "dcy");
	g_engineWrapper.addColVec(solveRHS[2].c_ptr(), vertCount + anchorCount, "dcz");

	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount, vertCount + anchorCount);
	matOptS.copyElements(matLs);
	for (int a = 0; a < anchorCount; ++a) {
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, 1.);
		matOptS.insertElem(anchorIndex[a] + 1, vertCount + a + 1, 1.);
	}

	if (!matOptS.testSymmetric()) {
		std::cout << "ERROR: MatOptS is not symmetric!" << std::endl;
		return;
	}

	g_engineWrapper.addSparseMat(matOptS, "matOptS");

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	g_engineWrapper.eval("lsx=matOptS\\dcx;");
	g_engineWrapper.eval("lsy=matOptS\\dcy;");
	g_engineWrapper.eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = g_engineWrapper.getDblVariablePtr("lsx");
	double *lsy = g_engineWrapper.getDblVariablePtr("lsy");
	double *lsz = g_engineWrapper.getDblVariablePtr("lsz");

	MeshCoordinates newCoord(oldMeshCoord);
	newCoord.add(lsx, lsy, lsz);
	mMesh->setVertCoordinates(newCoord);
}

void ShapeEditor::deformBiLaplacian()
{
	CStopWatch timer;	
	timer.startTimer();

	int anchorCount(0);
	std::vector<int> anchorIndex;
	std::vector<ZGeom::Vec3d> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();
    MeshCoordinates oldCoord = mMesh->getVertCoordinates();
	const ZGeom::SparseMatrix<double>& matLs = mMeshHelper->getMeshLaplacian(SymCot).getLS();
	ZGeom::SparseMatrix<double> matBiL;
	ZGeom::mulMatMat(matLs, matLs, matBiL);	
	ZGeom::SparseMatVecMultiplier mulBiLs(matBiL, true);	
	ZGeom::VecNd diffCoord[3];
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(vertCount);
		mulBiLs.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
	}
	
	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount, 0);
		solveRHS[i].copyElements(diffCoord[i], 0);
		for (int l = 0; l < anchorCount; ++l) {
			solveRHS[i][vertCount + l] = anchorWeight * anchorPos[l][i];
		}
	}
	g_engineWrapper.addColVec(solveRHS[0].c_ptr(), vertCount + anchorCount, "dcx");
	g_engineWrapper.addColVec(solveRHS[1].c_ptr(), vertCount + anchorCount, "dcy");
	g_engineWrapper.addColVec(solveRHS[2].c_ptr(), vertCount + anchorCount, "dcz");

	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount, vertCount);
	matOptS.copyElements(matBiL);
	for (int a = 0; a < anchorCount; ++a) 
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);
	g_engineWrapper.addSparseMat(matOptS, "matOptS");

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	g_engineWrapper.eval("lsx=matOptS\\dcx;");
	g_engineWrapper.eval("lsy=matOptS\\dcy;");
	g_engineWrapper.eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = g_engineWrapper.getDblVariablePtr("lsx");
	double *lsy = g_engineWrapper.getDblVariablePtr("lsy");
	double *lsz = g_engineWrapper.getDblVariablePtr("lsz");

	MeshCoordinates newCoord(vertCount, lsx, lsy, lsz);
	mMesh->setVertCoordinates(newCoord);
}

void ShapeEditor::deformMixedLaplacian(double ks, double kb)
{
	CStopWatch timer;	
	timer.startTimer();

	const int vertCount = mMesh->vertCount();
	const int anchorWeight = 1.0;

	int anchorCount(0);
	std::vector<int> anchorIndex;
	std::vector<ZGeom::Vec3d> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}
	
	std::vector<int> fixedVerts;
#if 0
	std::set<int> freeVerts;
	for (int i = 0; i < anchorCount; ++i) {
		std::set<int> sFree;
		mMesh->vertRingNeighborVerts(anchorIndex[i], 5, sFree, true);
		for (int v : sFree) freeVerts.insert(v);
	}
	for (int i = 0; i < vertCount; ++i) {
		if (freeVerts.find(i) == freeVerts.end())
			fixedVerts.push_back(i);
	}
#endif
	const int fixedCount = fixedVerts.size();

    MeshCoordinates oldCoord = mMesh->getVertCoordinates();
	g_engineWrapper.addColVec(oldCoord.getXCoord(), "ecx");
	g_engineWrapper.addColVec(oldCoord.getYCoord(), "ecy");
	g_engineWrapper.addColVec(oldCoord.getZCoord(), "ecz");

	const ZGeom::SparseMatrix<double>& matLs = mMeshHelper->getMeshLaplacian(SymCot).getLS();
	ZGeom::SparseMatrix<double> matBiL;
	ZGeom::mulMatMat(matLs, matLs, matBiL);
	g_engineWrapper.addSparseMat(matLs, "matL");
	g_engineWrapper.addSparseMat(matBiL, "matBiL");
	
	ZGeom::SparseMatrix<double> matL1(matLs);
	matL1.scale(-ks);
	ZGeom::SparseMatrix<double> matMixedL;
	ZGeom::addMatMat(matL1, matBiL, kb, matMixedL);	//matMixedL = -ks*matL + kb * matBiL
	g_engineWrapper.addSparseMat(matMixedL, "matMixedL");

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount + fixedCount, 0);
		for (int l = 0; l < anchorCount; ++l) {
			const ZGeom::Vec3d& oldPos = mMesh->vertPos(anchorIndex[l]);
			solveRHS[i][vertCount + l] = anchorWeight * (anchorPos[l][i] - oldPos[i]);
		}
	}
	g_engineWrapper.addColVec(solveRHS[0], "dcx");
	g_engineWrapper.addColVec(solveRHS[1], "dcy");
	g_engineWrapper.addColVec(solveRHS[2], "dcz");

	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount + fixedCount, vertCount);
	matOptS.copyElements(matMixedL);
	for (int a = 0; a < anchorCount; ++a) 
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);
	for (int a = 0; a < fixedCount; ++a)
		matOptS.insertElem(vertCount + anchorCount + a + 1, fixedVerts[a] + 1, anchorWeight);
	g_engineWrapper.addSparseMat(matOptS, "matOptS");

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	g_engineWrapper.eval("lsx=matOptS\\dcx;");
	g_engineWrapper.eval("lsy=matOptS\\dcy;");
	g_engineWrapper.eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = g_engineWrapper.getDblVariablePtr("lsx");
	double *lsy = g_engineWrapper.getDblVariablePtr("lsy");
	double *lsz = g_engineWrapper.getDblVariablePtr("lsz");

	MeshCoordinates newCoord(oldCoord);
	newCoord.add(lsx, lsy, lsz);

	mMesh->setVertCoordinates(newCoord);
}

void ShapeEditor::meanCurvatureFlow( double tMultiplier, int nRepeat /*= 1*/ )
{
	const int vertCount = mMesh->vertCount();
	double oriVol = mMesh->calVolume();

    ZGeom::SparseSymMatVecSolver heat_solver;
    computeHeatDiffuseMatrix(*mMesh, tMultiplier, heat_solver);
	const MeshLaplacian& laplacian = mMeshHelper->getMeshLaplacian(CotFormula);
	const ZGeom::SparseMatrix<double>& matW = laplacian.getW();
	const ZGeom::SparseMatrix<double>& matLc = laplacian.getLS();	// negative
	ZGeom::SparseMatVecMultiplier mulW(matW, true);
    
    MeshCoordinates oldCoord(vertCount);
    MeshCoordinates newCoord = mMesh->getVertCoordinates();	
	for (int n = 0; n < nRepeat; ++n) {
		for (int a = 0; a < 3; ++a) 
			mulW.mul(newCoord.getCoordFunc(a), oldCoord.getCoordFunc(a));				
		for (int a = 0; a < 3; ++a) 
            heat_solver.solve(oldCoord.getCoordFunc(a), newCoord.getCoordFunc(a));
	}

	mMesh->setVertCoordinates(newCoord);	
//	mMesh->scaleAreaToVertexNum();
// 	double newVol = mMesh->calVolume();
// 	double volPreserveScale = std::pow(oriVol/newVol, 1.0/3.0);
// 	mMesh->scaleAndTranslate(Vector3D(0,0,0), volPreserveScale);
}

void ShapeEditor::computeApproximations( const std::vector<ZGeom::VecNd>& vAtoms, 
										 ZGeom::SparseCoding* vApproxCoeff[3], 
										 int nReconstruct, 
										 std::vector<MeshCoordinates>& continuousCoords, 
										 MeshCoordinates& finalCoord )
{
	const int vertCount = mMesh->vertCount();
	const ZGeom::SparseCoding &vApproxX = *vApproxCoeff[0], &vApproxY = *vApproxCoeff[1], &vApproxZ = *vApproxCoeff[2];

	if (nReconstruct > vApproxX.size()) nReconstruct = (int)vApproxX.size();
	continuousCoords.resize(nReconstruct);
	finalCoord.resize(vertCount);

	for (int i = 0; i < nReconstruct; ++i) {
		finalCoord.getXCoord() += vApproxX[i].coeff() * vAtoms[vApproxX[i].index()];
		finalCoord.getYCoord() += vApproxY[i].coeff() * vAtoms[vApproxY[i].index()];
		finalCoord.getZCoord() += vApproxZ[i].coeff() * vAtoms[vApproxZ[i].index()];

		continuousCoords[i] = finalCoord;
	}
}

void ShapeEditor::updateEditBasis( const std::vector<ZGeom::VecNd>& vAtoms, const std::vector<int>& vSelectedIdx )
{
	mEditBasis.clear();
	for (int idx : vSelectedIdx) mEditBasis.push_back(vAtoms[idx]);
}

void ShapeEditor::evaluateApproximation( const MeshCoordinates& newCoord, const std::string leadText )
{
	const MeshCoordinates& oldCoord = getOldMeshCoord();
	const int vertCount = newCoord.size();
	std::cout << "** Evaluate " << leadText << "\n";
	double dif1 = oldCoord.difference(newCoord);
	MeshCoordinates lCoord1, lCoord2;
	computeGeometricLaplacianCoordinate(*mMesh, oldCoord, lCoord1);
	computeGeometricLaplacianCoordinate(*mMesh, newCoord, lCoord2);
	double dif2 = lCoord1.difference(lCoord2);
	std::cout << "-- Avg position error: " << dif1 / vertCount << '\n';
	std::cout << "-- Avg geometric error:    " << (dif1 + dif2)/2.0/vertCount << '\n';
}

void ShapeEditor::runTests()
{
    //testSurfaceArea();
	//testSparseCompression();		
	//testArtificialShapeMCA();
	//testArtificailShapeMCA2();
	//testArtificialShapeMCA3();
	//testDictionaryForDecomposition();
	//testSparseFeatureFinding();
    //testSparseInpainting();
    //testDictionaryCoherence();
    //testDenoisingDLRS();
    //testSparseDecomposition2();
    //testWaveletAnalysis();
    //testWaveletComputation();
    //fillHole();
    
    //bool skipExternalBoundary = false;
    //fillHoles(skipExternalBoundary);

//     mMesh = new CMesh;
//     mMesh->load("../../Data/favorite/eight_800.obj");
//    testSurfaceInpainting();    
}

//// Test partitioned approximation with graph Laplacian ////
//
void ShapeEditor::testSparseCompression()
{
	revertCoordinates();
	mStoredCoordinates.resize(5);
	
	std::cout << "\n======== Starting Sparse Compression test) ========\n";	

	const int totalVertCount = mMesh->vertCount();
	int eigenCount = -1;		// -1 means vertCount-1	
	int maxPatchSize = 1000;	// -1 means no segmentation
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	setStoredCoordinates(oldMeshCoord, 0);
	MeshCoordinates oldGeoCoord;
	computeGeometricLaplacianCoordinate(*mMesh, oldMeshCoord, oldGeoCoord);

	mShapeApprox.init(mMesh);
	mShapeApprox.doSegmentation(maxPatchSize);
	mSegmentPalette.generatePalette(mShapeApprox.partitionCount());	
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorSigAttr(StrAttrColorPartitions).attrValue().getColors();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
    emit meshSignatureAdded();
	mShapeApprox.doEigenDecomposition(Umbrella, eigenCount);	
	
	/* setup vectors of compress ratios and initialize vProgressiveCoords */
	const double maxCodingRatio = 1.0;
	double ratioToTest[] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5/*, 0.65, 0.8*/};
	std::vector<double> vCompressRatio(std::begin(ratioToTest), std::end(ratioToTest));

	std::vector<MeshCoordinates> vProgressiveCoords(vCompressRatio.size());
	std::ofstream ofs("output/approx_errors.txt");	
	for (int i = 0; i < vCompressRatio.size(); ++i) 
		ofs << vCompressRatio[i] << ((i < vCompressRatio.size()-1) ? '\t' : '\n');

	std::function<void(double,SparseApproxMethod,bool)> compressAndEvaluate = [&](double maxRatio, SparseApproxMethod method, bool useCompressionRatio) 
	{
		if (useCompressionRatio) {
			mShapeApprox.findSparseRepresentationByCompressionRatio(method, maxRatio);	
		} else {
			mShapeApprox.findSparseRepresentationByBasisRatio(method, maxRatio);
		}

		int ratiosCount = vCompressRatio.size();
		for (int i = 0; i < ratiosCount; ++i) {
			const double ratio = vCompressRatio[i];
			assert(ratio <= maxRatio);
			if (useCompressionRatio) {
				mShapeApprox.doSparseReconstructionByCompressionRatio(ratio, vProgressiveCoords[i]);
			} else {
				mShapeApprox.doSparseReconstructionByBasisRatio(ratio, vProgressiveCoords[i]);
			}
		}

		MeshCoordinates newGeoCoord;
		for (int i = 0; i < ratiosCount; ++i) {			
			computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
			double meshError = oldMeshCoord.difference(vProgressiveCoords[i]);
			double geoError = oldGeoCoord.difference(newGeoCoord);
			double mixedError = (meshError + geoError) / (2*totalVertCount);
			ofs << mixedError << ((i < ratiosCount-1) ? '\t' : '\n');
		}
	};

	std::function<void(SparseApproxMethod)> pursuitAndEvaluate = [&](SparseApproxMethod method) 
	{
		int ratiosCount = vCompressRatio.size();
		MeshCoordinates newGeoCoord;

		for (int i = 0; i < ratiosCount; ++i) {
			mShapeApprox.findSparseRepresentationByCompressionRatio(method, vCompressRatio[i]);			
			mShapeApprox.doSparseReconstructionBySize(-1, vProgressiveCoords[i]);
			computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
			double meshError = oldMeshCoord.difference(vProgressiveCoords[i]);
			double geoError = oldGeoCoord.difference(newGeoCoord);
			double mixedError = (meshError + geoError) / (2*totalVertCount);
			ofs << mixedError << ((i < ratiosCount-1) ? '\t' : '\n');
		}
	};

	// 1. low-pass approximate, MHB
#if 1
	printBeginSeparator("Low-pass Approx, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	compressAndEvaluate(maxCodingRatio, ZGeom::SA_Truncation, false);

	evaluateApproximation(vProgressiveCoords.back(), "1");
	setStoredCoordinates(vProgressiveCoords.back(), 1);
	mContReconstructCoords[0] = vProgressiveCoords;
	emit approxStepsChanged(0, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	// 2. SOMP, MHB
#if 0
	printBeginSeparator("SMP, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	compressAndEvaluate(maxCodingRatio, SA_SMP, true);

	evaluateApproximation(vProgressiveCoords.back(), "2");
	setStoredCoordinates(vProgressiveCoords.back(), 2);
	mContReconstructCoords[1] = vProgressiveCoords;
	emit approxStepsChanged(1, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	//2.5 SOMP, MHB+Spikes
#if 0
	printBeginSeparator("SOMP, MHB-Spikes", '-');
	mShapeApprox.constructDictionaries(DT_FourierSpikes);
	pursuitAndEvaluate(SA_SOMP);

	evaluateApproximation(vProgressiveCoords.back(), "2.5");
	printEndSeparator('-', 40);
#endif

	// 3. SOMP, SGW
#if 0
	printBeginSeparator("SOMP, SGW", '-');
	mShapeApprox.constructDictionaries(DT_SGW4);
	pursuitAndEvaluate(SA_SOMP);

	evaluateApproximation(vProgressiveCoords.back(), "3");
	setStoredCoordinates(vProgressiveCoords.back(), 3);
	mContReconstructCoords[2] = vProgressiveCoords;
	emit approxStepsChanged(2, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	// 4. SOMP, SGW + MHB
#if 1
	printBeginSeparator("SOMP, SGW-MHB", '-');
	mShapeApprox.constructDictionaries(DT_SGW4MHB);
	pursuitAndEvaluate(ZGeom::SA_SOMP);
	
	evaluateApproximation(vProgressiveCoords.back(), "4");
	setStoredCoordinates(vProgressiveCoords.back(), 4);
	mContReconstructCoords[3] = vProgressiveCoords;
	emit approxStepsChanged(3, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	ofs.close();
	std::cout << '\n';
	changeCoordinates(2);
	printEndSeparator('=', 40);
}

//// Test shape decomposition via signal separation
//
void ShapeEditor::testSparseDecomposition()
{
	/************************************************************************/
	/* 1. Support decomposition w or w/o segmentation                       */
	/* 2. Save the multilevel decomposition magnitude as color signatures   */
	/* 3. Optimized signal separation against MHB and SGW dictionary        */
	/*                                                                      */
	/************************************************************************/
	 
	revertCoordinates();
	mStoredCoordinates.resize(5);
	std::cout << "\n======== Starting sparseDcompositionTest) ========\n";	

	const int totalVertCount = mMesh->vertCount();
	int eigenCount = min(500, totalVertCount-1);
	int maxPatchSize = -1;

	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();	
	setStoredCoordinates(oldMeshCoord, 0);	// save original coordinates
	mShapeApprox.init(mMesh);
	mShapeApprox.doSegmentation(maxPatchSize);
	mSegmentPalette.generatePalette(mShapeApprox.partitionCount());	
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorSigAttr(StrAttrColorPartitions).attrValue().getColors();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
    emit meshSignatureAdded();
	
	mShapeApprox.doEigenDecomposition(Umbrella, eigenCount);	

	printBeginSeparator("SOMP, SGW", '-');

	// SOMP, SGW
	const int codingSize = 100;
	DictionaryType dictType = DT_SGW4MHB;
	SparseApproxMethod saMethod = ZGeom::SA_SOMP;

	MeshCoordinates reconstructedMeshCoord;
	
	mShapeApprox.constructDictionaries(dictType); 
	mShapeApprox.findSparseRepresentationBySize(saMethod, codingSize);
	mShapeApprox.doSparseReconstructionBySize(-1, reconstructedMeshCoord);
	
	std::cout << "Reconstruction error: " << oldMeshCoord.difference(reconstructedMeshCoord) << '\n';
	setStoredCoordinates(reconstructedMeshCoord, 1);
	

	const ZGeom::Dictionary& dict = mShapeApprox.mSubMeshApprox[0].getDict();
	const std::vector<ZGeom::SparseCodingItem>* vCoeff[3] = {
		&mShapeApprox.mSubMeshApprox[0].getSparseCoding(0), 
		&mShapeApprox.mSubMeshApprox[0].getSparseCoding(1), 
		&mShapeApprox.mSubMeshApprox[0].getSparseCoding(2)
	};
	int dictSize = dict.size();
	int actualCodingSize = int(vCoeff[0]->size());
	MeshCoordinates coordMHB(totalVertCount), coordSGW(totalVertCount);
	int mhbAtomCount(0), sgwAtomCount(0);
	
	for (int i = 0; i < actualCodingSize; ++i) {
		for (int c = 0; c < 3; ++c) {
			const ZGeom::SparseCodingItem& sc = (*vCoeff[c])[i];
			if (sc.index() < dictSize - totalVertCount) {
				coordSGW.getCoordFunc(c) += sc.coeff() * dict[sc.index()];
				if (c == 0) sgwAtomCount++;
			} else {
				coordMHB.getCoordFunc(c) += sc.coeff() * dict[sc.index()];
				if (c == 0) mhbAtomCount++;
			}
		}
	}

	setStoredCoordinates(coordMHB, 2);
	setStoredCoordinates(coordSGW, 3);
	std::cout << "#MHB atoms selected: " << mhbAtomCount << '\n';
	std::cout << "#SGW atoms selected: " << sgwAtomCount << '\n';

	MeshFeatureList* vSparseFeatures = new MeshFeatureList;
	for (int i = 0; i < actualCodingSize; ++i) {
		const ZGeom::SparseCodingItem& sc0 = (*vCoeff[0])[i];
		const ZGeom::SparseCodingItem& sc1 = (*vCoeff[1])[i];
		const ZGeom::SparseCodingItem& sc2 = (*vCoeff[2])[i];
		if (sc0.index() >= dictSize - totalVertCount) continue;

		int scale = sc0.index() / totalVertCount, idx = sc0.index() % totalVertCount;
		double coef = std::sqrt(sc0.coeff()*sc0.coeff() + sc1.coeff()*sc1.coeff() + sc2.coeff()*sc2.coeff());
		
		vSparseFeatures->addFeature(new MeshFeature(idx, scale));
		vSparseFeatures->back()->m_scalar1 = 3./*coef*/;
	}
	//delete vSparseFeatures;
	mMesh->addAttrMeshFeatures(*vSparseFeatures, StrAttrFeatureSparseSGW);


	MeshCoordinates coordResidual = oldMeshCoord.substract(reconstructedMeshCoord);
	setStoredCoordinates(coordResidual, 4);

	std::cout << '\n';	
	printEndSeparator('=', 40);
	changeCoordinates(1);
}

void ShapeEditor::testSparseDecomposition2()
{
    revertCoordinates();
    std::cout << "\n======== Starting testSparseDecomposition2 ========\n";
    const int totalVertCount = mMesh->vertCount();
    const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
    const double bbDiag = mMesh->getBoundingBox().length() * 2;
    const MeshCoordinates& coordOld = getOldMeshCoord();
    const vector<VecNd> vOriginalCoords = coordOld.to3Vec();

    /* Do Eigendecomposition */
    std::cout << "==== Do Eigendecomposition ====\n";
    int eigenCount = -1;					// -1 means full decomposition
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, -1);

    /* Computer Dictionary	*/
    std::cout << "\n==== Compute Dictionaries ====\n";
    Dictionary dictSpikes, dictMHB, dictSGW, dictMixed;
    computeDictionary(DT_UNIT, es, dictSpikes);
    computeDictionary(DT_Fourier, es, dictMHB);
    computeDictionary(DT_SGW3, es, dictSGW);
    combineDictionary(dictMHB, dictSGW, dictMixed);
    
    /* Compute OMP Approximation */
    std::cout << "\n==== Compute OMP Approximation ====\n";
    ZGeom::SparseApproximationOptions approxOpts;
    approxOpts.mCodingSize = 100;
    approxOpts.mApproxMethod = ZGeom::SA_OMP;
    {
        vector<SparseCoding> vOMPCodings;
        vector<VecNd> vApproximatedCoords;
        const Dictionary& dict = dictMHB;
        multiChannelSparseApproximate(vOriginalCoords, dict, vOMPCodings, approxOpts);
        multiChannelSparseReconstruct(dict, vOMPCodings, vApproximatedCoords);
        MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
        std::cout << "- dictMHB OMP Reconstruction error: " << coordOld.difference(coordApproxOMP) << '\n';
        setStoredCoordinates(coordApproxOMP, 1);
    }
    {
        vector<SparseCoding> vOMPCodings;
        vector<VecNd> vApproximatedCoords;
        const Dictionary& dict = dictSGW;
        multiChannelSparseApproximate(vOriginalCoords, dict, vOMPCodings, approxOpts);
        multiChannelSparseReconstruct(dict, vOMPCodings, vApproximatedCoords);
        MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
        std::cout << "- dictSGW OMP Reconstruction error: " << coordOld.difference(coordApproxOMP) << '\n';
        setStoredCoordinates(coordApproxOMP, 2);
    }

    /* Compute MCA Decomposition */
    std::cout << "\n==== Compute MCA Decomposition ====\n";
    vector<const Dictionary*> vDicts{ &dictMHB, &dictSGW };
    ZGeom::MCAoptions mcaOpts;
    mcaOpts.nIter = 100;
    mcaOpts.threshMode = ZGeom::MCAoptions::HARD_THRESH;
    mcaOpts.threshStrategy = ZGeom::MCAoptions::MIN_OF_MAX;
    vector<SparseCoding> vMCACodings[3];
    vector<VecNd> vCartoon(3), vTexture(3);
    for (int c = 0; c < 3; ++c) {
        singleChannelMCA3(vOriginalCoords[c], vDicts, vMCACodings[c], &mcaOpts);
        vCartoon[c] = singleChannelSparseReconstruct(*vDicts[0], vMCACodings[c][0]);
        vTexture[c] = singleChannelSparseReconstruct(*vDicts[1], vMCACodings[c][1]);
    }
    MeshCoordinates mcCartoon(totalVertCount, vCartoon);
    MeshCoordinates mcTexture(totalVertCount, vTexture);
    MeshCoordinates mcApproxMCA = mcCartoon.add(mcTexture);
    std::cout << "- dictMHB and dictSGW MCA Reconstruction error: " << coordOld.difference(mcApproxMCA) << '\n';
    setStoredCoordinates(mcCartoon, 3);
    setStoredCoordinates(mcApproxMCA, 4);

    VecNd vNormSGW(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) vNormSGW[i] = mcTexture.getVertCoordinate(i).length();
    vector<ZGeom::Colorf> vColorsSGW(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) vColorsSGW[i].falseColor(vNormSGW[i] / (bbDiag * 0.02), 1.f, ZGeom::CM_COOL);
    addColorSignature("color_mca_sgw", vColorsSGW);

	printEndSeparator('=', 40);	
}

void ShapeEditor::testDictionaryForDecomposition()
{
	using ZGeom::combineDictionary;
	revertCoordinates();
	std::cout << "\n======== Starting testDictionaryForDecomposition ========\n";
	const int totalVertCount = mMesh->vertCount();
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	std::vector<ZGeom::VecNd> vOriginalCoords{ oldMeshCoord.getXCoord(), oldMeshCoord.getYCoord(), oldMeshCoord.getZCoord() };

	/* Compute Laplacian and eigendecomposition */
	std::cout << "==== Do Eigendecomposition ====\n";
	int eigenCount = -1;					// -1 means full decomposition
	MeshLaplacian graphLaplacian, cotLaplacian, aniso1Laplacian, symCotLaplacian;
	ZGeom::EigenSystem esGraph, esAniso, esCot, esSymCot;
	graphLaplacian.constructUmbrella(mMesh);
	graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esGraph);
	aniso1Laplacian.constructAnisotropic2(mMesh);
	aniso1Laplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esAniso);
	cotLaplacian.constructCotFormula(mMesh);
	cotLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esCot);
	symCotLaplacian.constructSymCot(mMesh);	
	symCotLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esSymCot);
	
	/* Computer Dictionary */
	std::cout << "\n==== Compute Dictionaries ====\n";
	std::vector<Dictionary> dict(10);
	computeDictionary(DT_Fourier, esGraph, dict[0]);
	computeDictionary(DT_SGW4, esGraph, dict[1]);
	combineDictionary(dict[0], dict[1], dict[2]);
	computeDictionary(DT_SGW4, esAniso, dict[3]);
	combineDictionary(dict[0], dict[3], dict[4]);
	computeDictionary(DT_SGW4, esCot, dict[5]);
	combineDictionary(dict[0], dict[5], dict[6]);
	computeDictionary(DT_Fourier, esSymCot, dict[7]);
	computeDictionary(DT_SGW4, esSymCot, dict[8]);
	combineDictionary(dict[7], dict[8], dict[9]);

	/* Compute OMP Approximation */
	std::cout << "\n==== Compute OMP Approximation ====\n";
	std::vector<int> vCodeNum{ 20, 50, 100 };
	std::vector<std::vector<double> > vErrors(dict.size());
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mCodingSize = vCodeNum.back();
	approxOpts.mApproxMethod = ZGeom::SA_OMP;

	std::cout << "- OMP reconstruction error: \n";
	for (int i = 0; i < (int)dict.size(); ++i) {
		if (dict[i].empty()) continue;
		std::vector<ZGeom::SparseCoding> vOMPCodings;
		std::vector<VecNd> vApproximatedCoords;
		for (int j = 0; j < (int)vCodeNum.size(); ++j) {
			approxOpts.mCodingSize = vCodeNum[j];
			multiChannelSparseApproximate(vOriginalCoords, dict[i], vOMPCodings, approxOpts);
			multiChannelSparseReconstruct(dict[i], vOMPCodings, vApproximatedCoords);
			MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
			vErrors[i].push_back(oldMeshCoord.difference(coordApproxOMP));
		}		
		std::cout << " dict " << i << ": ";
		for (double e : vErrors[i]) std::cout << e << ", ";
		std::cout << std::endl;
	}
	
	std::cout << '\n';
	printEndSeparator('=', 40);
}

void ShapeEditor::testArtificialShapeMCA()
{
	revertCoordinates();
	std::cout << "\n======== Starting testArtificialShapeMCA ========\n";
	const int totalVertCount = mMesh->vertCount();
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	std::vector<ZGeom::VecNd> vOriginalCoords{ oldMeshCoord.getXCoord(), oldMeshCoord.getYCoord(), oldMeshCoord.getZCoord() };

	/// Do Eigendecomposition
	//
	std::cout << "==== Do Eigendecomposition ====\n";
	int eigenCount = -1;					// -1 means full decomposition
	MeshLaplacian graphLaplacian;
	graphLaplacian.constructUmbrella(mMesh);
	ZGeom::EigenSystem es;
	graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, es);


	/// Computer Dictionary
	//
	std::cout << "\n==== Compute Dictionaries ====\n";
	Dictionary dict1, dict2;
	computeDictionary(DT_Fourier, es, dict1);
	computeDictionary(DT_SGW1, es, dict2);

	/// Compute OMP Approximation
	//
	std::cout << "\n==== Compute OMP Approximation ====\n";
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mCodingSize = 50;
	approxOpts.mApproxMethod = ZGeom::SA_SOMP;
	std::vector<ZGeom::SparseCoding> vOMPCodings;
	std::vector<ZGeom::VecNd> vApproximatedCoords;
	multiChannelSparseApproximate(vOriginalCoords, dict1, vOMPCodings, approxOpts);
	multiChannelSparseReconstruct(dict1, vOMPCodings, vApproximatedCoords);
	MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
	setStoredCoordinates(coordApproxOMP, 0);	// set as base mesh coordinates
	changeCoordinates(0);
	std::cout << "- S-OMP Reconstruction error: " << oldMeshCoord.difference(coordApproxOMP) << '\n';

	/// Add SGW components
	//
	std::cout << "\n==== Mix with extra SGW components ====\n";
	std::vector<ZGeom::VecNd> vAlteredCoords = vApproximatedCoords;
	int nnz2 = 15;
	std::vector<ZGeom::SparseCoding> vAddedCoeff(3);
	std::mt19937 engine(0);
	std::uniform_int_distribution<int> distNZ(0, dict2.size() - 1);
	std::normal_distribution<double> distCoeff(0, bbDiag * 0.01);
	for (int i = 0; i < nnz2; ++i) {
		//int selectedNNZ = dict2.size() * (rand() / (double)RAND_MAX);
		int selectedNNZ = distNZ(engine);
		for (int c = 0; c < 3; ++c) {
			double coeff = distCoeff(engine);
			vAlteredCoords[c] += coeff * dict2[selectedNNZ];
			vAddedCoeff[c].addItem(selectedNNZ, coeff);
		}
	}
	MeshCoordinates coordAltered(totalVertCount, vAlteredCoords);
	setStoredCoordinates(coordAltered, 1);

	std::cout << "- #MHB: " << vOMPCodings[0].sortByCoeff().size() 
		      << "\t#SGW: " << vAddedCoeff[0].sortByCoeff().size() << '\n';
	std::cout << "-- MHB coefficients: ";
	for (auto p : vOMPCodings[0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- SGW coefficients: ";
	for (auto p : vAddedCoeff[0].getApproxItems()) std::cout << '(' << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	// compute and display color signatures
	const VecNd &vX0 = coordApproxOMP.getXCoord(), &vX1 = coordAltered.getXCoord();
	auto px1 = vX0.min_max_element(), px2 = vX1.min_max_element();
	//double x_min = min(px1.first, px2.first), x_max = max(px1.second, px2.second);
	double x_min = px1.first, x_max = px1.second;
	std::string colorStr0 = "color_x_coord_base";
	std::vector<ZGeom::Colorf> vColors0(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) vColors0[i].falseColor((vX0[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr0, vColors0);
	std::string colorStr1 = "color_x_coord_altered";
	std::vector<ZGeom::Colorf> vColors1(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) 
		vColors1[i].falseColor((vX1[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr1, vColors1);
	std::string colorStr2 = "color_x_coord_diff";
	std::vector<ZGeom::Colorf> vColors2(totalVertCount);
	std::vector<double> vDiff = (vX0 - vX1).toStdVector();
	for (auto& d : vDiff) d = std::fabs(d);
	double maxDiff = *std::max_element(vDiff.begin(), vDiff.end());
	for (int i = 0; i < totalVertCount; ++i)
		vColors2[i].falseColor(min(float(vDiff[i] / (bbDiag*0.02)), 1.f), 1.f, ZGeom::CM_COOL);
	addColorSignature(colorStr2, vColors2);

	/// Compute MCA decomposition
	//
	std::cout << "\n==== Compute MCA Decomposition ====\n";
	std::vector<const Dictionary*> vDicts{ &dict1, &dict2 };
	ZGeom::MCAoptions mcaOpts;
	mcaOpts.nIter = 50;
	mcaOpts.threshMode = ZGeom::MCAoptions::HARD_THRESH;
	mcaOpts.threshStrategy = ZGeom::MCAoptions::MIN_OF_MAX;
	std::vector<ZGeom::SparseCoding> vMCACodings[3];
	for (int c = 0; c < 3; ++c) {
		singleChannelMCA(vAlteredCoords[c], vDicts, vMCACodings[c], &mcaOpts);
	}
	std::vector<VecNd> vCartoon(3), vTexture(3);
	for (int c = 0; c < 3; ++c) {
		vCartoon[c] = singleChannelSparseReconstruct(dict1, vMCACodings[c][0]);
		vTexture[c] = singleChannelSparseReconstruct(dict2, vMCACodings[c][1]);
	}
	MeshCoordinates mcCartoon(totalVertCount, vCartoon);
	setStoredCoordinates(mcCartoon, 2);

	VecNd& vXC = vCartoon[0];
	VecNd& vXT = vTexture[0];
	vMCACodings[0][0].sortByCoeff(); vMCACodings[0][1].sortByCoeff();
	std::cout << "- #MCA1: " << vMCACodings[0][0].size() << "\t#MCA2: " << vMCACodings[0][1].size() << '\n';
	std::cout << "-- MCA1 coefficients: ";
	for (auto p : vMCACodings[0][0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- MCA2 coefficients: ";
	for (auto p : vMCACodings[0][1].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	std::string colorStrMCA1 = "color_x_coord_mca1", colorStrMCA2 = "color_x_coord_mca2";
	std::vector<ZGeom::Colorf> vColorsMCA1(totalVertCount), vColorsMCA2(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) {
		vColorsMCA1[i].falseColor((vXC[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
		vColorsMCA2[i].falseColor(min(float(fabs(vXT[i] / (0.02*bbDiag))), 1.f), 1.f, ZGeom::CM_COOL);
	}
	addColorSignature(colorStrMCA1, vColorsMCA1);
	addColorSignature(colorStrMCA2, vColorsMCA2);

	std::cout << '\n';
	printEndSeparator('=', 40);
}

void ShapeEditor::testArtificailShapeMCA2()
{
	revertCoordinates();
	std::cout << "\n======== Starting testArtificialShapeMCA ========\n";
	const int totalVertCount = mMesh->vertCount();
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& coordOld = getOldMeshCoord();
	std::vector<ZGeom::VecNd> vOriginalCoords{ coordOld.getXCoord(), coordOld.getYCoord(), coordOld.getZCoord() };

	/* Do Eigendecomposition */
	std::cout << "==== Do Eigendecomposition ====\n";
	int eigenCount = -1;					// -1 means full decomposition
	MeshLaplacian graphLaplacian;
	graphLaplacian.constructUmbrella(mMesh);
	ZGeom::EigenSystem es;
	graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, es);

	/* Computer Dictionary */
	std::cout << "\n==== Compute Dictionaries ====\n";
	Dictionary dict1, dict2;
	computeDictionary(DT_Fourier, es, dict1);
	double sgw_dict_time = time_call_sec([&](){ computeDictionary(DT_SGW1, es, dict2); });
	std::cout << "** sgw_dict_time: " << sgw_dict_time << '\n';
	

	/* Compute OMP Approximation */
	std::cout << "\n==== Compute OMP Approximation ====\n";
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mCodingSize = 50;
	approxOpts.mApproxMethod = ZGeom::SA_SOMP;
	std::vector<ZGeom::SparseCoding> vOMPCodings;
	std::vector<ZGeom::VecNd> vApproximatedCoords;
	multiChannelSparseApproximate(vOriginalCoords, dict1, vOMPCodings, approxOpts);
	multiChannelSparseReconstruct(dict1, vOMPCodings, vApproximatedCoords);
	MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
	setStoredCoordinates(coordOld, 0);	// set as base mesh coordinates
	changeCoordinates(0);
	std::cout << "- S-OMP Reconstruction error: " << coordOld.difference(coordApproxOMP) << '\n';

	/* Add SGW components */
	std::cout << "\n==== Mix with extra SGW components ====\n";
	std::vector<ZGeom::VecNd> vAlteredCoords = vOriginalCoords;
	int nnz2 = 15;
	std::vector<ZGeom::SparseCoding> vAddedCoeff(3);
	std::mt19937 engine(0);
	std::uniform_int_distribution<int> distNZ(0, dict2.size() - 1);
	std::normal_distribution<double> distCoeff(0, bbDiag * 0.015);
	for (int i = 0; i < nnz2; ++i) {
		int selectedNNZ = distNZ(engine);
		for (int c = 0; c < 3; ++c) {
			double coeff = distCoeff(engine);
			vAlteredCoords[c] += coeff * dict2[selectedNNZ];
			vAddedCoeff[c].addItem(selectedNNZ, coeff);
		}
	}
	MeshCoordinates coordAltered(totalVertCount, vAlteredCoords);
	setStoredCoordinates(coordAltered, 1);

	std::cout << "- #MHB: " << vOMPCodings[0].size() << "\t#SGW: " << vAddedCoeff[0].size() << '\n';
	vOMPCodings[0].sortByCoeff(); vAddedCoeff[0].sortByCoeff();
	std::cout << "-- MHB coefficients: ";
	for (auto p : vOMPCodings[0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- SGW coefficients: ";
	for (auto p : vAddedCoeff[0].getApproxItems()) std::cout << '(' << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	// compute and display color signatures
	const VecNd& vX0 = coordOld.getXCoord(), &vX1 = coordAltered.getXCoord();
	auto px1 = vX0.min_max_element(), px2 = vX1.min_max_element();
	double x_min = px1.first, x_max = px1.second;
	std::string colorStr0 = "color_x_coord_base";
	std::vector<ZGeom::Colorf> vColors0(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) vColors0[i].falseColor((vX0[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr0, vColors0);
	std::string colorStr1 = "color_x_coord_altered";
	std::vector<ZGeom::Colorf> vColors1(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i)
		vColors1[i].falseColor((vX1[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr1, vColors1);
	std::string colorStr2 = "color_x_coord_diff";
	std::vector<ZGeom::Colorf> vColors2(totalVertCount);
	std::vector<double> vDiff = (vX0 - vX1).toStdVector();
	for (auto& d : vDiff) d = std::fabs(d);
	double maxDiff = *std::max_element(vDiff.begin(), vDiff.end());
	for (int i = 0; i < totalVertCount; ++i)
		vColors2[i].falseColor(min(float(vDiff[i] / (bbDiag*0.02)), 1.f), 1.f, ZGeom::CM_COOL);
	addColorSignature(colorStr2, vColors2);

	/// Compute MCA decomposition
	//
	std::cout << "\n==== Compute MCA Decomposition ====\n";
	std::vector<const Dictionary*> vDicts{ &dict1, &dict2 };
	ZGeom::MCAoptions mcaOpts;
	mcaOpts.nIter = 50;
	mcaOpts.threshMode = ZGeom::MCAoptions::HARD_THRESH;
	mcaOpts.threshStrategy = ZGeom::MCAoptions::MIN_OF_MAX;
	std::vector<ZGeom::SparseCoding> vMCACodings[3];
	for (int c = 0; c < 3; ++c) {
		singleChannelMCA(vAlteredCoords[c], vDicts, vMCACodings[c], &mcaOpts);
	}
	std::vector<VecNd> vCartoon(3), vTexture(3);
	for (int c = 0; c < 3; ++c) {
		vCartoon[c] = singleChannelSparseReconstruct(dict1, vMCACodings[c][0]);
		vTexture[c] = singleChannelSparseReconstruct(dict2, vMCACodings[c][1]);
	}
	MeshCoordinates mcCartoon(totalVertCount, vCartoon);
	setStoredCoordinates(mcCartoon, 2);

	VecNd& vXC = vCartoon[0];
	VecNd& vXT = vTexture[0];
	std::cout << "- #MCA1: " << vMCACodings[0][0].sortByCoeff().size() 
		      << "\t#MCA2: " << vMCACodings[0][1].sortByCoeff().size() << '\n';
	std::cout << "-- MCA1 coefficients: ";
	for (auto p : vMCACodings[0][0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- MCA2 coefficients: ";
	for (auto p : vMCACodings[0][1].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	std::string colorStrMCA1 = "color_x_coord_mca1", colorStrMCA2 = "color_x_coord_mca2";
	std::vector<ZGeom::Colorf> vColorsMCA1(totalVertCount), vColorsMCA2(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) {
		vColorsMCA1[i].falseColor((vXC[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
		vColorsMCA2[i].falseColor(min(float(fabs(vXT[i] / (0.02*bbDiag))), 1.f), 1.f, ZGeom::CM_COOL);
	}
	addColorSignature(colorStrMCA1, vColorsMCA1);
	addColorSignature(colorStrMCA2, vColorsMCA2);

	/* add the origin of selected SGW as feature point */
	MeshFeatureList vSparseFeatures;
	for (auto c : vMCACodings[0][1].getApproxItems()) {
		int scale = c.index() / totalVertCount, idx = c.index() % totalVertCount;
		double coef = c.coeff();
		vSparseFeatures.addFeature(new MeshFeature(idx, scale));
		vSparseFeatures.back()->m_scalar1 = coef;
	}
	mMesh->addAttrMeshFeatures(vSparseFeatures, StrAttrFeatureMCA);

	std::cout << '\n';
	printEndSeparator('=', 40);
}

void ShapeEditor::testArtificialShapeMCA3()
{
	revertCoordinates();
	std::cout << "\n======== Starting testArtificialShapeMCA ========\n";
	const int totalVertCount = mMesh->vertCount();
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& coordOld = getOldMeshCoord();
	vector<VecNd> vOriginalCoords{ coordOld.getXCoord(), coordOld.getYCoord(), coordOld.getZCoord() };

	/* Do Eigendecomposition	*/
	std::cout << "==== Do Eigendecomposition ====\n";
	int eigenCount = -1;					// -1 means full decomposition
	MeshLaplacian graphLaplacian;
	graphLaplacian.constructUmbrella(mMesh);
	ZGeom::EigenSystem es;
	graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, es);


	/* Computer Dictionary	*/
	std::cout << "\n==== Compute Dictionaries ====\n";
	Dictionary dictMHB, dictSGW, dictHK, dictSpikes;
	computeDictionary(DT_UNIT, es, dictSpikes);
	computeDictionary(DT_Fourier, es, dictMHB);
	computeDictionary(DT_SGW1, es, dictSGW);
	//calHKDict(es, 5., dictHK);

	/* Compute OMP Approximation */
	std::cout << "\n==== Compute OMP Approximation ====\n";
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mCodingSize = 50;
	approxOpts.mApproxMethod = ZGeom::SA_OMP;
	std::vector<ZGeom::SparseCoding> vOMPCodings;
	std::vector<ZGeom::VecNd> vApproximatedCoords;
	multiChannelSparseApproximate(vOriginalCoords, dictMHB, vOMPCodings, approxOpts);
	multiChannelSparseReconstruct(dictMHB, vOMPCodings, vApproximatedCoords);
	MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
	setStoredCoordinates(coordOld, 0);	// set as base mesh coordinates
	changeCoordinates(0);
	std::cout << "- S-OMP Reconstruction error: " << coordOld.difference(coordApproxOMP) << '\n';

	/* Add random noise */	
	std::cout << "\n==== Mix with extra SGW components ====\n";
	std::vector<ZGeom::VecNd> vAlteredCoords = vOriginalCoords;

	int nnz2 = 20;
	std::vector<ZGeom::SparseCoding> vAddedCoeff(3);
	std::mt19937 engine(0);
	std::uniform_int_distribution<int> distNZ(0, dictSpikes.size() - 1);
	std::normal_distribution<double> distCoeff(0, bbDiag * 0.02);
	for (int i = 0; i < nnz2; ++i) {
		int selectedNNZ = distNZ(engine);
		for (int c = 0; c < 3; ++c) {
			double coeff = distCoeff(engine);
			vAlteredCoords[c] += coeff * dictSpikes[selectedNNZ];
			vAddedCoeff[c].addItem(selectedNNZ, coeff);
		}
	}

	std::cout << "- #MHB: " << vOMPCodings[0].size() << "\t#HK: " << vAddedCoeff[0].size() << '\n';
	vOMPCodings[0].sortByCoeff(); vAddedCoeff[0].sortByCoeff();
	std::cout << "-- MHB coefficients: ";
	for (auto p : vOMPCodings[0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- SGW coefficients: ";
	for (auto p : vAddedCoeff[0].getApproxItems()) std::cout << '(' << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	MeshCoordinates coordAltered(totalVertCount, vAlteredCoords);
	setStoredCoordinates(coordAltered, 1);


	/* compute and display color signatures */
	const VecNd& vX0 = coordOld.getXCoord(), &vX1 = coordAltered.getXCoord();
	auto px1 = vX0.min_max_element(), px2 = vX1.min_max_element();
	double x_min = px1.first, x_max = px1.second;
	std::string colorStr0 = "color_x_coord_base";
	std::vector<ZGeom::Colorf> vColors0(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) vColors0[i].falseColor((vX0[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr0, vColors0);
	std::string colorStr1 = "color_x_coord_altered";
	std::vector<ZGeom::Colorf> vColors1(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i)
		vColors1[i].falseColor((vX1[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr1, vColors1);
	std::string colorStr2 = "color_x_coord_diff";
	std::vector<ZGeom::Colorf> vColors2(totalVertCount);
	std::vector<double> vDiff = (vX0 - vX1).toStdVector();
	for (auto& d : vDiff) d = std::fabs(d);
	double maxDiff = *std::max_element(vDiff.begin(), vDiff.end());
	for (int i = 0; i < totalVertCount; ++i)
		vColors2[i].falseColor(min(float(vDiff[i] / (bbDiag*0.02)), 1.f), 1.f, ZGeom::CM_COOL);
	addColorSignature(colorStr2, vColors2);

	/* Compute MCA decomposition */
	std::cout << "\n==== Compute MCA Decomposition ====\n";
	std::vector<const Dictionary*> vDicts{ &dictMHB, &dictSGW };
	ZGeom::MCAoptions mcaOpts;
	mcaOpts.nIter = 50;
	mcaOpts.threshMode = ZGeom::MCAoptions::HARD_THRESH;
	mcaOpts.threshStrategy = ZGeom::MCAoptions::MEAN_OF_MAX;
	std::vector<ZGeom::SparseCoding> vMCACodings[3];
	std::vector<VecNd> vCartoon(3), vTexture(3);
	for (int c = 0; c < 3; ++c) {
		singleChannelMCA(vAlteredCoords[c], vDicts, vMCACodings[c], &mcaOpts);
		vCartoon[c] = singleChannelSparseReconstruct(*vDicts[0], vMCACodings[c][0]);
		vTexture[c] = singleChannelSparseReconstruct(*vDicts[1], vMCACodings[c][1]);
	}
	MeshCoordinates mcCartoon(totalVertCount, vCartoon);
	MeshCoordinates mcTexture(totalVertCount, vTexture);
	setStoredCoordinates(mcCartoon, 2);
	setStoredCoordinates(coordOld.substract(mcTexture), 3);

	VecNd& vXC = vCartoon[0];
	VecNd& vXT = vTexture[0];
	std::cout << "- #MCA1: " << vMCACodings[0][0].sortByCoeff().size()
			  << "\t#MCA2: " << vMCACodings[0][1].sortByCoeff().size() << '\n';
	std::cout << "-- MCA1 coefficients: ";
	for (auto p : vMCACodings[0][0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- MCA2 coefficients: ";
	for (auto p : vMCACodings[0][1].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	std::string colorStrMCA1 = "color_x_coord_mca1", colorStrMCA2 = "color_x_coord_mca2";
	std::vector<ZGeom::Colorf> vColorsMCA1(totalVertCount), vColorsMCA2(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) {
		vColorsMCA1[i].falseColor((vXC[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
		vColorsMCA2[i].falseColor(min(float(fabs(vXT[i] / (0.02*bbDiag))), 1.f), 1.f, ZGeom::CM_COOL);
	}
	addColorSignature(colorStrMCA1, vColorsMCA1);
	addColorSignature(colorStrMCA2, vColorsMCA2);

	/* add the origin of selected SGW as feature point */
	MeshFeatureList vSparseFeatures;
	for (auto c : vMCACodings[0][1].getApproxItems()) {
		int scale = c.index() / totalVertCount, idx = c.index() % totalVertCount;
		double coef = c.coeff();
		vSparseFeatures.addFeature(new MeshFeature(idx, scale));
		vSparseFeatures.back()->m_scalar1 = coef;
	}
	mMesh->addAttrMeshFeatures(vSparseFeatures, StrAttrFeatureMCA);

	printEndSeparator('=', 40);
}

//// Test approximation for feature analysis  ////
//
void ShapeEditor::testSparseFeatureFinding()
{
	revertCoordinates();
	CStopWatch timer;
	std::cout << "\n======== Starting testSparseFeatureFinding ========\n";

	const int totalVertCount = mMesh->vertCount();
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	vector<VecNd> vOriginalCoords{ oldMeshCoord.getXCoord(), oldMeshCoord.getYCoord(), oldMeshCoord.getZCoord() };
	setStoredCoordinates(oldMeshCoord, 0);

	/* Compute Laplacian and eigendecomposition */
	std::cout << "==== Do Eigendecomposition ====\n";
	int eigenCount = std::min(1500, totalVertCount - 1); // -1 means vertCount -1 
	eigenCount = -1;
	MeshLaplacian graphLaplacian, cotLaplacian, aniso1Laplacian, aniso2Laplacian, symCotLaplacian;
	graphLaplacian.constructUmbrella(mMesh);
	cotLaplacian.constructCotFormula(mMesh);
	aniso1Laplacian.constructAnisotropic1(mMesh);
	aniso2Laplacian.constructAnisotropic2(mMesh);
	symCotLaplacian.constructSymCot(mMesh);
	const ZGeom::EigenSystem& esGraph = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
	//const ZGeom::EigenSystem& esCot = mProcessor->prepareEigenSystem(cotLaplacian, eigenCount);
	const ZGeom::EigenSystem& esAniso = mMeshHelper->prepareEigenSystem(aniso1Laplacian, eigenCount);
	
	/* Computer Dictionary */
	std::cout << "\n==== Compute Dictionaries ====\n";
	std::vector<Dictionary> dict(10);
	computeDictionary(DT_Fourier, esGraph, dict[0]);
	double sgw_dict_time = time_call_sec([&](){ computeDictionary(DT_SGW4, esGraph, dict[1]); });
	std::cout << "** sgw_dict_time (s): " << sgw_dict_time << '\n';
	combineDictionary(dict[0], dict[1], dict[2]);
	//computeDictionary(DT_SGW4, esAniso, dict[3]);
	//combineDictionary(dict[0], dict[3], dict[4]);
	computeDictionary(DT_FourierSpikes, esGraph, dict[5]);

	auto sparseApproximate = [&](
		const MeshCoordinates& coordOld,
		const Dictionary& sparseDict,
		const ZGeom::SparseApproximationOptions& opts,
		vector<SparseCoding>& vCodings,
		MeshCoordinates& coordApprox) 
	{
		vector<VecNd> vOldCoords = oldMeshCoord.to3Vec();
		vector<VecNd> vApproximatedCoords;
		multiChannelSparseApproximate(vOldCoords, sparseDict, vCodings, opts);
		multiChannelSparseReconstruct(sparseDict, vCodings, vApproximatedCoords);
		coordApprox = MeshCoordinates(totalVertCount, vApproximatedCoords);
		std::cout << "-- reconstruction error: " << coordOld.difference(coordApprox) << '\n';
	};

	auto addFeatures = [&](
		const ZGeom::SparseCoding& sc, 
		int numGlobalAtom, 
		std::string featureName) 
	{
		MeshFeatureList vSparseFeatures;
		for (auto c : sc.getApproxItems()) {
			int scale = (c.index() - numGlobalAtom) / totalVertCount,
				idx = (c.index() - numGlobalAtom) % totalVertCount;
			if (scale > 0) vSparseFeatures.addFeature(new MeshFeature(idx, scale));
		}
		mMesh->addAttrMeshFeatures(vSparseFeatures, featureName);
		std::cout << "-- detected " << featureName << " (" << vSparseFeatures.size() << "): ";
		for (auto f : vSparseFeatures.getFeatureVector()) std::cout << f->m_index << ", ";
		std::cout << "\n\n";
	};
	
	enum AtomKeep { Before, After };
	auto separateReconstruct = [&](
		const Dictionary& sparseDict, 
		const vector<SparseCoding>& vCodings,
		int numAtomKeep, AtomKeep whichToKeep)
		-> MeshCoordinates
	{
		vector<SparseCoding> vTruncatedCodings(3);
		for (int l = 0; l < 3; ++l) {
			if (whichToKeep == Before) {
				for (auto c : vCodings[l].getApproxItems())
					if (c.index() < numAtomKeep) vTruncatedCodings[l].addItem(c);
			}
			else {
				for (auto c : vCodings[l].getApproxItems())
					if (c.index() >= numAtomKeep) vTruncatedCodings[l].addItem(c);
			}
		}
		std::vector<VecNd> vApproximatedCoords;
		multiChannelSparseReconstruct(sparseDict, vTruncatedCodings, vApproximatedCoords);
		return MeshCoordinates(totalVertCount, vApproximatedCoords);
	};

	auto visualizeCoordDiff = [&](const MeshCoordinates& coordDiff, std::string colorSigNuame) 
	{
		vector<ZGeom::Colorf> vColors(totalVertCount);
		double diffMax = bbDiag * 0.1;
		for (int i = 0; i < totalVertCount; ++i) {
			double diff = fabs(coordDiff.getVertCoordinate(i).length());
			vColors[i].falseColor(diff / diffMax, ZGeom::CM_JET);
		}
		addColorSignature(colorSigNuame, vColors);
	};

	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mApproxMethod = ZGeom::SA_SOMP;
	approxOpts.mCodingSize = 30;
	vector<SparseCoding> vApproxCodings[10];
	MeshCoordinates coordApprox[10];
	int coordCount = 1;

	{
		sparseApproximate(oldMeshCoord, dict[1], approxOpts, vApproxCodings[0], coordApprox[0]);
		setStoredCoordinates(coordApprox[0], coordCount++);
		addFeatures(vApproxCodings[0][0], 0, "mesh_feature_SOMP_dict1");
	}
	{
		sparseApproximate(oldMeshCoord, dict[2], approxOpts, vApproxCodings[1], coordApprox[1]);
		auto vInitCodingFront = multiExtractFront(vApproxCodings[1], dict[0].size());
		std::cout << "- initial #dict1: " << vInitCodingFront[0].size() << "\n";
		setStoredCoordinates(coordApprox[1], coordCount++);
		addFeatures(vApproxCodings[1][0], dict[0].size(), "mesh_feature_SOMP_dict2");
		MeshCoordinates coordGlobal = separateReconstruct(dict[2], vApproxCodings[1], dict[0].size(), Before);
		setStoredCoordinates(coordGlobal, coordCount++);
		visualizeCoordDiff(oldMeshCoord.substract(coordGlobal), "color_coord_dict2_diff");
		visualizeCoordDiff(separateReconstruct(dict[2], vApproxCodings[1], dict[0].size(), After), "color_coord_dict2_sgw");
	}
// 	{
// 		sparseApproximate(oldMeshCoord, dict[4], approxOpts, vApproxCodings[2], coordApprox[2]);
// 		setStoredCoordinates(coordApprox[2], coordCount++);
// 		addFeatures(vApproxCodings[2][0], dict[0].size(), "mesh_feature_SOMP3");
// 		MeshCoordinates coordGlobal = separateReconstruct(dict[4], vApproxCodings[2], dict[0].size(), Before);
// 		setStoredCoordinates(coordGlobal, coordCount++);
// 		visualizeCoordDiff(oldMeshCoord.substract(coordGlobal), "color_coord_diff2");
// 		visualizeCoordDiff(separateReconstruct(dict[4], vApproxCodings[2], dict[0].size(), After), "color_coord_diff_sgw2");
// 	}

	{
		vector< vector<SparseCoding> > vSeparateCodings;
		const Dictionary &dictMCA1 = dict[0], &dictMCA2 = dict[1];

		multiDictSparseDecomposition(
			oldMeshCoord,
			vector < const Dictionary* > {&dictMCA1, &dictMCA2},
			vector < int > {20, 40},
			vSeparateCodings);
		vector<VecNd> vCoordDict1, vCoordDict2;
		multiChannelSparseReconstruct(dictMCA1, vSeparateCodings[0], vCoordDict1);
		multiChannelSparseReconstruct(dictMCA2, vSeparateCodings[1], vCoordDict2);
		MeshCoordinates coordMCA1(totalVertCount, vCoordDict1);
		setStoredCoordinates(coordMCA1, coordCount++);
		visualizeCoordDiff(oldMeshCoord.substract(coordMCA1), "color_coord_MCA1_diff");
		visualizeCoordDiff(MeshCoordinates(totalVertCount, vCoordDict2), "color_coord_MCA2");
		addFeatures(vSeparateCodings[1][0], 0, "mesh_feature_MCA2");
	}
	
//  sparseApproximate(oldMeshCoord, dict[5], approxOpts, vApproxCodings[2], coordApprox[2]);
//  setStoredCoordinates(coordApprox[2], coordCount++);
//  addFeatures(vApproxCodings[2][0], dict[0].size(), "mesh_feature_SOMP3");

	changeCoordinates(1);
	printEndSeparator('=', 40);
}

void ShapeEditor::testSparseInpainting()
{
	// 1. Prepare dictionary
	// 2. Prepare coordinate data
	// 3. Perform inpainting 
	// 4. Visualize results

    revertCoordinates();
    std::cout << "\n======== Starting testSparseInpainting ========\n";
    const int totalVertCount = mMesh->vertCount();
    const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
    const double bbDiag = mMesh->getBoundingBox().length() * 2;
    const MeshCoordinates& coordOld = getOldMeshCoord();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();
    std::mt19937 pnEngine(0);
    std::uniform_int_distribution<int> distBlank(0, totalVertCount - 1);

    std::cout << "==== Do Eigendecomposition ====\n";
    int eigenCount = -1;					// -1 means full decomposition
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, -1);

    std::cout << "\n==== Compute Dictionaries ====\n";
    Dictionary dictMHB, dictSGW, dictMixed;
    computeDictionary(DT_Fourier, es, dictMHB);
//  computeDictionary(DT_SGW3, es, dictSGW);
//  combineDictionary(dictMHB, dictSGW, dictMixed);

    int countCoord = 1;
    auto inpaintingRecover = [&](const vector<VecNd>& vecCoords, double ratioMissing, const Dictionary& dictInpaint, std::string featureName) {
        vector<int> vMask(totalVertCount, 1);
        for (int i = 0; i < totalVertCount * ratioMissing; ++i) vMask[distBlank(pnEngine)] = 0;
        MeshFeatureList vFeatureMissing;
        for (int idx = 0; idx < totalVertCount; ++idx) {
            if (vMask[idx] == 0)
                vFeatureMissing.addFeature(new MeshFeature(idx, 0));
        }
        if (featureName != "") {
            mMesh->addAttrMeshFeatures(vFeatureMissing, featureName);

            vector<ZGeom::Colorf> vColorMissing(totalVertCount);
            for (int idx = 0; idx < totalVertCount; ++idx) {
                if (vMask[idx] == 0) vColorMissing[idx] = ZGeom::ColorRed;
                else vColorMissing[idx] = ZGeom::ColorMesh1;
            }
            addColorSignature("color_missing_vertex_0." + boost::lexical_cast<std::string>((int)std::round(ratioMissing*100)), vColorMissing);
        }            

        vector<VecNd> vCoordRecovered(3);
        vector<SparseCoding> vCodingInpaint(3);
        for (int i = 0; i < 3; ++i)
            vCoordRecovered[i] = singleChannelSparseInpaint(vecCoords[i], vMask, dictInpaint, vCodingInpaint[i]);

        MeshCoordinates coordRecovered(totalVertCount, vCoordRecovered);
        setStoredCoordinates(coordRecovered, countCoord++);
        std::cout << "-- Inpainting error: " << coordOld.difference(coordRecovered) << "\n";
    };

    inpaintingRecover(vOriginalCoords, .25, dictMHB, "feature_missing_vertex_0.25");
    inpaintingRecover(vOriginalCoords, .5, dictMHB, "feature_missing_vertex_0.50");    
    inpaintingRecover(vOriginalCoords, .75, dictMHB, "feature_missing_vertex_0.75");

    printEndSeparator('=', 40);
}

void ShapeEditor::testDenoisingDLRS()
{
    printBeginSeparator("Starting testDenoisingDLRS", '=');
    const int totalVertCount = mMesh->vertCount();
    const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
    const double bbDiag = mMesh->getBoundingBox().length() * 2;
    revertCoordinates();
    const MeshCoordinates& coordOld = getOldMeshCoord();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();
    MeshCoordinates coordNoisy = getNoisyCoord(0.005);
    setStoredCoordinates(coordNoisy, 1);    

    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    vector<VecNd> vDenoisedCoord = DLRS(g_engineWrapper, graphLaplacian.getLS(), 0.8, coordNoisy.to3Vec());
    MeshCoordinates coordDenoised(totalVertCount, vDenoisedCoord);
    setStoredCoordinates(coordDenoised, 2);
    MeshCoordinates coordResidual = coordDenoised.substract(coordNoisy);
    vector<ZGeom::Colorf> vColorDiff = colorCoordDiff(coordResidual, 0.05*bbDiag, ZGeom::CM_JET);
    addColorSignature("color_diff_DLRS", vColorDiff);

    auto vNormals = ZGeom::getMeshVertNormals(*mMesh);
    VecNd vResB(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i)
        vResB[i] = coordResidual.getVertCoordinate(i).dot((ZGeom::Vec3d)vNormals[i]);
    vector<Colorf> vColorRes(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) vColorRes[i].falseColor(fabs(vResB[i] / (0.05*bbDiag)));
    addColorSignature("color_residual_DLRS", vColorRes);

    vector<std::pair<int, double> > vResIdx;
    for (int i = 0; i < totalVertCount; ++i) vResIdx.push_back(std::make_pair(i, fabs(vResB[i])));
    using std::pair;
    std::sort(vResIdx.begin(), vResIdx.end(), [](pair<int, double> p1, pair<int, double> p2){ return p1.second > p2.second; });
    MeshFeatureList vFeatureMaxRes;
    for (int i = 0; i < 30; ++i) vFeatureMaxRes.addFeature(new MeshFeature(vResIdx[i].first, vResIdx[i].second));
    mMesh->addAttrMeshFeatures(vFeatureMaxRes, "feature_res_max");

    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, -1);
    Dictionary dictMHB, dictSGW, dictMixed;
    computeDictionary(DT_Fourier, es, dictMHB);
    computeDictionary(DT_SGW3, es, dictSGW);    
    //combineDictionary(dictMHB, dictSGW, dictMixed);
    ZGeom::SparseApproximationOptions opts;
    opts.mCodingSize = 30;
    opts.mApproxMethod = ZGeom::SA_OMP;
    opts.mMatlabEngine = &g_engineWrapper;
    vector<VecNd> vApproximatedCoords;
    SparseCoding sc;
    ZGeom::singleChannelSparseApproximate(vResB, dictSGW, sc, opts);
    VecNd vResApprox = ZGeom::singleChannelSparseReconstruct(dictSGW, sc);
    sc.sortByCoeff();
    MeshFeatureList vFeatureOMP;
    for (auto f : sc.getApproxItems())
        vFeatureOMP.addFeature(new MeshFeature(f.index() % totalVertCount, f.coeff()));
    mMesh->addAttrMeshFeatures(vFeatureOMP, "feature_res_omp");    

    vector<Colorf> vColorResApprox(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) vColorResApprox[i].falseColor(fabs(vResApprox[i] / (0.05*bbDiag)));
    addColorSignature("color_residual_approx", vColorResApprox);

    printEndSeparator('=', 40);
}

void ShapeEditor::testDictionaryCoherence()
{
    printBeginSeparator("Starting testDictionaryCoherence", '=');

    const int totalVertCount = mMesh->vertCount();
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, -1);
    Dictionary dictMHB, dictSGW3, dictMixed, dictSGW1, dictSGW2, dictSGW3MHB, dictSpikes;
    computeDictionary(DT_UNIT, es, dictSpikes);
    computeDictionary(DT_Fourier, es, dictMHB);
    computeDictionary(DT_SGW1, es, dictSGW1);
    computeDictionary(DT_SGW3, es, dictSGW3);
    combineDictionary(dictMHB, dictSpikes, dictMixed);
    
    std::cout << "Coherence of dictSpikes: " << dictSpikes.calCoherence() << '\n';
    std::cout << "Coherence of dictMHB: " << dictMHB.calCoherence() << '\n';
    std::cout << "Coherence of dictSGW3: " << dictSGW3.calCoherence() << '\n';
    std::cout << "Coherence of dictMixed: " << dictMixed.calCoherence() << '\n';
//  double co, tco;
//  tco = time_call([&](){co = dictMixed.calCoherence2(); });
//  std::cout << "Time to calCoherence: " << tco / 1000. << "\n";
//  std::cout << "Coherence of dictMixed: " << co << '\n';
//  std::cout << "Mutual coherence of dictMHB & dictSpikes: " << dictMHB.mutualCoherence(dictSpikes) << "\n";
//  std::cout << "Mutual coherence of dictMHB & dictSGW1: " << dictMHB.mutualCoherence(dictSGW1) << "\n";
//  std::cout << "Mutual coherence of dictMHB & dictSGW3: " << dictMHB.mutualCoherence(dictSGW3) << "\n";

    printEndSeparator('=', 40);
}

void ShapeEditor::testWaveletAnalysis()
{
    using namespace concurrency;

    printBeginSeparator("Starting testWaveletAnalysis", '=');
    const int totalVertCount = mMesh->vertCount();
    const double bbDiag = mMesh->getBoundingBox().length() * 2;
    revertCoordinates();
    const MeshCoordinates& coordOld = getOldMeshCoord();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();
    MeshCoordinates coordNoisy = getNoisyCoord(0.005);
    setStoredCoordinates(coordNoisy, 1);

    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, -1);
    Dictionary dictMHB, dictSGW;
    computeDictionary(DT_Fourier, es, dictMHB);
    computeDictionary(DT_SGW4, es, dictSGW);

    int sgwScales = dictSGW.atomCount() / totalVertCount;
    vector<array<VecNd,3>> vTransforms(sgwScales);
    vector<VecNd> vCoeffNorm(sgwScales), vCoeffLaplace(sgwScales);
    {
        parallel_for(0, sgwScales, [&](int s) {
            for (int c = 0; c < 3; ++c) {
                const VecNd& vSig = coordOld.getCoordFunc(c);
                vTransforms[s][c].resize(totalVertCount);
                for (int i = 0; i < totalVertCount; ++i) {
                    const VecNd& atom = dictSGW[s*totalVertCount + i];
                    vTransforms[s][c][i] = vSig.dot(atom);
                }
            }
            vCoeffNorm[s].resize(totalVertCount);
            for (int i = 0; i < totalVertCount; ++i)
                vCoeffNorm[s][i] = Vec3d(vTransforms[s][0][i], vTransforms[s][1][i], vTransforms[s][2][i]).length();

            vCoeffLaplace[s] = ZGeom::mulMatVec<double>(graphLaplacian.getLS(), vCoeffNorm[s], true);
            for (double &d : vCoeffLaplace[s]) d = fabs(d);
        });
        for (int s = 0; s < sgwScales; ++s) {
            vector<Colorf> vColor(totalVertCount);
            const VecNd& vSignal = vCoeffNorm[s];
            double coeffMax = vSignal.max_element() + 1e-3;
            for (int i = 0; i < totalVertCount; ++i)
                vColor[i].falseColor(vSignal[i] / coeffMax);
            addColorSignature("color_ori_coord_sgw_l" + Int2String(s), vColor);
        }
    }

    {
        parallel_for(0, sgwScales, [&](int s) {
            for (int c = 0; c < 3; ++c) {
                const VecNd& vSig = coordNoisy.getCoordFunc(c);
                vTransforms[s][c].resize(totalVertCount);
                for (int i = 0; i < totalVertCount; ++i) {
                    const VecNd& atom = dictSGW[s*totalVertCount + i];
                    vTransforms[s][c][i] = vSig.dot(atom);
                }
            }
            vCoeffNorm[s].resize(totalVertCount);
            for (int i = 0; i < totalVertCount; ++i)
                vCoeffNorm[s][i] = Vec3d(vTransforms[s][0][i], vTransforms[s][1][i], vTransforms[s][2][i]).length();
        });
        for (int s = 0; s < sgwScales; ++s) {
            vector<Colorf> vColor(totalVertCount);
            const VecNd& vSignal = vCoeffNorm[s];
            double coeffMax = vSignal.max_element() + 1e-3;
            for (int i = 0; i < totalVertCount; ++i)
                vColor[i].falseColor(vSignal[i] / coeffMax);
            addColorSignature("color_noisy_coord_sgw_l" + Int2String(s), vColor);
        }
    }
    
    vector<VecNd> vDenoisedCoord = DLRS(g_engineWrapper, graphLaplacian.getLS(), 1.0, coordOld.to3Vec());
    MeshCoordinates coordDenoised(totalVertCount, vDenoisedCoord);
    MeshCoordinates coordResidual = coordDenoised.substract(coordOld);
    setStoredCoordinates(coordDenoised, 2);
    changeCoordinates(2);
    auto vNormals = ZGeom::getMeshVertNormals(*mMesh);
    VecNd vResB(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i)
        vResB[i] = coordResidual.getVertCoordinate(i).dot((ZGeom::Vec3d)vNormals[i]);
    vector<Colorf> vColorRes(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) vColorRes[i].falseColor(fabs(vResB[i] / (0.05*bbDiag)));
    addColorSignature("color_residual_DLRS", vColorRes);
    {
        VecNd vResBLaplace = ZGeom::mulMatVec<double>(graphLaplacian.getLS(), vResB, true);
        for (double &d : vResBLaplace) d = fabs(d);
        double resLapMax = vResBLaplace.max_element();
        vector<Colorf> vColorResLaplace(totalVertCount);
        for (int i = 0; i < totalVertCount; ++i) vColorResLaplace[i].falseColor(vResBLaplace[i] / resLapMax);
        addColorSignature("color_residual_DLRS_laplace", vColorResLaplace);
    }
     
    printEndSeparator('=', 40);
}

void ShapeEditor::testWaveletComputation()
{
    const int totalVertCount = mMesh->vertCount();
    const double bbDiag = mMesh->getBoundingBox().length() * 2;
    revertCoordinates();
    const MeshCoordinates& coordOld = getOldMeshCoord();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();

    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);

    CStopWatch timer;
    timer.startTimer();
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, 500);
    timer.stopTimer("Decomposition time: ");

    Dictionary dictMHB, dictSGW;    
    timer.startTimer();
    computeDictionary(DT_SGW4, es, dictSGW);
    timer.stopTimer("SGW dict time: ");
}

void ShapeEditor::fillHoles(bool skipExternalBoundary)
{
    ZGeom::identifyMeshBoundaries(*mMesh);
    int numBoundaries = ZGeom::getMeshBoundaryLoops(*mMesh).size();
    if (skipExternalBoundary) numBoundaries--;
    if (numBoundaries <= 0) {
        std::cout << "No hole detected!" << std::endl;
        return;
    }

    auto boundaryLoopEdges = ZGeom::getMeshBoundaryLoopHalfEdges(*mMesh);
    for (int i = 0; i < numBoundaries; ++i) {
        vector<int>& boundaryEdges = boundaryLoopEdges[i];
        fillBoundedHole(boundaryEdges);
    }

    mMesh->clearNonEssentialAttributes();
    ZGeom::gatherMeshStatistics(*mMesh);
    this->resetStoredCoordinates();
    visualizeBoundaries();
}

void ShapeEditor::fillBoundedHole(const std::vector<int>& boundaryEdgeIdx)
{
    using namespace std;
    using namespace ZGeom;
    int nOldVerts = mMesh->vertCount();
    int nOldEdges = mMesh->halfEdgeCount();
    int nOldFaces = mMesh->faceCount();

    int N = (int)boundaryEdgeIdx.size();    // number of boundary vertices
    vector<int> boundaryVertIdx(N);
    vector<CVertex*> boundaryVertPtr(N);
    vector<Vec3d> boundaryVertPos(N);
    for (int i = 0; i < N; ++i) {
        boundaryVertIdx[i] = mMesh->getHalfEdge(boundaryEdgeIdx[i])->getVertIndex(0);
        boundaryVertPtr[i] = mMesh->vert(boundaryVertIdx[i]);
        boundaryVertPos[i] = boundaryVertPtr[i]->pos();
    }

    //////////////////////////////////////////////////////////////////////////

    /* triangulating hole */
    vector<vector<WeightSet>> weights(N);
    vector<vector<int>> midIdx(N);
    for (int i = 0; i < N; ++i) {
        weights[i].resize(N);
        midIdx[i].resize(N, 0);
    }
    for (int i = 0; i < N - 2; ++i) {
        int vIdx1 = boundaryVertIdx[i], vIdx2 = boundaryVertIdx[i + 1], vIdx3 = boundaryVertIdx[i + 2];
        CVertex *v1 = mMesh->vert(vIdx1), *v2 = mMesh->vert(vIdx2), *v3 = mMesh->vert(vIdx3);
        CHalfEdge *e12 = v1->adjacentTo(v2), *e23 = v2->adjacentTo(v3);
        Vec3d vn = triNormal(v1->pos(), v3->pos(), v2->pos());
        Vec3d vn12 = e12->getAttachedFace()->calNormal(), vn23 = e23->getAttachedFace()->calNormal();
        weights[i][i + 2].angle = max(vecAngle(vn, vn12), vecAngle(vn, vn23));
        weights[i][i + 2].area = triArea(boundaryVertPos[i], boundaryVertPos[i + 1], boundaryVertPos[i + 2]);
        midIdx[i][i + 2] = i + 1;
    }

    for (int j = 3; j < N; ++j) {
        for (int i = 0; i < N - j; ++i) {
            int k = i + j;
            CVertex *vi = boundaryVertPtr[i], *vk = boundaryVertPtr[k];
            int minMidIdx = -1;
            WeightSet minWeight_ik(4, 1e30);

            for (int m = i + 1; m < k; ++m) {
                CVertex *vm = boundaryVertPtr[m];
                WeightSet triWeight;
                triWeight.area = ZGeom::triArea(boundaryVertPos[i], boundaryVertPos[m], boundaryVertPos[k]);

                Vec3d vn_ikm = triNormal(boundaryVertPos[i], boundaryVertPos[k], boundaryVertPos[m]);
                Vec3d vn_im, vn_mk;
                if (m == i + 1)
                    vn_im = vi->adjacentTo(vm)->getAttachedFace()->calNormal();
                else {
                    int l1 = midIdx[i][m];
                    vn_im = triNormal(boundaryVertPos[i], boundaryVertPos[m], boundaryVertPos[l1]);
                }
                if (k == m + 1)
                    vn_mk = vm->adjacentTo(vk)->getAttachedFace()->calNormal();
                else {
                    int l2 = midIdx[m][k];
                    vn_mk = triNormal(boundaryVertPos[m], boundaryVertPos[k], boundaryVertPos[l2]);
                }
                triWeight.angle = max(vecAngle(vn_ikm, vn_im), vecAngle(vn_ikm, vn_mk));
                if (k == i + N - 1) {
                    Vec3d vn_ik = vk->adjacentTo(vi)->getAttachedFace()->calNormal();
                    triWeight.angle = max(triWeight.angle, vecAngle(vn_ikm, vn_ik));
                }

                WeightSet weight_imk = weights[i][m] + weights[m][k] + triWeight;
                if (weight_imk < minWeight_ik) {
                    minWeight_ik = weight_imk;
                    minMidIdx = m;
                }
            }
            midIdx[i][k] = minMidIdx;
            weights[i][k] = minWeight_ik;
        }
    }

    MeshLineList triangulationLines;
    vector<vector<int>> patchTri;
    queue<pair<int, int>> tracePairs;
    tracePairs.push(make_pair(0, N - 1));
    while (!tracePairs.empty()) {
        pair<int, int> lineIdx = tracePairs.front();
        tracePairs.pop();
        int i = lineIdx.first, k = lineIdx.second;
        if (i + 2 == k) {
            patchTri.push_back(vector < int > { i, i + 1, k });
            triangulationLines.push_back(LineSegment(boundaryVertPos[i], boundaryVertPos[k]));
        }
        else {
            int o = midIdx[i][k];
            if (o > i + 1) tracePairs.push(make_pair(i, o));
            patchTri.push_back(vector < int > { i, o, k });
            triangulationLines.push_back(LineSegment(boundaryVertPos[i], boundaryVertPos[k]));
            if (o < k - 1) tracePairs.push(make_pair(o, k));
        }
    }

    cout << "boundary vert count: " << boundaryVertIdx.size() << endl;
    cout << "patching face count: " << patchTri.size() << endl;

    vector<CHalfEdge*> patchEdges;
    vector<CFace*> patchFaces;
    for (vector<int> tri : patchTri) {
        CVertex *vi = mMesh->vert(boundaryVertIdx[tri[0]]),
                *vj = mMesh->vert(boundaryVertIdx[tri[1]]),
                *vk = mMesh->vert(boundaryVertIdx[tri[2]]);
        // construct new face (vi->vk->vj->vi)
        CFace *f = new CFace(3);
        CHalfEdge *eik = new CHalfEdge(), *ekj = new CHalfEdge(), *eji = new CHalfEdge();   // be careful with the clockwise
        eik->setVertOrigin(vi);
        ekj->setVertOrigin(vk);
        eji->setVertOrigin(vj);         vi->addHalfEdge(eik); vj->addHalfEdge(eji); vk->addHalfEdge(ekj);
        CMesh::makeFace(eik, ekj, eji, f);
        patchEdges.push_back(eik); patchEdges.push_back(ekj); patchEdges.push_back(eji);
        patchFaces.push_back(f);
    }
    for (CHalfEdge* he : patchEdges) {
        if (he->twinHalfEdge() == NULL) {
            CVertex *v1 = he->vert(0), *v2 = he->vert(1);
            CHalfEdge* e21 = v2->adjacentTo(v1);
            CHalfEdge::makeTwins(he, e21);
        }
    }
    for (size_t i = 0; i < patchEdges.size(); ++i) {
        patchEdges[i]->m_eIndex = nOldEdges + i;
        mMesh->m_vHalfEdges.push_back(patchEdges[i]);
    }
    for (size_t i = 0; i < patchFaces.size(); ++i) {
        patchFaces[i]->m_fIndex = nOldFaces + i;
        mMesh->m_vFaces.push_back(patchFaces[i]);
    }
    /*  end of triangulating hole  */

    //////////////////////////////////////////////////////////////////////////

    /* Delaunay-like refinement */
    set<CHalfEdge*> boundaryHalfEdges;
    for (int ei : boundaryEdgeIdx) boundaryHalfEdges.insert(mMesh->getHalfEdge(ei));

    map<int, double> lengthAttr;
    for (int i = 0; i < N; ++i) {
        CVertex* pv = boundaryVertPtr[i];
        double avgLen(0);
        for (CHalfEdge *he : pv->m_HalfEdges) avgLen += he->length();
        avgLen /= pv->outValence();
        lengthAttr[pv->getIndex()] = avgLen;
    }

    double alpha = 3.; //sqrt(2.);
    bool noFaceSplit(true), noEdgeSwap(true);
    while (true)
    {
        vector<CFace*> splitCandidates;
        for (int i = nOldFaces; i < mMesh->faceCount(); ++i) {
            CFace* face = mMesh->getFace(i);
            Vec3d vc = mMesh->getFace(i)->calBarycenter();
            CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
            CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
            double vcLenAttr = (lengthAttr[fv[0]->getIndex()] + lengthAttr[fv[1]->getIndex()] + lengthAttr[fv[2]->getIndex()]) / 3.0;
            bool splitTest = true;
            for (int m = 0; m < 3; ++m) {
                double a = alpha * (fv[m]->pos() - vc).length();
                if (a <= max(lengthAttr[fv[m]->getIndex()], vcLenAttr)) {
                    splitTest = false; break;
                }
            }
            // do faceSplit if splitTest passed
            if (splitTest) splitCandidates.push_back(face);
        }

        if (splitCandidates.empty()) break;
        else {
            for (CFace* face : splitCandidates) {
                CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                double vcLenAttr = (lengthAttr[fv[0]->getIndex()] + lengthAttr[fv[1]->getIndex()] + lengthAttr[fv[2]->getIndex()]) / 3.0;
                CVertex* centroid = mMesh->faceSplit3(face->getFaceIndex());
                lengthAttr[centroid->getIndex()] = vcLenAttr;
            }
        }

        // relaxing all interior half-edges
        int swapCount = 0;
        while (true && swapCount++ < 100) {
            noEdgeSwap = true;
            for (int i = nOldEdges; i < mMesh->m_vHalfEdges.size(); ++i) {
                if (boundaryHalfEdges.find(mMesh->getHalfEdge(i)->twinHalfEdge()) != boundaryHalfEdges.end()) continue;
                if (mMesh->relaxEdge(mMesh->getHalfEdge(i))) {
                    noEdgeSwap = false; // break;
                }
            }
            if (noEdgeSwap) break;
        }
    }
    /*  end of Delaunay refinement  */
   
    MeshRegion newHole;
    newHole.vert_on_boundary = boundaryVertIdx;
    for (int i = nOldVerts; i < mMesh->vertCount(); ++i)
        newHole.vert_inside.push_back(i);
    this->filled_boundaries.push_back(newHole);
}

void ShapeEditor::visualizeBoundaries()
{
    MeshLineList boundaryLines;
    for (MeshRegion &fhv : filled_boundaries) {
        int nEdges = fhv.he_on_boundary.size();
        for (int i = 0; i < nEdges; ++i) {
            CHalfEdge* he = mMesh->getHalfEdge(fhv.he_on_boundary[i]);
            LineSegment ls(he->vert(0)->pos(), he->vert(1)->pos(), false);
            ls.color1 = ZGeom::ColorAzure;
            boundaryLines.push_back(ls);
        }
    }
    boundaryLines.lineWidthScale = 2.0;
    mMesh->addAttrLines(boundaryLines, "hole_boundaries");
    emit meshLineFeatureChanged();
}

void ShapeEditor::holeFairing()
{
    holeFairingFourierOMP();
}

void ShapeEditor::holeFairingFourierOMP()
{
    const int totalVertCount = mMesh->vertCount();
    const MeshCoordinates& coordOld = getOldMeshCoord();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();

    std::cout << "==== Do Eigendecomposition ====\n";
    int eigenCount = 500;					// -1 means full decomposition
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
    Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    vector<int> vMask(totalVertCount, 1);
    for (MeshRegion &fhv : filled_boundaries) {
        for (int vi : fhv.vert_inside) vMask[vi] = 0;
    }

    vector<VecNd> vCoordRecovered(3);
    vector<SparseCoding> vCodingInpaint(3);
    for (int i = 0; i < 3; ++i) {
        vCoordRecovered[i] = singleChannelSparseInpaint(vOriginalCoords[i], vMask, dictMHB, vCodingInpaint[i]);
        for (int j = 0; j < totalVertCount; ++j) {
            if (vMask[j]) vCoordRecovered[i][j] = vOriginalCoords[i][j];
        }
    }
    MeshCoordinates coordRecovered(totalVertCount, vCoordRecovered);
    addCoordinate(coordRecovered, "faired_OMP");
    std::cout << "OMP fairing finished!" << std::endl;
    changeCoordinates("faired_OMP");
}

ZGeom::DenseMatrixd inpaintLARS(const DenseMatrixd& matCoord, const DenseMatrixd& matDict, const std::vector<int>& vMissingIdx, double eps)
{
    g_engineWrapper.addDenseMat(matCoord, "coord");
    g_engineWrapper.addDenseMat(matDict, "dict");

    ZGeom::VecNd vecMissing(vMissingIdx.size());
    for (int i = 0; i < vecMissing.size(); ++i) vecMissing[i] = double(vMissingIdx[i] + 1);
    g_engineWrapper.addColVec(vecMissing, "missing_idx");
    g_engineWrapper.addDoubleScalar(eps, "lambda");

    CStopWatch timer;
    timer.startTimer();
    g_engineWrapper.eval("[coord_est,~,err] =  zmesh_inpaint_lars(dict, coord, missing_idx, lambda);");
    timer.stopTimer("LARS time: ");
    ZGeom::DenseMatrixd result = g_engineWrapper.getDenseMat("coord_est");

    return result;
}

ZGeom::DenseMatrixd inpaintL1LS(const DenseMatrixd& matCoord, const DenseMatrixd& matDict, const std::vector<int>& vMissingIdx, double lambda, double tol)
{
    g_engineWrapper.addDenseMat(matCoord, "coord");
    g_engineWrapper.addDenseMat(matDict, "dict");

    ZGeom::VecNd vecMissing(vMissingIdx.size());
    for (int i = 0; i < vecMissing.size(); ++i) vecMissing[i] = double(vMissingIdx[i] + 1);
    g_engineWrapper.addColVec(vecMissing, "missing_idx");
    g_engineWrapper.addDoubleScalar(lambda, "lambda");
    g_engineWrapper.addDoubleScalar(tol, "tol");

    g_engineWrapper.eval("[coord_est,~,err] =  zmesh_inpaint_l1ls(dict, coord, missing_idx, lambda, tol);");
    ZGeom::DenseMatrixd result = g_engineWrapper.getDenseMat("coord_est");

    return result;
}

void ShapeEditor::holeFairingFourierLARS()
{
    const int totalVertCount = mMesh->vertCount();
    const MeshCoordinates& coordOld = getOldMeshCoord();

    std::cout << "==== Do Eigendecomposition ====\n";
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    int eigenCount = -1;    // -1 means full decomposition
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
    Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    vector<int> inpaintVert;
    for (MeshRegion &fhv : filled_boundaries) {
        for (int vi : fhv.vert_inside) inpaintVert.push_back(vi);
    }
    ZGeom::DenseMatrixd matCoordOld = coordOld.toDenseMatrix();
    ZGeom::DenseMatrixd matDict = dictMHB.toDenseMatrix();
    double eps = 0.0001;
    
    ZGeom::DenseMatrixd matCoordInpainted = inpaintLARS(matCoordOld, matDict, inpaintVert, eps);

    MeshCoordinates coordInpainted(totalVertCount);
    coordInpainted.fromDenseMatrix(matCoordInpainted);

    addCoordinate(coordInpainted, "faired_LARS");
    std::cout << "LARS fairing finished!" << std::endl;
    changeCoordinates("faired_LARS");
}

void ShapeEditor::holeFairingFourierLS()
{
}

void ShapeEditor::holeFairingLeastSquare()
{
    using namespace std;
    revertCoordinates();
    const int totalVertCount = mMesh->vertCount();
    const MeshCoordinates& coordOld = getOldMeshCoord();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();
    MeshLaplacian matLaplacian;
    matLaplacian.constructTutte(mMesh);
    set<int> controlVerts;
    for (int i = 0; i < totalVertCount; ++i) controlVerts.insert(i);
    for (MeshRegion &fhv : filled_boundaries) {
        for (int vi : fhv.vert_inside) controlVerts.erase(vi);
    }
    int controlVertCount = controlVerts.size();
    vector<int> vControlVerts(controlVerts.begin(), controlVerts.end());

    double controlWeight = 1.0;
    vector<int> rowIdxA, colIdxA;
    vector<double> valsA;
    matLaplacian.getSparseMatrix().convertToCOO(rowIdxA, colIdxA, valsA, ZGeom::MAT_FULL);
    for (int i = 0; i < controlVertCount; ++i) {
        rowIdxA.push_back(totalVertCount + i + 1);
        colIdxA.push_back(vControlVerts[i] + 1);
        valsA.push_back(controlWeight);
    }
    SparseMatrix<double> matA(totalVertCount + controlVertCount, totalVertCount);
    matA.convertFromCOO(totalVertCount + controlVertCount, totalVertCount, rowIdxA, colIdxA, valsA);

    DenseMatrixd matB(totalVertCount + controlVertCount, 3);
    for (int i = 0; i < controlVertCount; ++i) {
        for (int j = 0; j < 3; ++j)
            matB(totalVertCount + i, j) = controlWeight * vOriginalCoords[j][vControlVerts[i]];
    }

    g_engineWrapper.addSparseMat(matA, "matA");
    g_engineWrapper.addDenseMat(matB, "matB");
    g_engineWrapper.eval("matX=matA\\matB;");
    DenseMatrixd matX = g_engineWrapper.getDenseMat("matX");

    vector<VecNd> vCoordRecovered = vOriginalCoords;
    for (int i = 0; i < totalVertCount; ++i) {
        if (controlVerts.find(i) == controlVerts.end()) {
            // replace vertices inside holes with recovered coordinate
            for (int j = 0; j < 3; ++j)  vCoordRecovered[j][i] = matX(i, j);     
        }
    }    
    MeshCoordinates coordRecovered(totalVertCount, vCoordRecovered);
    addCoordinate(coordRecovered, "faired_least_square");
    std::cout << "Least Square Fairing finished!\n";
    changeCoordinates("faired_least_square");
}

void ShapeEditor::holeEstimateCurvature()
{
    revertCoordinates();
    int totalVertCount = mMesh->vertCount();
    vector<int> affectedVert;
    for (MeshRegion &fhv : filled_boundaries) {
        for (int vi : fhv.vert_inside) affectedVert.push_back(vi);
        for (int vi : fhv.vert_on_boundary) affectedVert.push_back(vi);
    }
    vector<int> vMask(totalVertCount, 1);
    for (int i : affectedVert) vMask[i] = 0;

    int eigenCount = 300;					// -1 means full decomposition
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    std::cout << "==== Do Eigendecomposition ====\n";
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
    std::cout << "\n==== Compute Dictionaries ====\n";
    Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    MeshLineList vNormalLines;

    /* normals calculated from filled model */
    ZGeom::calMeshAttrVertNormals(*mMesh);
    const vector<Vec3d>& vNormals1 = ZGeom::getMeshVertNormals(*mMesh);
    for (int vi : affectedVert) {
        LineSegment ls(mMesh->vertPos(vi), (Vec3d)vNormals1[vi], true);
        ls.color1 = ZGeom::ColorBlue; ls.color2 = ZGeom::ColorBlue;
        vNormalLines.push_back(ls);
    }

    /* inpaint affected normals */
    vector<VecNd> vNormalOld(3), vNormalRecovered(3);
    for (int i = 0; i < 3; ++i) {
        vNormalOld[i].resize(totalVertCount);
        for (int j = 0; j < totalVertCount; ++j)
            vNormalOld[i][j] = vNormals1[j][i];
    }

    vector<SparseCoding> vCodingInpaint(3);
    for (int i = 0; i < 3; ++i) {
        vNormalRecovered[i] = singleChannelSparseInpaint(vNormalOld[i], vMask, dictMHB, vCodingInpaint[i]);
    }
    vector<ZGeom::Vec3d> vNormals2(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) {
        if (vMask[i]) vNormals2[i] = vNormals1[i];
        else {
            vNormals2[i] = ZGeom::Vec3d(vNormalRecovered[0][i], vNormalRecovered[1][i], vNormalRecovered[2][i]);
            vNormals2[i].normalize();
        }
    }
    for (int vi : affectedVert) {
        LineSegment ls(mMesh->vertPos(vi), (Vec3d)vNormals2[vi], true);
        ls.color1 = ZGeom::ColorRed; ls.color2 = ZGeom::ColorRed;
        vNormalLines.push_back(ls);
    }

    /* inpaint mean curvature values */
    ZGeom::calMeshAttrMeanGaussCurvatures(*mMesh);
    VecNd oldCurvature = VecNd(ZGeom::getMeshMeanCurvatures(*mMesh));
    VecNd newCurvature = singleChannelSparseInpaint(oldCurvature, vMask, dictMHB, vCodingInpaint[0]);
    for (int vi = 0; vi < totalVertCount; ++vi) {
        if (vMask[vi]) newCurvature[vi] = oldCurvature[vi];
    }
    auto minmax_curv = std::minmax_element(oldCurvature.c_ptr(), oldCurvature.c_ptr_end());
    double min_curv = *minmax_curv.first, max_curv = *minmax_curv.second;
    vector<Colorf> colorOldCurv(totalVertCount), colorNewCurv(totalVertCount);
    for (int i = 0; i < totalVertCount; ++i) {
        colorOldCurv[i].falseColor((oldCurvature[i] - min_curv) / (max_curv - min_curv));
        colorNewCurv[i].falseColor((newCurvature[i] - min_curv) / (max_curv - min_curv));
    }

    mMesh->addAttrLines(vNormalLines, "hole_affected_vert_normals");
    emit meshLineFeatureChanged();
    addColorSignature("color_old_curvature", colorOldCurv);
    addColorSignature("color_estimated_curvature", colorNewCurv);
    std::cout << "Finished!\n";
}

void ShapeEditor::holeEstimateNormals()
{

}

void ShapeEditor::generateNoise(const std::vector<int>& selectedVerts, double sigma /*= 0.02*/)
{
    std::mt19937 engine(0);
    std::normal_distribution<double> distCoeff(0, sigma);
    MeshCoordinates noisyCoord = getOldMeshCoord();
    for (int vIdx : selectedVerts) {
        for (int c = 0; c < 3; ++c) {
            noisyCoord(vIdx, c) += distCoeff(engine);
        }
    }
    
    addCoordinate(noisyCoord, "coord_hole_noisy");
    changeCoordinates("coord_hole_noisy");
}

void ShapeEditor::inpaintHolesLARS(const std::vector<int>& inpaintVert, double eps)
{
    const int totalVertCount = mMesh->vertCount();
    const MeshCoordinates& coordOld = getOldMeshCoord();

    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    int eigenCount = -1;    // -1 means full decomposition
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
    Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    ZGeom::DenseMatrixd matCoordOld = coordOld.toDenseMatrix();
    ZGeom::DenseMatrixd matDict = dictMHB.toDenseMatrix();

    ZGeom::DenseMatrixd matCoordInpainted = inpaintLARS(matCoordOld, matDict, inpaintVert, eps);

    MeshCoordinates coordInpainted(totalVertCount);
    coordInpainted.fromDenseMatrix(matCoordInpainted);

    addCoordinate(coordInpainted, "hole_inpainted_LARS");
    std::cout << "LARS hole inpainting finished!" << std::endl;
    std::cout << "LARS inpainting error: " << ZGeom::compareCoordRMSE(coordOld, coordInpainted, inpaintVert) << std::endl;

    changeCoordinates("hole_inpainted_LARS");
}

void ShapeEditor::inpaintHolesL1LS(const std::vector<int>& selctedVerts, double lambda /*= 1e-3*/, double tol /*= 1e-3*/)
{
    const int totalVertCount = mMesh->vertCount();
    MeshCoordinates coordOld = getOldMeshCoord();

    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);
    int eigenCount = -1;    // -1 means full decomposition
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
    Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    ZGeom::DenseMatrixd matCoordOld = coordOld.toDenseMatrix();
    ZGeom::DenseMatrixd matDict = dictMHB.toDenseMatrix();

    ZGeom::DenseMatrixd matCoordInpainted = inpaintL1LS(matCoordOld, matDict, selctedVerts, lambda, tol);

    MeshCoordinates coordInpainted(totalVertCount);
    coordInpainted.fromDenseMatrix(matCoordInpainted);

    addCoordinate(coordInpainted, "hole_inpainted_L1LS");
    std::cout << "L1LS hole inpainting finished!" << std::endl;
    std::cout << "L1LS inpainting error: " << ZGeom::compareCoordRMSE(coordOld, coordInpainted, selctedVerts) << std::endl;

    changeCoordinates("hole_inpainted_L1LS");
}

void ShapeEditor::testSurfaceInpainting()
{
    using namespace std;

    ofstream ofs(string("output/inpainting_output_") + mMesh->getMeshName() + ".txt");
    CStopWatch timer;

    int N = mMesh->vertCount();
    ofs << "Mesh size: " << N << endl;
    
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(mMesh);    
    int eigenCount = -1;    // -1 means full decomposition
    if (N > 1000) eigenCount = 1000;
    timer.startTimer();
    const ZGeom::EigenSystem& es = mMeshHelper->prepareEigenSystem(graphLaplacian, eigenCount);
    timer.stopTimer("Decomposition time: ");
    ofs << "Decomposition time: " << timer.getElapsedTime() << endl;
    Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    MeshCoordinates coordOld = getOldMeshCoord();
    ZGeom::DenseMatrixd matCoordOld = coordOld.toDenseMatrix();
    ZGeom::DenseMatrixd matDict = dictMHB.toDenseMatrix();

    default_random_engine generator((unsigned int)time(NULL));
    std::random_device rd;
    std::mt19937 g(rd());
    vector<double> missing_ratios{ 0.05, 0.2, 0.5 };
    
    int count(0);    
    for (double missing_ratio : missing_ratios) {
    // for each missing_ratio
        int nHoleVerts = std::round(N * missing_ratio);
        ofs << "missing ratio: " << missing_ratio << endl;
        vector<int> seed_counts{ nHoleVerts, 5, 1 };
        for (int nSeed : seed_counts) {
        // for different number of hole seeds
            ofs << "Hole seed count: " << nSeed << endl;
            
            // do three tests and use best
            vector<pair<double, double>> testResults;
            for (int k = 0; k < 1; ++k) {
                std::vector<int> seedVerts(N);
                for (int i = 0; i < N; ++i) seedVerts[i] = i;
                std::shuffle(seedVerts.begin(), seedVerts.end(), g);
                seedVerts = std::vector < int > {seedVerts.begin(), seedVerts.begin() + nSeed};
                MeshRegion hole = ZGeom::generateRandomMeshRegion(*mMesh, seedVerts, nHoleVerts);

                double lambda = 1e-2;
                double tol = 1e-3;


                timer.startTimer();
                ZGeom::DenseMatrixd matCoordInpainted = inpaintL1LS(matCoordOld, matDict, hole.vert_inside, lambda, tol);
                timer.stopTimer();
                double inpaintTime = timer.getElapsedTime();

                MeshCoordinates coordInpainted(N);
                coordInpainted.fromDenseMatrix(matCoordInpainted);
                double mse = ZGeom::compareCoordRMSE(coordOld, coordInpainted, hole.vert_inside);

                testResults.push_back(make_pair(mse, inpaintTime));

                if (++count % 5 == 0) cout << "test " << count << " finished!" << endl;
            }

            sort(testResults.begin(), testResults.end(), 
                [](const pair<double, double>& t1, const pair<double, double>& t2) { return t1.first < t2.first; });
            ofs << "Error: " << testResults.front().first << "; Time: " << testResults.front().second << endl;
            ofs << endl;
        } 
        ofs << endl;
    }
    cout << "Finished!" << endl;
}
