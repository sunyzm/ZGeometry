#include "ShapeEditor.h"
#include <iostream>
#include <random>
#include <functional>
#include <iomanip>
#include <ppl.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/util.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/Approximation.h>
#include <ZGeom/SparseSolver.h>
#include <ZGeom/MCA.h>
#include "SpectralGeometry.h"
#include "global.h"

using ZGeom::Dictionary;
using ZGeom::SparseApproxMethod;
using ZGeom::VecNd;
using std::vector;

Vector3D toVector3D(const ZGeom::Vec3d& v) { return Vector3D(v[0], v[1], v[2]); }
ZGeom::Vec3d toVec3d(const Vector3D& v) { return ZGeom::Vec3d(v.x, v.y, v.z); }

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

void ShapeEditor::init( DifferentialMeshProcessor* processor )
{
	mProcessor = processor; 
	mMesh = processor->getMesh();
	mTotalScales = 0;
	mCurCoordID = 0;
	processor->getMesh()->retrieveVertCoordinates(mOriginalCoord); 
	mStoredCoordinates.resize(1);
	mStoredCoordinates[0] = mOriginalCoord;

	std::cout << "Shape editor is initialized!" << std::endl;
}

void ShapeEditor::revertCoordinates()
{
	mMesh->setVertCoordinates(getOldMeshCoord());
}

void ShapeEditor::nextCoordinates()
{
	changeCoordinates((mCurCoordID + 1) % mStoredCoordinates.size());
}

void ShapeEditor::changeCoordinates( int coordID )
{
	if (coordID < 0 || coordID >= mStoredCoordinates.size()) return;
	mCurCoordID = coordID;
	qout.output("Current coordinate ID: " + Int2String(mCurCoordID), OUT_STATUS);

	if (mStoredCoordinates[mCurCoordID].empty())
		qout.output("!! Selected coordinate is empty!", OUT_CONSOLE);
	else 
		mMesh->setVertCoordinates(mStoredCoordinates[mCurCoordID]);
}

void ShapeEditor::setStoredCoordinates( const MeshCoordinates& newCoord, int idx )
{
	if (idx >= mStoredCoordinates.size()) mStoredCoordinates.resize(idx + 1);
	mStoredCoordinates[idx] = newCoord;
}

MeshCoordinates& ShapeEditor::getStoredCoordinate( int idx )
{
	if (idx >= mStoredCoordinates.size()) mStoredCoordinates.resize(idx + 1);
	return mStoredCoordinates[idx];
}

void ShapeEditor::prepareAnchors( int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos ) const
{
	const std::map<int, Vector3D>& anchors = mProcessor->getHandles();
	anchorCount = anchors.size();
	anchorIndex.clear();
	anchorPos.clear();
	for (auto a : anchors) {
		anchorIndex.push_back(a.first);
		anchorPos.push_back(a.second);
	}
}

void ShapeEditor::addNoise( double phi )
{
	assert(phi > 0 && phi < 1);
	const int vertCount = mMesh->vertCount();
	const double avgLen = mMesh->getAvgEdgeLength();
	MeshCoordinates newCoord;
	mMesh->retrieveVertCoordinates(newCoord);

	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0, phi);

	for (int vIdx = 0; vIdx < vertCount; ++vIdx) {
		for (int a = 0; a < 3; ++a) {
			double noise = avgLen * distribution(generator);
			newCoord.getCoordFunc(a)[vIdx] += noise;
		}
	}

	mMesh->setVertCoordinates(newCoord);
}

void ShapeEditor::addColorSignature(const std::string& colorSigName, const std::vector<ZGeom::Colorf>& vColors)
{
	std::vector<ZGeom::Colorf>& vNewColors = mMesh->addColorAttr(colorSigName).attrValue();
	vNewColors = vColors;
	emit signatureComputed(QString(colorSigName.c_str()));
}

void ShapeEditor::continuousReconstruct( int selected, int contCoordIdx )
{
	if (selected < 0 || selected >= 4 || contCoordIdx < 0 || 
		contCoordIdx >= mContReconstructCoords[selected].size()) return;
	if (mContReconstructCoords[selected].empty()) return;

	mMesh->setVertCoordinates(mContReconstructCoords[selected][contCoordIdx]);

	/* visualize the coordinate difference against the original coordinate */
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrAttrColorPosDiff).attrValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrAttrColorPosDiff.c_str()));
	emit coordinateSelected(selected, contCoordIdx);

	/* show SGW approximation features */
#if 0
	if (selected == 2) {
		const int vertCount = mMesh->vertCount();
		MeshProperty* curAtomLocation = mProcessor->retrievePropertyByID(FEATURE_SGW_SOMP);
		if (NULL == curAtomLocation) return;
		MeshFeatureList *mfl = dynamic_cast<MeshFeatureList*>(curAtomLocation);
		mfl->clear();
		int atomIdx = mApproxCoeff[selected][atomCount].index();
		mfl->addFeature(new MeshFeature(atomIdx % vertCount, atomIdx / vertCount));
	}
#endif
}

void ShapeEditor::fourierReconstruct( int nEig )
{
	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord;
	mMesh->retrieveVertCoordinates(oldCoord);
	MeshCoordinates newCoord(vertCount);
	LaplacianType lapType = SymCot;
		
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ZGeom::EigenSystem &mhb = mProcessor->getMHB(lapType);
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
	std::vector<Vector3D> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int vertCount = mMesh->vertCount();
	int hIdx = anchorIndex[0];		// only use the first handle to deform
	Vector3D handleTrans = anchorPos[0] - mMesh->getVertex(hIdx)->getPosition();

	std::vector<int> vFreeIdx;
	mMesh->vertRingNeighborVerts(hIdx, 5, vFreeIdx, true);
//	vFreeIdx.resize(vertCount);
//	for (int i = 0; i < vertCount; ++i) vFreeIdx[i] = i;

	const int freeVertCount = vFreeIdx.size();
	std::vector<double> vDist2Handle(freeVertCount);
	concurrency::parallel_for (0, freeVertCount, [&](int i) {
		double dist = mMesh->calGeodesic(hIdx, vFreeIdx[i]);
		vDist2Handle[i] = dist;
	});
	double distMax = *std::max_element(vDist2Handle.begin(), vDist2Handle.end());

	std::vector<Vector3D> vDeformedPos(freeVertCount);
	for (int i = 0; i < freeVertCount; ++i) {
		CVertex* pv = mMesh->getVertex(vFreeIdx[i]);
		vDeformedPos[i] = pv->getPosition() + handleTrans * (1.0 - vDist2Handle[i]/distMax);
	}

	MeshCoordinates newCoord;
	mMesh->retrieveVertCoordinates(newCoord);
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
	std::vector<Vector3D> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();

	MeshCoordinates oldCoord;
	mMesh->retrieveVertCoordinates(oldCoord);

	const MeshLaplacian& ml = mProcessor->getMeshLaplacian(SymCot);
	const ZGeom::SparseMatrix<double>& matLs = ml.getLS();
	ZGeom::SparseMatVecMultiplier mulLs(matLs, true);	
	ZGeom::VecNd diffCoord[3];
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(vertCount);
		mulLs.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
	}

	std::vector<ZGeom::VecNd> solveRHS(3);
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
	std::vector<Vector3D> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();

	const MeshCoordinates& oldMeshCoord = this->getOldMeshCoord();

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(SymCot).getLS();

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
	std::vector<Vector3D> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();

	MeshCoordinates oldCoord;
	mMesh->retrieveVertCoordinates(oldCoord);

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(SymCot).getLS();
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
	std::vector<Vector3D> anchorPos;
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

	MeshCoordinates oldCoord;
	mMesh->retrieveVertCoordinates(oldCoord);
	g_engineWrapper.addColVec(oldCoord.getXCoord(), "ecx");
	g_engineWrapper.addColVec(oldCoord.getYCoord(), "ecy");
	g_engineWrapper.addColVec(oldCoord.getZCoord(), "ecz");

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(SymCot).getLS();
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
			const Vector3D& oldPos = mMesh->getVertexPosition(anchorIndex[l]);
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

	mProcessor->computeHeatDiffuseMat(tMultiplier);
	ZGeom::SparseSymMatVecSolver& solver = mProcessor->getHeatSolver();
	const MeshLaplacian& laplacian = mProcessor->getMeshLaplacian(CotFormula);
	const ZGeom::SparseMatrix<double>& matW = laplacian.getW();
	const ZGeom::SparseMatrix<double>& matLc = laplacian.getLS();	// negative
	ZGeom::SparseMatVecMultiplier mulW(matW, true);

	MeshCoordinates oldCoord, newCoord;
	mMesh->retrieveVertCoordinates(newCoord);
	oldCoord.resize(vertCount);
	
	for (int n = 0; n < nRepeat; ++n) {
		for (int a = 0; a < 3; ++a) 
			mulW.mul(newCoord.getCoordFunc(a), oldCoord.getCoordFunc(a));				
		for (int a = 0; a < 3; ++a) 
			solver.solve(oldCoord.getCoordFunc(a), newCoord.getCoordFunc(a));
	}

	mMesh->setVertCoordinates(newCoord);
	
	mMesh->scaleAreaToVertexNum();
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
	//testSparseCompression();		
	//sparseDecompositionTest();
	//sparseDecompositionTest2();
	//testArtificialShapeMCA();
	//testArtificailShapeMCA2();
	//testArtificialShapeMCA3();
	//testDictionaryForDecomposition();
	testSparseFeatureFinding();
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
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrAttrColorPartitions).attrValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrAttrColorPartitions.c_str()));
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

//// Test approximation for feature analysis  ////
//
void ShapeEditor::testSparseFeatureFinding()
{
	using ZGeom::combineDictionary;
	revertCoordinates();
	CStopWatch timer;
	std::cout << "\n======== Starting testSparseFeatureFinding ========\n";

	const int totalVertCount = mMesh->vertCount();
	int eigenCount = min(1500, totalVertCount - 1); // -1 means vertCount -1 
	eigenCount = -1;
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	std::vector<ZGeom::VecNd> vOriginalCoords{ oldMeshCoord.getXCoord(), oldMeshCoord.getYCoord(), oldMeshCoord.getZCoord() };
	setStoredCoordinates(oldMeshCoord, 0);

	/* Compute Laplacian and eigendecomposition */
	std::cout << "==== Do Eigendecomposition ====\n";
	MeshLaplacian graphLaplacian, cotLaplacian, aniso1Laplacian, symCotLaplacian;
	ZGeom::EigenSystem esGraph, esAniso, esCot, esSymCot;
	graphLaplacian.constructUmbrella(mMesh);
	graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esGraph);
	aniso1Laplacian.constructAnisotropic2(mMesh);
	//aniso1Laplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esAniso);
	cotLaplacian.constructCotFormula(mMesh);
	//cotLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esCot);
	symCotLaplacian.constructSymCot(mMesh);
	//symCotLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, esSymCot);

	/* Computer Dictionary */
	std::cout << "\n==== Compute Dictionaries ====\n";
	std::vector<Dictionary> dict(10);
	computeDictionary(DT_Fourier, esGraph, dict[0]);
	computeDictionary(DT_SGW4, esGraph, dict[1]);
	combineDictionary(dict[0], dict[1], dict[2]);
	//computeDictionary(DT_SGW4, esAniso, dict[3]);
	//combineDictionary(dict[0], dict[3], dict[4]);
	computeDictionary(DT_FourierSpikes, esGraph, dict[5]);

	/* Compute OMP Approximation */
	std::cout << "\n==== Compute OMP Approximation ====\n";
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mApproxMethod = ZGeom::SA_SOMP;
	approxOpts.mCodingSize = 100;
	std::vector<ZGeom::SparseCoding> vOMPCodings;
	std::vector<VecNd> vApproximatedCoords;
	multiChannelSparseApproximate(vOriginalCoords, dict[1], vOMPCodings, approxOpts);
	multiChannelSparseReconstruct(dict[1], vOMPCodings, vApproximatedCoords);
	MeshCoordinates coordApproxOMP(totalVertCount, vApproximatedCoords);
	std::cout << "reconstruction error: " << oldMeshCoord.difference(coordApproxOMP) << '\n';
	setStoredCoordinates(coordApproxOMP, 1);
	changeCoordinates(1);

	/* add the origin of selected SGW as feature point */
	auto addFeatures = [&](const ZGeom::SparseCoding& sc, int numGlobalAtom, std::string featureName) {
		MeshFeatureList vSparseFeatures;
		for (auto c : sc.getApproxItems()) {
			int scale = (c.index() - numGlobalAtom) / totalVertCount,
				idx = (c.index() - numGlobalAtom) % totalVertCount;
			if (scale > 0) vSparseFeatures.addFeature(new MeshFeature(idx, scale));
		}
		mMesh->addAttrMeshFeatures(vSparseFeatures, featureName);
		std::cout << "Detected " << featureName << " (" << vSparseFeatures.size() << "): ";
		for (auto f : vSparseFeatures.getFeatureVector()) std::cout << f->m_index << ", ";
		std::cout << '\n';
	};
	addFeatures(vOMPCodings[0], 0, "mesh_feature_SOMP");

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
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrAttrColorPartitions).attrValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrAttrColorPartitions.c_str()));
	
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
	using ZGeom::VecNd;
	/************************************************************************/
	/* 1. Do eigendecomposition and compute dictionary MHB and SGW;         */
	/* 2. Compute decomposition of shapes against dictionaries via S-OMP    */
	/*    and MCA;															*/	
	/* 3. Visualize the magnitude of separated signals as well as           */
	/*    decomposed shapes (seemingly not meaningful)                      */
	/* 4. Focus on non-partitioned shape for now                            */
	/************************************************************************/

	revertCoordinates();
	std::cout << "\n======== Starting sparseDecompositionTest2 ========\n";
	const int totalVertCount = mMesh->vertCount();	
	const double originalAvgEdgeLen = mMesh->getAvgEdgeLength();
	const double bbDiag = mMesh->getBoundingBox().length() * 2;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	std::vector<ZGeom::VecNd> vOriginalCoords { oldMeshCoord.getXCoord(), oldMeshCoord.getYCoord(), oldMeshCoord.getZCoord() };	

	/// Do Eigendecomposition
	//
	std::cout << "==== Do Eigendecomposition ====\n";
	mShapeApprox.init(mMesh);
	mShapeApprox.doSegmentation(-1);		// -1 means no segmentation 
	int eigenCount = -1;					// -1 means full decomposition
	mShapeApprox.doEigenDecomposition(Umbrella, eigenCount);	// do full decomposition
	const ZGeom::EigenSystem& es = mShapeApprox.getSubEigenSystem(0);

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

	std::cout << "- #MHB: " << vOMPCodings[0].size() << "\t#SGW: " << vAddedCoeff[0].size() << '\n';
	vOMPCodings[0].sortByCoeff(); vAddedCoeff[0].sortByCoeff();
	std::cout << "-- MHB coefficients: ";
	for (auto p : vOMPCodings[0].getApproxItems()) std::cout << "(" << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";
	std::cout << "-- SGW coefficients: ";
	for (auto p : vAddedCoeff[0].getApproxItems()) std::cout << '(' << p.index() << ", " << p.coeff() << ") ";
	std::cout << "\n";

	// compute and display color signatures
	const VecNd& vX0 = coordApproxOMP.getXCoord(), vX1 = coordAltered.getXCoord();
	auto px1 = vX0.min_max_element(), px2 = vX1.min_max_element();
	//double x_min = min(px1.first, px2.first), x_max = max(px1.second, px2.second);
	double x_min = px1.first, x_max = px1.second;
	std::string colorStr0 = "color_x_coord_base";
	std::vector<ZGeom::Colorf> vColors0(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) vColors0[i].falseColor((vX0[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr0, vColors0);
	std::string colorStr1 = "color_x_coord_altered";
	std::vector<ZGeom::Colorf> vColors1(totalVertCount);
	for (int i = 0; i < totalVertCount; ++i) vColors1[i].falseColor((vX1[i] - x_min) / (x_max - x_min), 1.f, ZGeom::CM_JET);
	addColorSignature(colorStr1, vColors1);
	std::string colorStr2 = "color_x_coord_diff";
	std::vector<ZGeom::Colorf> vColors2(totalVertCount);
	std::vector<double> vDiff = (vX0-vX1).toStdVector();
	for (auto& d : vDiff) d = std::fabs(d);
	double maxDiff = *std::max_element(vDiff.begin(), vDiff.end());
	for (int i = 0; i < totalVertCount; ++i) vColors2[i].falseColor(min(float(vDiff[i]/(bbDiag*0.02)), 1.f), 1.f, ZGeom::CM_JET);
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
		vColorsMCA2[i].falseColor(min(float(fabs(vXT[i]/(0.02*bbDiag))), 1.f), 1.f, ZGeom::CM_JET);
	}
	addColorSignature(colorStrMCA1, vColorsMCA1);
	addColorSignature(colorStrMCA2, vColorsMCA2);

#if 0
	std::vector<ZGeom::SparseCoding> vCodings;
	std::vector<ZGeom::VecNd> vReconstructedCoords;

	Dictionary dict;
	DictionaryType dictType = DT_SGW4MHB;
	computeDictionary(dictType, es, dict);
	int nDictSize = dict.size();
	int nEigenCount = es.eigVecCount();
	int nWaveletLevels = 4;
	
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mCodingSize = 300;
	approxOpts.mApproxMethod = ZGeom::SA_SOMP;
	multiChannelSparseApproximate(vOriginalCoords, dict, vCodings, approxOpts);
	multiChannelSparseReconstruct(dict, vCodings, vReconstructedCoords);
	MeshCoordinates coordReconstructedOMP(totalVertCount, vReconstructedCoords);
	setStoredCoordinates(coordReconstructedOMP, 1);
	std::cout << "S-OMP Reconstruction error: " << oldMeshCoord.difference(coordReconstructedOMP) << '\n';
	changeCoordinates(1);
#if 0
	std::vector<ZGeom::SparseCoding> codingMHB, codingSGW;
	multiChannelSplitSparseCoding(vCodings, codingSGW, codingMHB, 
		[=](const ZGeom::SparseCodingItem& ci) { return ci.index() < nWaveletLevels * totalVertCount; }
	);
	std::cout << "#MHB atoms selected: " << codingMHB[0].size() << '\n';
	std::cout << "#SGW atoms selected: " << codingSGW[0].size() << '\n';
	std::vector<ZGeom::VecNd> vSGWCoords, vMHBCoords;
	multiChannelSparseReconstruct(dict, codingSGW, vSGWCoords);
	multiChannelSparseReconstruct(dict, codingMHB, vMHBCoords);
	MeshCoordinates coordSGW(totalVertCount, vSGWCoords), coordMHB(totalVertCount, vMHBCoords);
	setStoredCoordinates(coordMHB, 2);
	setStoredCoordinates(coordSGW, 3);
#endif
#endif 

#if 0 
	/// SOMP
	int codingSize = 300;
	mShapeApprox.constructDictionaries(DT_SGW4MHB);
	ZGeom::SparseApproxMethod saMethod = ZGeom::SA_SOMP;
	
	mShapeApprox.findSparseRepresentationBySize(saMethod, codingSize);
		
	MeshCoordinates reconstructedMeshCoord;
	mShapeApprox.findSparseRepresentationBySize(saMethod, codingSize);
	mShapeApprox.doSparseReconstructionBySize(-1, reconstructedMeshCoord);

	const ZGeom::Dictionary& dictOMP = mShapeApprox.getSubDictionary(0);
	const std::vector<ZGeom::SparseCodingItem>* vCoeff[3] = {
		&mShapeApprox.mSubMeshApprox[0].getSparseCoding(0),
		&mShapeApprox.mSubMeshApprox[0].getSparseCoding(1),
		&mShapeApprox.mSubMeshApprox[0].getSparseCoding(2)
	};
	int dictSize = dictOMP.size();
	int actualCodingSize = int(vCoeff[0]->size());
	MeshCoordinates coordMHB(totalVertCount), coordSGW(totalVertCount);
	int mhbAtomCount(0), sgwAtomCount(0);
	for (int i = 0; i < actualCodingSize; ++i) {
		for (int c = 0; c < 3; ++c) {
			const ZGeom::SparseCodingItem& sc = (*vCoeff[c])[i];
			if (sc.index() < dictSize - es.eigVecCount()) {
				coordSGW.getCoordFunc(c) += sc.coeff() * dictOMP[sc.index()];
				if (c == 0) sgwAtomCount++;
			}
			else {
				coordMHB.getCoordFunc(c) += sc.coeff() * dictOMP[sc.index()];
				if (c == 0) mhbAtomCount++;
			}
		}
	}
	std::cout << "S-OMP Reconstruction error: " << oldMeshCoord.difference(reconstructedMeshCoord) << '\n';
	setStoredCoordinates(reconstructedMeshCoord, 1);
	setStoredCoordinates(coordMHB, 2);
	setStoredCoordinates(coordSGW, 3);
	std::cout << "#MHB atoms selected: " << mhbAtomCount << '\n';
	std::cout << "#SGW atoms selected: " << sgwAtomCount << '\n';

	{
		std::vector<double> sgwRatio(totalVertCount);
		for (int i = 0; i < totalVertCount; ++i) {
			sgwRatio[i] = coordSGW[i].length() / coordMHB[0].length();
		}
		double maxRatio = *std::max_element(sgwRatio.begin(), sgwRatio.end());
		for (double& d : sgwRatio) d /= maxRatio;
		std::string colorStr = "color_omp_sgw_ratio";
		std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(colorStr).attrValue();
		vColors.resize(totalVertCount);
		for (int i = 0; i < totalVertCount; ++i) vColors[i].falseColor(sgwRatio[i], 1.0f, ZGeom::CM_JET);
		emit signatureComputed(QString(colorStr.c_str()));
	}

	/// MCA
	Dictionary dictMHB, dictSGW4;
	computeDictionary(DT_Fourier, es, dictMHB);
	computeDictionary(DT_SGW4, es, dictSGW4);
	std::vector<const Dictionary*> vDicts{ &dictMHB, &dictSGW4 };
	std::vector<ZGeom::VecNd> vOriCoords{ oldMeshCoord.getCoordFunc(0), oldMeshCoord.getCoordFunc(1), oldMeshCoord.getCoordFunc(2) };
	std::vector<std::vector<ZGeom::SparseCoding> > vCoordCodings(3);
	ZGeom::MCAoptions mcaOpts;
	mcaOpts.nIter = 500;
	mcaOpts.threshMode = ZGeom::MCAoptions::HARD_THRESH;
	mcaOpts.threshStrategy = ZGeom::MCAoptions::MIN_OF_MAX;

	std::vector<ZGeom::VecNd> vMHBCoords(3), vSGWCoords(3);
	for (int c = 0; c < 3; ++c) {
		singleChannelMCA(vOriCoords[c], vDicts, vCoordCodings[c], &mcaOpts);
		vMHBCoords[c] = SparseReconstructSingleChannel(dictMHB, vCoordCodings[c][0]);
		vSGWCoords[c] = SparseReconstructSingleChannel(dictSGW4, vCoordCodings[c][1]);
	}
	
	MeshCoordinates coordMCA(totalVertCount), coordMCA_MHB(totalVertCount), coordMCA_SGW(totalVertCount);
	for (int c = 0; c < 3; ++c) {
		coordMCA.getCoordFunc(c) = vMHBCoords[c] + vSGWCoords[c];
		coordMCA_MHB.getCoordFunc(c) = vMHBCoords[c];
		coordMCA_SGW.getCoordFunc(c) = vSGWCoords[c];
	}
	std::cout << "MCA reconstruction error: " << oldMeshCoord.difference(coordMCA) << '\n';
	setStoredCoordinates(coordMCA, 4);
	setStoredCoordinates(coordMCA_MHB, 5);
	setStoredCoordinates(coordMCA_SGW, 6);

	{
		std::vector<double> sgwRatio(totalVertCount);
		for (int i = 0; i < totalVertCount; ++i) {
			sgwRatio[i] = coordMCA_SGW[i].length() / coordMCA_MHB[0].length();
		}
		double maxRatio = *std::max_element(sgwRatio.begin(), sgwRatio.end());
		for (double& d : sgwRatio) d /= maxRatio;
		std::string colorStr = "color_mca_sgw_ratio";
		std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(colorStr).attrValue();
		vColors.resize(totalVertCount);
		for (int i = 0; i < totalVertCount; ++i) vColors[i].falseColor(sgwRatio[i], 1.0f, ZGeom::CM_JET);
		emit signatureComputed(QString(colorStr.c_str()));
	}
#endif

	std::cout << '\n';
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
	setStoredCoordinates(coordOld, 0);	// set as base mesh coordinates
	changeCoordinates(0);
	std::cout << "- S-OMP Reconstruction error: " << coordOld.difference(coordApproxOMP) << '\n';

	/// Add SGW components
	//
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
	std::vector<ZGeom::VecNd> vOriginalCoords{ coordOld.getXCoord(), coordOld.getYCoord(), coordOld.getZCoord() };

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
	Dictionary dictMHB, dictSGW, dictHK, dictSpikes;
	computeDictionary(DT_UNIT, es, dictSpikes);
	computeDictionary(DT_Fourier, es, dictMHB);
	computeDictionary(DT_SGW1, es, dictSGW);
	calHKDict(es, 5., dictHK);

	/// Compute OMP Approximation
	//
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

	/// Add SGW components
	//
	std::cout << "\n==== Mix with extra SGW components ====\n";
	std::vector<ZGeom::VecNd> vAlteredCoords = vOriginalCoords;
	int nnz2 = 15;
	std::vector<ZGeom::SparseCoding> vAddedCoeff(3);
	std::mt19937 engine(0);
	std::uniform_int_distribution<int> distNZ(0, dictSpikes.size() - 1);
	std::normal_distribution<double> distCoeff(0, bbDiag * 0.01);
	for (int i = 0; i < nnz2; ++i) {
		int selectedNNZ = distNZ(engine);
		for (int c = 0; c < 3; ++c) {
			double coeff = distCoeff(engine);
			vAlteredCoords[c] += coeff * dictSpikes[selectedNNZ];
			vAddedCoeff[c].addItem(selectedNNZ, coeff);
		}
	}
	MeshCoordinates coordAltered(totalVertCount, vAlteredCoords);
	setStoredCoordinates(coordAltered, 1);

	std::cout << "- #MHB: " << vOMPCodings[0].size() << "\t#HK: " << vAddedCoeff[0].size() << '\n';
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
	std::vector<const Dictionary*> vDicts{ &dictMHB, &dictSGW };
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
		vCartoon[c] = singleChannelSparseReconstruct(dictMHB, vMCACodings[c][0]);
		vTexture[c] = singleChannelSparseReconstruct(dictSGW, vMCACodings[c][1]);
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
