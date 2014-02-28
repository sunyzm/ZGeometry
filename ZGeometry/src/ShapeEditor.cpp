#include "ShapeEditor.h"
#include <iostream>
#include <random>
#include <functional>
#include <iomanip>
#include <ppl.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/Approximation.h>
#include <ZUtil/timer.h>
#include <ZUtil/zassert.h>
#include <ZUtil/zutil_io.h>
#include "global.h"

using ZGeom::Dictionary;

Vector3D toVector3D(const ZGeom::Vec3d& v) { return Vector3D(v[0], v[1], v[2]); }
ZGeom::Vec3d toVec3d(const Vector3D& v) { return ZGeom::Vec3d(v.x, v.y, v.z); }

bool readPursuits(const std::string& pursuitFile, ZGeom::FunctionApproximation *vPursuits[]) 
{
	if (!ZUtil::fileExist(pursuitFile)) 
		return false;

	std::ifstream ifs(pursuitFile.c_str());
	int atomCount;
	ifs >> atomCount;
	
	for (int b = 0; b < 3; ++b) {
		vPursuits[b]->resize(atomCount);
		for (int i = 0; i < atomCount; ++i) {
			ZGeom::ApproxItem& item = (*vPursuits[b])[i];
			ifs >> item.res() >> item.index() >>item.coeff();
		}
	}

	return true;
}

void writePursuits(const std::string& pursuitFile, ZGeom::FunctionApproximation *vPursuits[]) {
	std::ofstream ofs(pursuitFile.c_str());
	int atomCount = vPursuits[0]->size();
	ofs << atomCount << std::endl;
	for (int b = 0; b < 3; ++b) {
		for (int i = 0; i < atomCount; ++i) {
			const ZGeom::ApproxItem& item = (*vPursuits[b])[i];
			ofs << item.res() << ' ' << item.index() << ' ' << item.coeff() << std::endl;
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

void ShapeEditor::init( DifferentialMeshProcessor* processor )
{
	mProcessor = processor; 
	mMesh = processor->getMesh();
	mTotalScales = 0;
	mCurCoordID = 0;
	mStoredCoordinates.resize(4);
	processor->getMesh()->getVertCoordinates(mOriginalCoord); 
	std::cout << "Shape editor is initialized!" << std::endl;
}

void ShapeEditor::revertCoordinates()
{
	mMesh->setVertCoordinates(getOldMeshCoord());
	qout.output("Coordinate reverted", OUT_STATUS);
}

void ShapeEditor::nextCoordinates()
{
	changeCoordinates((mCurCoordID + 1) % mStoredCoordinates.size());
}

void ShapeEditor::changeCoordinates( int coordID )
{
	if (coordID < 0 || coordID >= mStoredCoordinates.size()) return;
	mCurCoordID = coordID;
	std::cout << "-- Current coordinate ID: " << mCurCoordID << '\n';

	if (mStoredCoordinates[mCurCoordID].empty())
		std::cout << "!! Selected coordinate is empty!\n";
	else mMesh->setVertCoordinates(mStoredCoordinates[mCurCoordID]);
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
	mMesh->getVertCoordinates(newCoord);

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

void ShapeEditor::continuousReconstruct( int selected, int contCoordIdx )
{
	if (selected < 0 || selected >= 4 || contCoordIdx < 0 || 
		contCoordIdx >= mContReconstructCoords[selected].size()) return;
	if (mContReconstructCoords[selected].empty()) return;

	mMesh->setVertCoordinates(mContReconstructCoords[selected][contCoordIdx]);

	/* visualize the coordinate difference against the original coordinate */
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrColorPosDiff).getValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrColorPosDiff.c_str()));
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
	mMesh->getVertCoordinates(oldCoord);
	MeshCoordinates newCoord(vertCount);
	LaplacianType lapType = SymCot;
		
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ZGeom::EigenSystem &mhb = mProcessor->getMHB(lapType);
	ZUtil::logic_assert(nEig <= mhb.eigVecCount(), "Insufficient eigenvectors in mhb");

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
	mMesh->getVertCoordinates(newCoord);
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
	mMesh->getVertCoordinates(oldCoord);

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(SymCot).getLS();
	ZGeom::SparseMatVecMultiplier mulLs(matLs, true);	
	ZGeom::VecNd diffCoord[3];
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(vertCount);
		mulLs.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
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
	matOptS.copyElements(matLs);
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
	mMesh->getVertCoordinates(oldCoord);

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
	mMesh->getVertCoordinates(oldCoord);
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

void ShapeEditor::deformSpectralWavelet()
{
	CStopWatch timer;	
	timer.startTimer();

	int anchorCount(0);
	std::vector<int> anchorIndex;
	std::vector<Vector3D> anchorPos;
	prepareAnchors(anchorCount, anchorIndex, anchorPos);
	if (anchorCount == 0) {
		//std::cout << "At least one anchor need to be picked!" << std::endl;
//		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();

	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	g_engineWrapper.addColVec(oldCoord.getXCoord(), "ecx");
	g_engineWrapper.addColVec(oldCoord.getYCoord(), "ecy");
	g_engineWrapper.addColVec(oldCoord.getZCoord(), "ecz");
	
	ZGeom::DenseMatrixd& matSGW = mProcessor->getWaveletMat();
	if (matSGW.empty()) mProcessor->computeSGW1(CotFormula);
	const int waveletCount = matSGW.rowCount();
	g_engineWrapper.addDenseMat(matSGW, "matSGW");
	
#if 1
	ZGeom::VecNd diffCoord[3];
	ZGeom::DenseMatVecMultiplier mulMixedW(matSGW);	
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(waveletCount);
		mulMixedW.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
	}
#endif

	std::vector<int> fixedVerts;
	std::set<int> freeVerts;
	for (int i = 0; i < anchorCount; ++i) {
		std::set<int> sFree;
		mMesh->vertRingNeighborVerts(anchorIndex[i], 3, sFree, true);
		for (int v : sFree) freeVerts.insert(v);
	}
#if 0
	for (int i = 0; i < vertCount; ++i) {
		if (freeVerts.find(i) == freeVerts.end())
			fixedVerts.push_back(i);
	}
#endif
	const int fixedCount = fixedVerts.size();

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(waveletCount + anchorCount + fixedCount, 0);

		solveRHS[i].copyElements(diffCoord[i], 0);

		for (int l = 0; l < anchorCount; ++l) {
			solveRHS[i][waveletCount + l] = anchorWeight * anchorPos[l][i];
		}

		for (int l = 0; l < fixedCount; ++l) {
			solveRHS[i][waveletCount + anchorCount + l] = anchorWeight * mMesh->getVertexPosition(fixedVerts[l])[i];
		}
	}
	g_engineWrapper.addColVec(solveRHS[0], "dcx");
	g_engineWrapper.addColVec(solveRHS[1], "dcy");
	g_engineWrapper.addColVec(solveRHS[2], "dcz");

	ZGeom::DenseMatrixd matOpt(matSGW);
	matOpt.expand(waveletCount + anchorCount + fixedCount, vertCount);
	for (int l = 0; l < anchorCount; ++l)
		matOpt(waveletCount + l, anchorIndex[l]) = anchorWeight;
	for (int l = 0; l < fixedCount; ++l)
		matOpt(waveletCount + anchorCount + l, fixedVerts[l]) = anchorWeight;
	g_engineWrapper.addDenseMat(matOpt, "matOpt");

	timer.stopTimer("Prepare deformation time: ");

	timer.startTimer();	
	//g_engineWrapper.eval("lsx=cgls(matOpt, dcx);");
	//g_engineWrapper.eval("lsy=cgls(matOpt, dcy);");
	//g_engineWrapper.eval("lsz=cgls(matOpt, dcz);");
	g_engineWrapper.eval("[lsx,flagx,resx]=lsqr(matOpt, dcx);");
	g_engineWrapper.eval("lsy=lsqr(matOpt, dcy);");
	g_engineWrapper.eval("lsz=lsqr(matOpt, dcz);");
	timer.stopTimer("Deformation time: ");

	double *lsx = g_engineWrapper.getDblVariablePtr("lsx");
	double *lsy = g_engineWrapper.getDblVariablePtr("lsy");
	double *lsz = g_engineWrapper.getDblVariablePtr("lsz");

	MeshCoordinates newCoord(vertCount, lsx, lsy, lsz);
	//MeshCoordinates newCoord(oldCoord);
	//newCoord.add(lsx, lsy, lsz);
	mMesh->setVertCoordinates(newCoord);

	evaluateApproximation(newCoord, "SGW deform");	
}

void ShapeEditor::reconstructSpectralWavelet()
{
	CStopWatch timer;	
	timer.startTimer();

	const int vertCount = mMesh->vertCount();	
	const ZGeom::DenseMatrixd& matW = mProcessor->getWaveletMat();
	if (matW.empty()) {
		mProcessor->computeSGW1();
		g_engineWrapper.addDenseMat(matW, "matSGW");
	}
	const int waveletCount = matW.rowCount();

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(waveletCount, 0);
	}
	g_engineWrapper.addColVec(solveRHS[0], "dcx");
	g_engineWrapper.addColVec(solveRHS[1], "dcy");
	g_engineWrapper.addColVec(solveRHS[2], "dcz");
	g_engineWrapper.addDenseMat(matW, "matOpt");

	timer.stopTimer("Prepare deformation time: ");

	timer.startTimer();	
	//g_engineWrapper.eval("lsx=cgls(matOpt, dcx);");
	//g_engineWrapper.eval("lsy=cgls(matOpt, dcy);");
	//g_engineWrapper.eval("lsz=cgls(matOpt, dcz);");
	g_engineWrapper.eval("[lsx,flagx,resx]=lsqr(matOpt, dcx);");
	g_engineWrapper.eval("lsy=lsqr(matOpt, dcy);");
	g_engineWrapper.eval("lsz=lsqr(matOpt, dcz);");
	timer.stopTimer("Deformation time: ");

	double *lsx = g_engineWrapper.getDblVariablePtr("lsx");
	double *lsy = g_engineWrapper.getDblVariablePtr("lsy");
	double *lsz = g_engineWrapper.getDblVariablePtr("lsz");

	MeshCoordinates newCoord;
	mMesh->getVertCoordinates(newCoord);
	newCoord.add(lsx, lsy, lsz);
	mMesh->setVertCoordinates(newCoord);

	evaluateApproximation(newCoord, "SGW reconstruction");	
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
	mMesh->getVertCoordinates(newCoord);
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
										 ZGeom::FunctionApproximation* vApproxCoeff[3], 
										 int nReconstruct, 
										 std::vector<MeshCoordinates>& continuousCoords, 
										 MeshCoordinates& finalCoord )
{
	const int vertCount = mMesh->vertCount();
	const ZGeom::FunctionApproximation &vApproxX = *vApproxCoeff[0], &vApproxY = *vApproxCoeff[1], &vApproxZ = *vApproxCoeff[2];

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
	DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, oldCoord, lCoord1);
	DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, newCoord, lCoord2);
	double dif2 = lCoord1.difference(lCoord2);
	std::cout << "-- Avg position error: " << dif1 / vertCount << '\n';
	std::cout << "-- Avg geometric error:    " << (dif1 + dif2)/2.0/vertCount << '\n';
}

void ShapeEditor::runTests()
{
	//monolithicApproximationTest1(true, true);
	//monolithicApproximationTest2(false, false);
	partitionedApproximationTest2();
	//partitionedApproximationTest1();
}

//// Test with graph Laplacian ////
//
void ShapeEditor::monolithicApproximationTest1( bool computeSGW, bool sgwCoding )
{
	using ZGeom::VecNd;
	std::cout << "\n--------" << " Starting MonoApproxTest1 (Graph Laplacian) "
		      << std::setfill('-') << std::setw(8) << '\n';
	
	CStopWatch timer;
	/* prepare prerequisite data */
	const int vertCount = mMesh->vertCount();
	const MeshCoordinates& oldCoord = getOldMeshCoord();	
	const ZGeom::EigenSystem &graphMHB = mProcessor->getMHB(Umbrella);

	std::vector<double> wDiag;
	mProcessor->getMeshLaplacian(Umbrella).getW().getDiagonal(wDiag);
	ZGeom::InnerProdcutFunc innerProdDiagW = [&](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2) {
		double *y = new double[vertCount];
		vdmul(&vertCount, v1.c_ptr(), &wDiag[0], &y[0]); 
		double res = cblas_ddot(vertCount, y, 1, v2.c_ptr(), 1);
		delete []y;
		return res;
	};

	const ZGeom::VecNd &vx = oldCoord.getCoordFunc(0), &vy = oldCoord.getCoordFunc(1), &vz = oldCoord.getCoordFunc(2);
	std::vector<ZGeom::VecNd> vSignals;
	vSignals.push_back(vx); vSignals.push_back(vy); vSignals.push_back(vz);

	ZGeom::FunctionApproximation vApproxX, vApproxY, vApproxZ;
	std::vector<ZGeom::FunctionApproximation*> vApproxCoeff;
	vApproxCoeff.push_back(&vApproxX); vApproxCoeff.push_back(&vApproxY); vApproxCoeff.push_back(&vApproxZ);

	/************************************************************************/
	/*  In each approximation test, the following is required
	/*    (0) Clear previously computed results;
	/*    (1) Compute atoms (eigenbasis, wavelet, etc.);
	/*    (2) Compute multi-channel approximation coefficients;
	/*    (3) Reconstruct with the computed coefficients; Evaluate results	                                        */
	/************************************************************************/ 

	const int nEigTotal = graphMHB.eigVecCount();
	int nAtomSel = 100;
	std::vector<VecNd>& vAtoms = mAtoms;

	/* Test 1, Fourier approximation */
	{	
		std::cout << '\n';
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();
		for (int i = 0; i < nEigTotal; ++i) 
			vAtoms.push_back(graphMHB.getEigVec(i));
		GeneralizedSimultaneousFourierApprox(vSignals, vAtoms, nAtomSel, vApproxCoeff, innerProdDiagW);
		computeApproximations(vAtoms, &vApproxCoeff[0], nAtomSel, mContReconstructCoords[0], getStoredCoordinate(0));
		//std::cout << "Reconstruct error (1): " << oldCoord.difference(mCoords[1]) << "\n";
		evaluateApproximation(getStoredCoordinate(0), "1");
	}	
	//////////////////////////////////////////////////////////////////////////


	/* Test 2, Simultaneous Fourier Matching Pursuit */
	{
		std::cout << '\n';
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();		
		for (int i = 0; i < nEigTotal; ++i) 
			vAtoms.push_back(graphMHB.getEigVec(i));
		timer.startTimer();
		ZGeom::GeneralizedSimultaneousMP(vSignals, vAtoms, nAtomSel, vApproxCoeff, innerProdDiagW, 2.);
		timer.stopTimer("Time to compute Fourier SMP: ");
		computeApproximations(vAtoms, &vApproxCoeff[0], nAtomSel, mContReconstructCoords[1], getStoredCoordinate(1));
		
		std::vector<int> vSelectedAtomIdx = vApproxX.getAllAtomIndex();	
		int selectedFromTop = std::count_if(vSelectedAtomIdx.begin(), vSelectedAtomIdx.end(), [&](int idx) { return idx < vSelectedAtomIdx.size();} );
		std::cout << "Ratio of MP overlapping: " << (double)selectedFromTop / (double)vSelectedAtomIdx.size() << '\n';
		evaluateApproximation(getStoredCoordinate(1), "2");
	}

	//////////////////////////////////////////////////////////////////////////
	
	/* Test 3, Wavelet Simultaneous-OMP */
	if (computeSGW) 
	{
		std::cout << '\n';
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();	
		std::cout << "To compute wavelet matching pursuit..\n";
		ZGeom::DenseMatrixd& matSGW = mProcessor->getWaveletMat();
		mProcessor->computeSGW1(Umbrella);	
		mTotalScales = matSGW.rowCount() / vertCount;
		//mProcessor->computeMixedAtoms1(CotFormula);
		for (int i = 0; i < matSGW.rowCount(); ++i) {
			ZGeom::VecNd newBasis = matSGW.getRowVec(i);
			newBasis.normalize(ZGeom::RegularProductFunc);
			vAtoms.push_back(newBasis);
		}

		if (sgwCoding) {
			nAtomSel = 100;
			timer.startTimer();
			ZGeom::SimultaneousOMP(vSignals, vAtoms, nAtomSel, vApproxCoeff, 2);
			timer.stopTimer("Time to compute Wavelet SOMP: ");

			std::ofstream ofcoding("output/coding1.txt");
			for (int i = 0; i < nAtomSel; ++i) {
				for (int c = 0; c < 3; ++c) {
					auto& sc = (*vApproxCoeff[c])[i];
					ofcoding << '(' << sc.index() << ',' << sc.coeff() << ") "; 
				}
				ofcoding << '\n';
			}

			computeApproximations(vAtoms, &vApproxCoeff[0], nAtomSel, mContReconstructCoords[2], getStoredCoordinate(2));
			evaluateApproximation(getStoredCoordinate(2), "3");
			changeCoordinates(2);
#if 0
			std::vector<int> vSelectedAtomIdx = vApproxX.getAllAtomIndex();	
			updateEditBasis(vAtoms, vSelectedAtomIdx);
			int nClass = vAtoms.size() / vertCount;
			std::vector<int> vScaleStat(nClass, 0);
			for (int idx : vSelectedAtomIdx) {
				vScaleStat[idx/vertCount] += 1;
			}
			std::cout << "Selected atom stats: \n";
			for (int s = 0; s < nClass; ++s) 
				std::cout << "  scale " << s << ": " << vScaleStat[s] << '\n';
			std::cout << "    total: " << vSelectedAtomIdx.size() << "\n";
			MeshFeatureList *mfl = new MeshFeatureList;
			for (int i = 0; i < 50; ++i) {
				int atomIdx = vSelectedAtomIdx[i];
				mfl->addFeature(new MeshFeature(atomIdx % vertCount, atomIdx / vertCount));
			}
			mfl->setIDandName(FEATURE_SGW_SOMP, "Feature_SGW_SOMP");
			mProcessor->addProperty(mfl);
			mProcessor->setActiveFeaturesByID(FEATURE_SGW_SOMP);
#endif
		}
	} //end of wavelet OMP
	//////////////////////////////////////////////////////////////////////////
	
	std::cout << "MonoApproxTest1 completed!\n";
	std::cout << std::setfill('+') << std::setw(40) << '\n';
}

//// test with CotFormula Laplacian ////
//
void ShapeEditor::monolithicApproximationTest2( bool doWavelet /*= true*/ )
{
	using ZGeom::VecNd;
	std::cout << "\n--------" << " Starting MonoApproxTest2 (Cot Laplacian) "
		      << std::setfill('-') << std::setw(8) << '\n';
	
	CStopWatch timer;
	/* prepare prerequisite data */
	const int vertCount = mMesh->vertCount();
	const MeshCoordinates& oldCoord = getOldMeshCoord();
	const ZGeom::EigenSystem &cotMHB = mProcessor->getMHB(CotFormula);

	std::vector<double> wDiag;
	mProcessor->getMeshLaplacian(CotFormula).getW().getDiagonal(wDiag);
	ZGeom::InnerProdcutFunc innerProdDiagW = [&](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2) {
		double *y = new double[vertCount];
		vdmul(&vertCount, v1.c_ptr(), &wDiag[0], &y[0]); 
		double res = cblas_ddot(vertCount, y, 1, v2.c_ptr(), 1);
		delete []y;
		return res;
	};

	const ZGeom::VecNd &vx = oldCoord.getCoordFunc(0), &vy = oldCoord.getCoordFunc(1), &vz = oldCoord.getCoordFunc(2);
	std::vector<ZGeom::VecNd> vSignals;
	vSignals.push_back(vx); vSignals.push_back(vy); vSignals.push_back(vz);

	ZGeom::FunctionApproximation vApproxX, vApproxY, vApproxZ;
	std::vector<ZGeom::FunctionApproximation*> vApproxCoeff;
	vApproxCoeff.push_back(&vApproxX); vApproxCoeff.push_back(&vApproxY); vApproxCoeff.push_back(&vApproxZ);

	/************************************************************************/
	/*  In each approximation test, the following is required
	/*    (0) Clear previously computed results;
	/*    (1) Compute atoms (eigenbasis, wavelet, etc.);
	/*    (2) Compute multi-channel approximation coefficients;
	/*    (3) Reconstruct with the computed coefficients; Evaluate results	                                        */
	/************************************************************************/ 

	const int nEigTotal = cotMHB.eigVecCount();
	int nAtomSel = 100;
	int nReconstruct = 100;
	std::vector<VecNd>& vAtoms = mAtoms;

	/* Test 1, Fourier approximation */
	{		
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();

		for (int i = 0; i < nEigTotal; ++i) 
			vAtoms.push_back(cotMHB.getEigVec(i));

		GeneralizedSimultaneousFourierApprox(vSignals, vAtoms, nAtomSel, vApproxCoeff, innerProdDiagW);
		
		computeApproximations(vAtoms, &vApproxCoeff[0], nReconstruct, mContReconstructCoords[0], getStoredCoordinate(0));
		
		std::cout << "Reconstruct error (1): " << oldCoord.difference(getStoredCoordinate(0)) << "\n\n";
	}	
	//////////////////////////////////////////////////////////////////////////


	/* Test 2, Simultaneous Fourier Matching Pursuit */
	{
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();
		
		for (int i = 0; i < nEigTotal; ++i) 
			vAtoms.push_back(cotMHB.getEigVec(i));

		timer.startTimer();
		ZGeom::GeneralizedSimultaneousMP(vSignals, vAtoms, nAtomSel, vApproxCoeff, innerProdDiagW, 2.);
		timer.stopTimer("Time to compute Fourier SMP: ");

		computeApproximations(vAtoms, &vApproxCoeff[0], nReconstruct, mContReconstructCoords[1], getStoredCoordinate(1));
		
		std::vector<int> vSelectedAtomIdx = vApproxX.getAllAtomIndex();	
		int selectedFromTop = std::count_if(vSelectedAtomIdx.begin(), vSelectedAtomIdx.end(), [&](int idx) { return idx < vSelectedAtomIdx.size();} );
		std::cout << "Ratio of MP overlapping: " << (double)selectedFromTop / (double)vSelectedAtomIdx.size() << '\n';
		std::cout << "Reconstruct error (2): " << oldCoord.difference(getStoredCoordinate(1)) << "\n\n";
	}
	//////////////////////////////////////////////////////////////////////////


	/* Test 3, Wavelet Simultaneous-OMP */
	if (doWavelet) {
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();
	
		std::cout << "To compute wavelet matching pursuit..\n";
		ZGeom::DenseMatrixd& matSGW = mProcessor->getWaveletMat();
		mProcessor->computeSGW1(CotFormula);	
		mTotalScales = matSGW.rowCount() / vertCount;
		//mProcessor->computeMixedAtoms1(CotFormula);
		for (int i = 0; i < matSGW.rowCount(); ++i) {
			ZGeom::VecNd newBasis = matSGW.getRowVec(i);
			newBasis.normalize(ZGeom::RegularProductFunc);
			vAtoms.push_back(newBasis);
		}
#if 0
		nAtomSel = 100;
		timer.startTimer();
		ZGeom::SimultaneousOMP(vSignals, vAtoms, nAtomSel, vApproxCoeff, 2);
		timer.stopTimer("Time to compute Wavelet SOMP: ");

		//nReconstruct = nAtomSel;
		computeApproximations(vAtoms, &vApproxCoeff[0], nReconstruct, mContReconstructCoords[2], mCoords[3]);
		std::cout << "Reconstruct error (3): " << oldCoord.difference(mCoords[3]) << "\n\n";
		mApproxCoeff[2] = vApproxX;

		std::vector<int> vSelectedAtomIdx = vApproxX.getAllAtomIndex();	
		updateEditBasis(vAtoms, vSelectedAtomIdx);

		int nClass = vAtoms.size() / vertCount;
		std::vector<int> vScaleStat(nClass, 0);
		for (int idx : vSelectedAtomIdx) {
			vScaleStat[idx/vertCount] += 1;
		}
		std::cout << "Selected atom stats: \n";
		for (int s = 0; s < nClass; ++s) 
			std::cout << "  scale " << s << ": " << vScaleStat[s] << '\n';
		std::cout << "    total: " << vSelectedAtomIdx.size() << "\n\n"; 

		MeshFeatureList *mfl = new MeshFeatureList;
		for (int i = 0; i < 50; ++i) {
			int atomIdx = vSelectedAtomIdx[i];
			mfl->addFeature(new MeshFeature(atomIdx % vertCount, atomIdx / vertCount));
		}
		mfl->setIDandName(FEATURE_SGW_SOMP, "Feature_SGW_SOMP");
		mProcessor->addProperty(mfl);
		mProcessor->setActiveFeaturesByID(FEATURE_SGW_SOMP);
#endif
	} //end of wavelet OMP
	//////////////////////////////////////////////////////////////////////////
	
	std::cout << "MonoApproxTest2 completed!\n";
	std::cout << std::setfill('+') << std::setw(40) << '\n';
}

//// Test partitioned approximation with graph Laplacian ////
//
void ShapeEditor::partitionedApproximationTest1()
{
	std::function<void(std::string, char c)> printBeginSeparator = [&](std::string s, char c) {
		std::cout << '\n';
		for (int i = 0; i < 8; ++i) std::cout << c;
		std::cout << ' ' << s << ' ';
		for (int i = 0; i < 8; ++i) std::cout << c;
		std::cout << '\n';
	};
	std::function<void(char, int)> printEndSeparator = [&](char c, int num) {
		std::cout << std::setfill(c) << std::setw(num) << '\n';
	};

	revertCoordinates();
	std::cout << "\n======== Starting PartitionedApproxTest 1 ========\n";	

	const int totalVertCount = mMesh->vertCount();
	int eigenCount = -1;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	MeshCoordinates oldGeoCoord;
	DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, oldMeshCoord, oldGeoCoord);

	mShapeApprox.init(mMesh);
	mShapeApprox.doSegmentation(1000);
	mSegmentPalette.generatePalette(mShapeApprox.partitionCount());	
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrColorPartitions).getValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrColorPartitions.c_str()));
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
			DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
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

			DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
			double meshError = oldMeshCoord.difference(vProgressiveCoords[i]);
			double geoError = oldGeoCoord.difference(newGeoCoord);
			double mixedError = (meshError + geoError) / (2*totalVertCount);
			ofs << mixedError << ((i < ratiosCount-1) ? '\t' : '\n');
		}
	};

	double ratioOverhead;
	// 1. low-pass approximate, MHB
	printBeginSeparator("Low-pass Approx, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	compressAndEvaluate(maxCodingRatio, SA_Truncation, false);

	evaluateApproximation(vProgressiveCoords.back(), "1");
	setStoredCoordinates(vProgressiveCoords.back(), 0);
	mContReconstructCoords[0] = vProgressiveCoords;
	emit approxStepsChanged(0, vProgressiveCoords.size());
	printEndSeparator('-', 40);
		
	// 2. SOMP, MHB
	printBeginSeparator("SMP, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	compressAndEvaluate(maxCodingRatio, SA_SMP, true);

	evaluateApproximation(vProgressiveCoords.back(), "2");
	setStoredCoordinates(vProgressiveCoords.back(), 1);
	mContReconstructCoords[1] = vProgressiveCoords;
	emit approxStepsChanged(1, vProgressiveCoords.size());
	printEndSeparator('-', 40);

	//2.5 SOMP, MHB+Spikes
#if 0
	printBeginSeparator("SOMP, MHB-Spikes", '-');
	mShapeApprox.constructDictionaries(DT_FourierSpikes);
	pursuitAndEvaluate(SA_SOMP);

	evaluateApproximation(vProgressiveCoords.back(), "2.5");
	printEndSeparator('-', 40);
#endif

	// 3. SOMP, SGW
#if 1
	printBeginSeparator("SOMP, SGW", '-');
	mShapeApprox.constructDictionaries(DT_SGW4);
	pursuitAndEvaluate(SA_SOMP);

	evaluateApproximation(vProgressiveCoords.back(), "3");
	setStoredCoordinates(vProgressiveCoords.back(), 2);
	mContReconstructCoords[2] = vProgressiveCoords;
	emit approxStepsChanged(2, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	// 4. SOMP, SGW + MHB
#if 1
	printBeginSeparator("SOMP, SGW-MHB", '-');
	mShapeApprox.constructDictionaries(DT_SGW4MHB);
	pursuitAndEvaluate(SA_SOMP);
	
	evaluateApproximation(vProgressiveCoords.back(), "4");
	setStoredCoordinates(vProgressiveCoords.back(), 3);
	mContReconstructCoords[3] = vProgressiveCoords;
	emit approxStepsChanged(3, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	ofs.close();
	std::cout << '\n';
	changeCoordinates(1);
	std::cout << "** PartitionedApproxTest 1 completed!\n";
	printEndSeparator('=', 40);
}

void ShapeEditor::partitionedApproximationTest2()
{
	std::function<void(std::string, char c)> printBeginSeparator = [&](std::string s, char c) {
		std::cout << '\n';
		for (int i = 0; i < 8; ++i) std::cout << c;
		std::cout << ' ' << s << ' ';
		for (int i = 0; i < 8; ++i) std::cout << c;
		std::cout << '\n';
	};
	std::function<void(char, int)> printEndSeparator = [&](char c, int num) {
		std::cout << std::setfill(c) << std::setw(num) << '\n';
	};

	revertCoordinates();
	std::cout << "\n======== Starting PartitionedApproxTest 2 ========\n";	

	const int totalVertCount = mMesh->vertCount();
	int eigenCount = -1;
	const double maxCodingRatio = 1.0;
	const MeshCoordinates& oldMeshCoord = getOldMeshCoord();
	MeshCoordinates oldGeoCoord;
	DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, oldMeshCoord, oldGeoCoord);

	mShapeApprox.init(mMesh);
	mShapeApprox.doSegmentation(1000);
	mSegmentPalette.generatePalette(mShapeApprox.partitionCount());	
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrColorPartitions).getValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrColorPartitions.c_str()));
	mShapeApprox.doEigenDecomposition(Umbrella, eigenCount);	

	/* setup vectors of compress ratios and initialize vProgressiveCoords */
	double ratioToTest[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.4/*, 0.5, 0.6, 0.7, 0.8*/};
	std::vector<double> vCompressRatio(std::begin(ratioToTest), std::end(ratioToTest));
	//for (int i = 1; i <= 8; ++i) vCompressRatio.push_back(0.1 * i);
	std::vector<MeshCoordinates> vProgressiveCoords(vCompressRatio.size());
	std::ofstream ofs("output/approx_errors.txt");	
	for (int i = 0; i < vCompressRatio.size()-1; ++i) ofs << vCompressRatio[i] << '\t';
	ofs << vCompressRatio.back() << '\n';

	std::function<void(double,bool)> compressAndEvaluate = [&](double overhead, bool exploitSparsity) 
	{
		int ratioSteps = vCompressRatio.size();
		for (int i = 0; i < ratioSteps; ++i) {
			double basisRatio = vCompressRatio[i] - overhead;
			mShapeApprox.doSparseReconstructionByRatio(basisRatio, vProgressiveCoords[i], exploitSparsity);
		}

		MeshCoordinates newGeoCoord;
		for (int i = 0; i < ratioSteps; ++i) {			
			DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
			double meshError = oldMeshCoord.difference(vProgressiveCoords[i]);
			double geoError = oldGeoCoord.difference(newGeoCoord);
			double mixedError = (meshError + geoError) / (2*totalVertCount);
			ofs << mixedError << (i < ratioSteps-1) ? '\t' : '\n';
		}		
	};

	std::function<void(double)> multiPursuitAndEvaluate = [&](double overhead) 
	{
		int ratioSteps = vCompressRatio.size();
		for (int i = 0; i < ratioSteps; ++i) {
			mShapeApprox.findSparseRepresentationByRatio(SA_SOMP, vCompressRatio[i] - overhead, true);
			mShapeApprox.doSparseReconstructionBySize(-1, vProgressiveCoords[i]);

			MeshCoordinates newGeoCoord;
			DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
			double meshError = oldMeshCoord.difference(vProgressiveCoords[i]);
			double geoError = oldGeoCoord.difference(newGeoCoord);
			double mixedError = (meshError + geoError) / (2*totalVertCount);
			ofs << mixedError << (i < ratioSteps-1) ? '\t' : '\n';
		}
	};

	std::function<void(double)> multiPursuitAndEvaluate2 = [&](double overhead) 
	{
		int ratioSteps = vCompressRatio.size();
		for (int i = 0; i < ratioSteps; ++i) {
			mShapeApprox.findSparseRepresentationByRatio(SA_SMP, vCompressRatio[i] - overhead, true);
			mShapeApprox.doSparseReconstructionBySize(-1, vProgressiveCoords[i]);

			MeshCoordinates newGeoCoord;
			DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, vProgressiveCoords[i], newGeoCoord);
			double meshError = oldMeshCoord.difference(vProgressiveCoords[i]);
			double geoError = oldGeoCoord.difference(newGeoCoord);
			double mixedError = (meshError + geoError) / (2*totalVertCount);
			ofs << mixedError << '\t';
		}
		ofs << '\n';
	};

	double ratioOverhead;
	// 1. low-pass approximate, MHB
	printBeginSeparator("Low-pass Approx, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	mShapeApprox.findSparseRepresentationByRatio(SA_Truncation, maxCodingRatio, false);	
	ratioOverhead = 0.;
	compressAndEvaluate(ratioOverhead, false);

	evaluateApproximation(vProgressiveCoords.back(), "1");
	setStoredCoordinates(vProgressiveCoords.back(), 0);
	mContReconstructCoords[0] = vProgressiveCoords;
	emit approxStepsChanged(0, vProgressiveCoords.size());
	printEndSeparator('-', 40);

	// 2. SOMP, MHB
	printBeginSeparator("SOMP, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	mShapeApprox.findSparseRepresentationByRatio(SA_SMP, maxCodingRatio, false);
	ratioOverhead = 1.0/96.0;
	compressAndEvaluate(ratioOverhead, true);

	evaluateApproximation(vProgressiveCoords.back(), "2");
	setStoredCoordinates(vProgressiveCoords.back(), 1);
	mContReconstructCoords[1] = vProgressiveCoords;
	emit approxStepsChanged(1, vProgressiveCoords.size());
	printEndSeparator('-', 40);

	// 3. SOMP, SGW
#if 0
	printBeginSeparator("SOMP, SGW3", '-');
	mShapeApprox.constructDictionaries(DT_SGW3);
	ratioOverhead = 4.0/96.0;
	multiPursuitAndEvaluate(ratioOverhead);
#endif
#if 1
	printBeginSeparator("SOMP, SGW4", '-');
	mShapeApprox.constructDictionaries(DT_SGW4);
	ratioOverhead = 5.0/96.0;
	multiPursuitAndEvaluate(ratioOverhead);
#endif
#if 0
	printBeginSeparator("SOMP, SGW5", '-');
	mShapeApprox.constructDictionaries(DT_SGW5);
	ratioOverhead = 6.0/96.0;
	multiPursuitAndEvaluate(ratioOverhead);
#endif
#if 1
	evaluateApproximation(vProgressiveCoords.back(), "3");
	setStoredCoordinates(vProgressiveCoords.back(), 2);
	mContReconstructCoords[2] = vProgressiveCoords;
	emit approxStepsChanged(2, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	// 4. SOMP, SGW + MHB
#if 0
	printBeginSeparator("SOMP, SGW3-MHB", '-');
	mShapeApprox.constructDictionaries(DT_SGW3MHB);
	ratioOverhead = 5.0/96.0;
	multiPursuitAndEvaluate(ratioOverhead);
#endif
#if 1
	printBeginSeparator("SOMP, SGW4-MHB", '-');
	mShapeApprox.constructDictionaries(DT_SGW4MHB);
	ratioOverhead = 6.0/96.0;
	multiPursuitAndEvaluate(ratioOverhead);
#endif
#if 0
	printBeginSeparator("SOMP, SGW5-MHB", '-');
	mShapeApprox.constructDictionaries(DT_SGW5MHB);
	ratioOverhead = 7.0/96.0;
	multiPursuitAndEvaluate(ratioOverhead);
#endif
#if 1	
	evaluateApproximation(vProgressiveCoords.back(), "4");
	setStoredCoordinates(vProgressiveCoords.back(), 3);
	mContReconstructCoords[3] = vProgressiveCoords;
	emit approxStepsChanged(3, vProgressiveCoords.size());
	printEndSeparator('-', 40);
#endif

	// 5. SOMP, SGW + MHB
#if 0
	printBeginSeparator("SMP, SGW4-MHB", '-');
	mShapeApprox.constructDictionaries(DT_SGW4MHB);
	ratioOverhead = 6.0/96.0;
	multiPursuitAndEvaluate2(ratioOverhead);

	evaluateApproximation(vProgressiveCoords.back(), "5");
	setStoredCoordinates(vProgressiveCoords.back(), 4);
	mContReconstructCoords[4] = vProgressiveCoords;
	printEndSeparator('-', 40);
#endif

	ofs.close();
	std::cout << '\n';
	changeCoordinates(1);
	std::cout << "** PartitionedApproxTest 2 completed!\n";
	printEndSeparator('=', 40);
}


/////// compute various eigenvectors indexed by Fiedler vector ////////////////
////
void ShapeEditor::spectrumTest1()	
{
	const int vertCount = mMesh->vertCount();
	LaplacianType lapType = SymCot;
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ZGeom::EigenSystem &mhb = mProcessor->getMHB(lapType);
	const int nEig = mhb.eigVecCount();

	using ZGeom::VecNd;
	const VecNd& fiedlerVec = mhb.getEigVec(1);
	std::vector<std::pair<int, double> > vertFiedler;
	for (int i = 0; i < vertCount; ++i) vertFiedler.push_back(std::make_pair(i, fiedlerVec[i]));

	std::sort(vertFiedler.begin(), vertFiedler.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) { return p1.second < p2.second; });
	
	std::ofstream ofs("output/SortedEigVec.txt");
	for (auto p : vertFiedler) {
		ofs << p.first;
		for (int ek = 1; ek <= 10; ek += 1) {
			ofs << ' ' << mhb.getEigVec(ek)[p.first];
		}
		ofs << std::endl;
	}
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
