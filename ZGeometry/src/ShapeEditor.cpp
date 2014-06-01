#include "ShapeEditor.h"
#include <iostream>
#include <random>
#include <functional>
#include <iomanip>
#include <ppl.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/Approximation.h>
#include <ZGeom/util.h>
#include "global.h"

using ZGeom::Dictionary;

Vector3D toVector3D(const ZGeom::Vec3d& v) { return Vector3D(v[0], v[1], v[2]); }
ZGeom::Vec3d toVec3d(const Vector3D& v) { return ZGeom::Vec3d(v.x, v.y, v.z); }

bool readPursuits(const std::string& pursuitFile, ZGeom::FunctionApproximation *vPursuits[]) 
{
	if (!fileExist(pursuitFile)) 
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
	processor->getMesh()->retrieveVertCoordinates(mOriginalCoord); 
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
	qout.output("Current coordinate ID: " + Int2String(mCurCoordID), OUT_STATUS);

	if (mStoredCoordinates[mCurCoordID].empty())
		qout.output("!! Selected coordinate is empty!", OUT_CONSOLE);
	else 
		mMesh->setVertCoordinates(mStoredCoordinates[mCurCoordID]);
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
	mMesh->retrieveVertCoordinates(oldCoord);
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
	mMesh->retrieveVertCoordinates(newCoord);
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

void ShapeEditor::runTests()
{
	//sparseCompressionTest();		// compression test
	//approximationTest2();
	sparseDecompositionTest();
}

//// compute various eigenvectors indexed by Fiedler vector /////
//
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

//// Test partitioned approximation with graph Laplacian ////
//
void ShapeEditor::sparseCompressionTest()
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
	DifferentialMeshProcessor::computeGeometricLaplacianCoordinate(*mMesh, oldMeshCoord, oldGeoCoord);

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

	// 1. low-pass approximate, MHB
#if 1
	printBeginSeparator("Low-pass Approx, MHB", '-');
	mShapeApprox.constructDictionaries(DT_Fourier);
	compressAndEvaluate(maxCodingRatio, SA_Truncation, false);

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
	pursuitAndEvaluate(SA_SOMP);
	
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
void ShapeEditor::sparseFeatureFindingTest1()
{
	printBeginSeparator("Starting ApproxTest2 (feature analysis)",'=');
	CStopWatch timer;

	LaplacianType lapType;
	DictionaryType dictType;
	SparseApproxMethod approxMethod;
	lapType = CotFormula;
	lapType = Umbrella;
	dictType = DT_SGW4;
	approxMethod = SA_SOMP; 
	//approxMethod = SA_LASSO; 
	
	// initializing, segmentation, coloring, and eigendecomposition
	int totalVertCount = mMesh->vertCount();
	int eigenCount = min(600, totalVertCount-1); // -1 means vertCount -1 
	int maxPatchSize = 1000;					 // -1 means no segmentation
	mShapeApprox.init(mMesh);
	mShapeApprox.doSegmentation(maxPatchSize);
	mSegmentPalette.generatePalette(mShapeApprox.partitionCount());	
	std::vector<ZGeom::Colorf>& vColors = mMesh->addColorAttr(StrAttrColorPartitions).attrValue();
	colorPartitions(mShapeApprox.mPartIdx, mSegmentPalette, vColors);
	emit signatureComputed(QString(StrAttrColorPartitions.c_str()));
	mShapeApprox.doEigenDecomposition(lapType, eigenCount);	

	// possible signals: coordinate or curvature functions
	MeshCoordinates meshCoord = mMesh->getVertCoordinates();
	std::vector<double> vMeanCurvatures = mMesh->getMeanCurvature();
	std::vector<double> vGaussCurvature = mMesh->getGaussCurvature();

	// construct dictionary
	timer.startTimer();
	mShapeApprox.mSubMeshApprox[0].constructDict(dictType);
	timer.stopTimer("Time to construct dictionary: ", "s");

	{ /** To print signals and dictionary data **/
		std::ofstream ofs;
		ofs.open("output/signals.txt");
		for(int i = 0; i < totalVertCount; ++i) {
			ofs << meshCoord.getXCoord()[i] << ' ' << meshCoord.getYCoord()[i] << ' ' 
				<< meshCoord.getZCoord()[i] << ' ' << vMeanCurvatures[i] << ' ' 
				<< vGaussCurvature[i] << '\n';
		}
		ofs.close();
		ofs.open("output/dictionary.txt");
		const ZGeom::Dictionary& dict = mShapeApprox.mSubMeshApprox[0].getDict();
		for (int i = 0; i < totalVertCount; ++i) {
			for (int k = 0; k < dict.size(); ++k) 
				ofs << dict[k][i] << ((k < dict.size() - 1) ? ' ' : '\n');
		}	
		ofs.close();
	}/** End of printing signals and dictionary data **/

	// compute sparse coding
	int dictSize = mShapeApprox.mSubMeshApprox[0].dictSize();
	ZGeom::FunctionApproximation vCoeff;
	//std::vector<double>& vSignal = meshCoord.getXCoord().toStdVector();
	std::vector<double>& vSignal = vMeanCurvatures;

	SparseCodingOptions opts;
	opts.mApproxMethod = approxMethod;
	opts.mCodingAtomCount = 323;
	opts.lambda1 = 0.15;

	timer.startTimer();
	mShapeApprox.mSubMeshApprox[0].computeSparseCoding(vSignal, opts, vCoeff);
	timer.stopTimer("Time to do sparse coding: ", "s");

	ZGeom::VecNd vReconstruct = ReconstructApproximationSingleChannel(mShapeApprox.mSubMeshApprox[0].getDict(), vCoeff);
	double residual = (ZGeom::VecNd(vSignal) - vReconstruct).norm2();	
	std::cout << "Approximation residual (" << vCoeff.size() << " basis): " << residual << '\n';
	
	std::vector<std::pair<int,double> > vPursuit;
	for (int i = 0; i < vCoeff.size(); ++i) 
		vPursuit.push_back(std::make_pair(vCoeff[i].index(), vCoeff[i].coeff()));
	std::sort(vPursuit.begin(), vPursuit.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2)->bool{
		return std::fabs(p1.second) > std::fabs(p2.second);
	});
	//for (auto p : vPursuit) std::cout << p.first + 1 << ' ' << p.second << '\n';

	// analysis of coding
	int nScales = dictSize / totalVertCount;
	std::vector<int> vAtomScaleCount(nScales, 0);
	for (auto c : vCoeff.getApproxItems()) {
		vAtomScaleCount[c.index() / totalVertCount]++;
	}
	for (int i = 0; i < nScales; ++i)
		std::cout << "scale " << i << ": " << vAtomScaleCount[i] << '\n';

	MeshFeatureList* vSparseFeatures = new MeshFeatureList;
	std::ofstream ofs("output/sparse_coding.csv");
	for (auto c : vCoeff.getApproxItems()) {
		int scl = c.index() / totalVertCount, idx = c.index() % totalVertCount;
		double coef = c.coeff();
		ofs << scl << ", " << idx << ", " << coef << '\n';
		
		vSparseFeatures->addFeature(new MeshFeature(idx, scl));
		vSparseFeatures->back()->m_scalar1 = coef;
	}
	ofs.close();
	mMesh->addAttrMeshFeatures(*vSparseFeatures, StrAttrFeatureSparseSGW);
	
	printEndSeparator('=', 40);
}

//// Test shape decomposition via signal separation
//
void ShapeEditor::sparseDecompositionTest()
{
	/************************************************************************/
	/* 1. Support decomposition w or w/o segmentation                       */
	/* 2. Save the multilevel decomposition magnitude as color signatures   */
	/* 3. Optimized signal separation against MHB and SGW dictionary        */
	/*                                                                      */
	/************************************************************************/
	 
	revertCoordinates();
	mStoredCoordinates.resize(5);
	std::cout << "\n======== Starting ApproxTest3 (Decomposition test) ========\n";	

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
	SparseApproxMethod saMethod = SA_SOMP;

	MeshCoordinates reconstructedMeshCoord;
	
	mShapeApprox.constructDictionaries(dictType); 
	mShapeApprox.findSparseRepresentationBySize(saMethod, codingSize);
	mShapeApprox.doSparseReconstructionBySize(-1, reconstructedMeshCoord);
	
	std::cout << "Reconstruction error: " << oldMeshCoord.difference(reconstructedMeshCoord) << '\n';
	setStoredCoordinates(reconstructedMeshCoord, 1);
	

	const ZGeom::Dictionary& dict = mShapeApprox.mSubMeshApprox[0].getDict();
	const std::vector<ZGeom::ApproxItem>* vCoeff[3] = {
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
			const ZGeom::ApproxItem& sc = (*vCoeff[c])[i];
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
		const ZGeom::ApproxItem& sc0 = (*vCoeff[0])[i];
		const ZGeom::ApproxItem& sc1 = (*vCoeff[1])[i];
		const ZGeom::ApproxItem& sc2 = (*vCoeff[2])[i];
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

