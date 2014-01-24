#include "ShapeEditor.h"
#include <iostream>
#include <random>
#include <functional>
#include <ppl.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/Approximation.h>
#include <ZUtil/timer.h>
#include <ZUtil/zassert.h>
#include <ZUtil/zutil_io.h>
#include "global.h"


Vector3D toVector3D(const ZGeom::Vec3d& v) { return Vector3D(v[0], v[1], v[2]); }
ZGeom::Vec3d toVec3d(const Vector3D& v) { return ZGeom::Vec3d(v.x, v.y, v.z); }

bool readPursuits(const std::string& pursuitFile, ZGeom::FunctionApproximation *vPursuits[]) {
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


void ShapeEditor::init( DifferentialMeshProcessor* processor )
{
	mProcessor = processor; 
	mMesh = processor->getMesh();
	mEngine = processor->getMatlabEngineWrapper();
	
	mTotalScales = 0;
	mCoordSelect = 0;
	mCoords.resize(4);
	processor->getMesh()->getVertCoordinates(mCoords[0]); 

	std::cout << "Shape editor is initialized!" << std::endl;

	//addNoise(0.1);
	reconstructionTest2();
	
	mCoordSelect = mCoords[3].empty() ? 0 : 3;
	mMesh->setVertCoordinates(mCoords[mCoordSelect]);

	//editTest2();
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

void ShapeEditor::revertCoordinates()
{
	mCoordSelect = 0;
	std::cout << "Selected coordinate: " << mCoordSelect << std::endl;
	mMesh->setVertCoordinates(mCoords[mCoordSelect]);
}

void ShapeEditor::changeCoordinates()
{
	int oldCoordSelect = mCoordSelect;
	mCoordSelect = (mCoordSelect + 1) % mCoords.size();
	std::cout << "Selected coordinate: " << mCoordSelect << std::endl;
	
	if (mCoords[mCoordSelect].empty()) {
		std::cout << "Selected coordinate is empty!" << std::endl;
	} else {
		mMesh->setVertCoordinates(mCoords[mCoordSelect]);
	}
}

void ShapeEditor::continuousReconstruct( int selected, int atomCount )
{
	if (selected < 0 || selected >= 3 || atomCount < 0 || atomCount >= mContReconstructCoords[selected].size()) return;
	mMesh->setVertCoordinates(mContReconstructCoords[selected][atomCount]);

	if (selected == 2) {
		const int vertCount = mMesh->vertCount();
		MeshProperty* curAtomLocation = mProcessor->retrievePropertyByID(FEATURE_SGW_SOMP);
		if (NULL == curAtomLocation) return;
		MeshFeatureList *mfl = dynamic_cast<MeshFeatureList*>(curAtomLocation);
		mfl->clear();
		int atomIdx = mApproxCoeff[selected][atomCount].index();
		mfl->addFeature(new MeshFeature(atomIdx % vertCount, atomIdx / vertCount));
	}
}

void ShapeEditor::manifoldHarmonicsReconstruct( int nEig )
{
	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	MeshCoordinates newCoord(vertCount);
	LaplacianType lapType = SymCot;
		
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ManifoldHarmonics &mhb = mProcessor->getMHB(lapType);
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
	mEngine->addColVec(solveRHS[0].c_ptr(), vertCount + anchorCount, "dcx");
	mEngine->addColVec(solveRHS[1].c_ptr(), vertCount + anchorCount, "dcy");
	mEngine->addColVec(solveRHS[2].c_ptr(), vertCount + anchorCount, "dcz");

	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount, vertCount);
	matOptS.copyElements(matLs);
	for (int a = 0; a < anchorCount; ++a) 
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);
	mEngine->addSparseMat(matOptS, "matOptS");
	
	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	mEngine->eval("lsx=matOptS\\dcx;");
	mEngine->eval("lsy=matOptS\\dcy;");
	mEngine->eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = mEngine->getDblVariablePtr("lsx");
	double *lsy = mEngine->getDblVariablePtr("lsy");
	double *lsz = mEngine->getDblVariablePtr("lsz");

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
	mEngine->addColVec(solveRHS[0].c_ptr(), vertCount + anchorCount, "dcx");
	mEngine->addColVec(solveRHS[1].c_ptr(), vertCount + anchorCount, "dcy");
	mEngine->addColVec(solveRHS[2].c_ptr(), vertCount + anchorCount, "dcz");

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

	mEngine->addSparseMat(matOptS, "matOptS");

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	mEngine->eval("lsx=matOptS\\dcx;");
	mEngine->eval("lsy=matOptS\\dcy;");
	mEngine->eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = mEngine->getDblVariablePtr("lsx");
	double *lsy = mEngine->getDblVariablePtr("lsy");
	double *lsz = mEngine->getDblVariablePtr("lsz");

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
	mEngine->addColVec(solveRHS[0].c_ptr(), vertCount + anchorCount, "dcx");
	mEngine->addColVec(solveRHS[1].c_ptr(), vertCount + anchorCount, "dcy");
	mEngine->addColVec(solveRHS[2].c_ptr(), vertCount + anchorCount, "dcz");

	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount, vertCount);
	matOptS.copyElements(matBiL);
	for (int a = 0; a < anchorCount; ++a) 
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);
	mEngine->addSparseMat(matOptS, "matOptS");

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	mEngine->eval("lsx=matOptS\\dcx;");
	mEngine->eval("lsy=matOptS\\dcy;");
	mEngine->eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = mEngine->getDblVariablePtr("lsx");
	double *lsy = mEngine->getDblVariablePtr("lsy");
	double *lsz = mEngine->getDblVariablePtr("lsz");

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
	mEngine->addColVec(oldCoord.getXCoord(), "ecx");
	mEngine->addColVec(oldCoord.getYCoord(), "ecy");
	mEngine->addColVec(oldCoord.getZCoord(), "ecz");

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(SymCot).getLS();
	ZGeom::SparseMatrix<double> matBiL;
	ZGeom::mulMatMat(matLs, matLs, matBiL);
	mEngine->addSparseMat(matLs, "matL");
	mEngine->addSparseMat(matBiL, "matBiL");
	
	ZGeom::SparseMatrix<double> matL1(matLs);
	matL1.scale(-ks);
	ZGeom::SparseMatrix<double> matMixedL;
	ZGeom::addMatMat(matL1, matBiL, kb, matMixedL);	//matMixedL = -ks*matL + kb * matBiL
	mEngine->addSparseMat(matMixedL, "matMixedL");

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount + fixedCount, 0);
		for (int l = 0; l < anchorCount; ++l) {
			const Vector3D& oldPos = mMesh->getVertexPosition(anchorIndex[l]);
			solveRHS[i][vertCount + l] = anchorWeight * (anchorPos[l][i] - oldPos[i]);
		}
	}
	mEngine->addColVec(solveRHS[0], "dcx");
	mEngine->addColVec(solveRHS[1], "dcy");
	mEngine->addColVec(solveRHS[2], "dcz");

	ZGeom::SparseMatrix<double> matOptS(vertCount + anchorCount + fixedCount, vertCount);
	matOptS.copyElements(matMixedL);
	for (int a = 0; a < anchorCount; ++a) 
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);
	for (int a = 0; a < fixedCount; ++a)
		matOptS.insertElem(vertCount + anchorCount + a + 1, fixedVerts[a] + 1, anchorWeight);
	mEngine->addSparseMat(matOptS, "matOptS");

	timer.stopTimer("Prepare deformation time: ");
	timer.startTimer();	
	mEngine->eval("lsx=matOptS\\dcx;");
	mEngine->eval("lsy=matOptS\\dcy;");
	mEngine->eval("lsz=matOptS\\dcz;");
	timer.stopTimer("Deformation time: ");

	double *lsx = mEngine->getDblVariablePtr("lsx");
	double *lsy = mEngine->getDblVariablePtr("lsy");
	double *lsz = mEngine->getDblVariablePtr("lsz");

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
	mEngine->addColVec(oldCoord.getXCoord(), "ecx");
	mEngine->addColVec(oldCoord.getYCoord(), "ecy");
	mEngine->addColVec(oldCoord.getZCoord(), "ecz");
	
	ZGeom::DenseMatrixd& matSGW = mProcessor->getWaveletMat();
	if (matSGW.empty()) mProcessor->computeSGW1(CotFormula);
	const int waveletCount = matSGW.rowCount();
	
	/*std::vector<double> diagW;
	mProcessor->getMeshLaplacian(CotFormula).getW().getDiagonal(diagW);
	Concurrency::parallel_for(0, waveletCount, [&](int i){
		double *pr = matSGW.raw_ptr() + vertCount * i;
		for (int j = 0; j < vertCount; ++j) pr[j] *= diagW[j];
	});*/

	mEngine->addDenseMat(matSGW, "matSGW");
	
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
	mEngine->addColVec(solveRHS[0], "dcx");
	mEngine->addColVec(solveRHS[1], "dcy");
	mEngine->addColVec(solveRHS[2], "dcz");

	ZGeom::DenseMatrixd matOpt(matSGW);
	matOpt.expand(waveletCount + anchorCount + fixedCount, vertCount);
	for (int l = 0; l < anchorCount; ++l)
		matOpt(waveletCount + l, anchorIndex[l]) = anchorWeight;
	for (int l = 0; l < fixedCount; ++l)
		matOpt(waveletCount + anchorCount + l, fixedVerts[l]) = anchorWeight;
	mEngine->addDenseMat(matOpt, "matOpt");

	timer.stopTimer("Prepare deformation time: ");

	timer.startTimer();	
	//mEngine->eval("lsx=cgls(matOpt, dcx);");
	//mEngine->eval("lsy=cgls(matOpt, dcy);");
	//mEngine->eval("lsz=cgls(matOpt, dcz);");
	mEngine->eval("[lsx,flagx,resx]=lsqr(matOpt, dcx);");
	mEngine->eval("lsy=lsqr(matOpt, dcy);");
	mEngine->eval("lsz=lsqr(matOpt, dcz);");
	timer.stopTimer("Deformation time: ");

	double *lsx = mEngine->getDblVariablePtr("lsx");
	double *lsy = mEngine->getDblVariablePtr("lsy");
	double *lsz = mEngine->getDblVariablePtr("lsz");

	MeshCoordinates newCoord(vertCount, lsx, lsy, lsz);
	//MeshCoordinates newCoord(oldCoord);
	//newCoord.add(lsx, lsy, lsz);
	mMesh->setVertCoordinates(newCoord);

	evalReconstruct(newCoord);
}

void ShapeEditor::reconstructSpectralWavelet()
{
	CStopWatch timer;	
	timer.startTimer();

	const int vertCount = mMesh->vertCount();	
	const ZGeom::DenseMatrixd& matW = mProcessor->getWaveletMat();
	if (matW.empty()) {
		mProcessor->computeSGW1();
		mEngine->addDenseMat(matW, "matSGW");
	}
	const int waveletCount = matW.rowCount();

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(waveletCount, 0);
	}
	mEngine->addColVec(solveRHS[0], "dcx");
	mEngine->addColVec(solveRHS[1], "dcy");
	mEngine->addColVec(solveRHS[2], "dcz");
	mEngine->addDenseMat(matW, "matOpt");

	timer.stopTimer("Prepare deformation time: ");

	timer.startTimer();	
	//mEngine->eval("lsx=cgls(matOpt, dcx);");
	//mEngine->eval("lsy=cgls(matOpt, dcy);");
	//mEngine->eval("lsz=cgls(matOpt, dcz);");
	mEngine->eval("[lsx,flagx,resx]=lsqr(matOpt, dcx);");
	mEngine->eval("lsy=lsqr(matOpt, dcy);");
	mEngine->eval("lsz=lsqr(matOpt, dcz);");
	timer.stopTimer("Deformation time: ");

	double *lsx = mEngine->getDblVariablePtr("lsx");
	double *lsy = mEngine->getDblVariablePtr("lsy");
	double *lsz = mEngine->getDblVariablePtr("lsz");

	MeshCoordinates newCoord;
	mMesh->getVertCoordinates(newCoord);
	newCoord.add(lsx, lsy, lsz);
	mMesh->setVertCoordinates(newCoord);

	evalReconstruct(newCoord);
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

void ShapeEditor::evalReconstruct( const MeshCoordinates& newCoord ) const
{
	double errorSum(0);
	const int vertCount = mMesh->vertCount();
	for (int i = 0; i < vertCount; ++i) {
		errorSum += (newCoord[i] - getOldMeshCoord()[i]).length();
	}

	double avgError = errorSum / vertCount / mMesh->getAvgEdgeLength();
	std::cout << "Average reconstruct error: " << avgError << std::endl;
}

void ShapeEditor::reconstructionTest1()
{
	const int vertCount = mMesh->vertCount();
	LaplacianType lapType = CotFormula;
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ManifoldHarmonics &mhb = mProcessor->getMHB(lapType);
	const int nEig = mhb.eigVecCount();
	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	const ZGeom::VecNd &vx = oldCoord.getCoordFunc(0), &vy = oldCoord.getCoordFunc(1), &vz = oldCoord.getCoordFunc(2);
	std::vector<double> wDiag;
	mProcessor->getMeshLaplacian(lapType).getW().getDiagonal(wDiag);

	ZGeom::InnerProdcutFunc innerProdMatW = [&](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2) { return ZGeom::innerProductSym(v1, matW, v2); };
		
	ZGeom::InnerProdcutFunc innerProdDiagW = [&](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2) {
		double *y = new double[vertCount];
		vdmul(&vertCount, v1.c_ptr(), &wDiag[0], &y[0]);
		//for (int i = 0; i < vertCount; ++i) y[i] = v1[i] * wDiag[i];
		double res = cblas_ddot(vertCount, y, 1, v2.c_ptr(), 1);
		delete []y;
		return res;
	};

	ZGeom::InnerProdcutFunc innerProdRegular = [=](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2) {
		return cblas_ddot(vertCount, v1.c_ptr(), 1, v2.c_ptr(), 1);
	};

	////////////////  Fourier approximation  /////////////////////////////////////////////////
	////
	{
		ZGeom::VecNd xCoord(vertCount, 0), yCoord(vertCount, 0), zCoord(vertCount, 0);
		std::vector<double> xCoeff(nEig), yCoeff(nEig), zCoeff(nEig);
		ZGeom::FunctionApproximation vPursuitX;

		for (int i = 0; i < 50; ++i) {
			const ZGeom::VecNd& eigVec = mhb.getEigVec(i);
			xCoeff[i] = innerProdDiagW(vx, eigVec);
			yCoeff[i] = innerProdDiagW(vy, eigVec);
			zCoeff[i] = innerProdDiagW(vz, eigVec);
			xCoord += xCoeff[i] * eigVec;
			yCoord += yCoeff[i] * eigVec;
			zCoord += zCoeff[i] * eigVec;

			vPursuitX.addItem((vx - xCoord).norm2(), i, xCoeff[i]); 
		}
		std::ofstream ofs1("output/fourier_approx.txt");
		for (auto t : vPursuitX.getApproxItems()) {
			ofs1 << t.res() << '\t' << t.index() << '\t' << t.coeff() << std::endl;
		}
		mCoords[1] = MeshCoordinates(vertCount, &xCoord[0], &yCoord[0], &zCoord[0]);
		std::cout << "Reconstruct error (1): " << oldCoord.difference(mCoords[1]) << "\n\n";
	}	
	//////////////////////////////////////////////////////////////////////////
	
	////////////////// Fourier matching pursuit //////////////////////////////
	////
	{
		ZGeom::FunctionApproximation vPursuitX, vPursuitY, vPursuitZ;
		std::vector<ZGeom::VecNd> vBasis;
		for (int i = 0; i < nEig; ++i) vBasis.push_back(mhb.getEigVec(i));
#if 1
		std::vector<ZGeom::VecNd> vSignals;
		std::vector<ZGeom::FunctionApproximation*> vPursuits;
		vSignals.push_back(vx); vSignals.push_back(vy); vSignals.push_back(vz);
		vPursuits.push_back(&vPursuitX); vPursuits.push_back(&vPursuitY); vPursuits.push_back(&vPursuitZ);
		ZGeom::GeneralizedSimultaneousMP(vSignals, vBasis, nEig, vPursuits, innerProdDiagW, 2.);
#else
		CStopWatch timer;
		timer.startTimer();
		ZGeom::GeneralizedMP(vx, vBasis, nEig, vPursuitX, innerProdDiagW);
		timer.stopTimer("Time to compute Fourier MP: ");
		ZGeom::GeneralizedMP(vy, vBasis, nEig, vPursuitY, innerProdDiagW);
		ZGeom::GeneralizedMP(vz, vBasis, nEig, vPursuitZ, innerProdDiagW);

		std::ofstream ofs2("output/fourier_mp_approx.txt");
		for (auto t : vPursuitX.getApproxItems()) {
			ofs2 << t.res() << '\t' << t.index() << '\t' << t.coeff() << std::endl;
		}
#endif		
		mCoords[2].resize(vertCount);
		for (int i = 0; i < 50; ++i) {
			mCoords[2].getXCoord() += vPursuitX[i].coeff() * vBasis[vPursuitX[i].index()];
			mCoords[2].getYCoord() += vPursuitY[i].coeff() * vBasis[vPursuitY[i].index()];
			mCoords[2].getZCoord() += vPursuitZ[i].coeff() * vBasis[vPursuitZ[i].index()];
		}

		std::cout << "Reconstruct error (2): " << oldCoord.difference(mCoords[2]) << "\n\n";
	}
	//////////////////////////////////////////////////////////////////////////
	
#if 1
	//////////////////////////////////////////////////////////////////////////
	//// Wavelet matching pursuit
	{
		bool loadCache = (g_configMgr.getConfigValueInt("LOAD_SGW_CACHE") == 1);

		std::cout << "To compute wavelet matching pursuit" << std::endl;
		ZGeom::DenseMatrixd& matSGW = mProcessor->getWaveletMat();
		if (matSGW.empty()) {
			std::string sgwFile = "cache/" + mMesh->getMeshName() + ".sgw";
			if (loadCache && ZUtil::fileExist(sgwFile)) {
				matSGW.read(sgwFile);
			} else
			{
				mProcessor->computeSGW1(lapType);
				matSGW.write(sgwFile);
			}			
		}

		ZGeom::InnerProdcutFunc& innerProdSelected = /*innerProdDiagW;*/innerProdRegular;
		std::vector<ZGeom::VecNd>& vBasis = mEditBasis;
		
		for (int i = 0; i < matSGW.rowCount(); ++i) {
			ZGeom::VecNd newBasis = matSGW.getRowVec(i);
			newBasis.normalize(innerProdSelected);
			vBasis.push_back(newBasis);
		}

		ZGeom::FunctionApproximation vPursuitX, vPursuitY, vPursuitZ;
		ZGeom::FunctionApproximation *vPursuits[3] = {&vPursuitX, &vPursuitY, &vPursuitZ};

		const std::string waveltPursuitFile = "cache/" + mMesh->getMeshName() + ".omp";
		if (loadCache && readPursuits(waveltPursuitFile, vPursuits)) ;
		else 
		{
			std::vector<ZGeom::VecNd> vSignals;
			std::vector<ZGeom::FunctionApproximation*> vApprox;
			vSignals.push_back(vx); vSignals.push_back(vy); vSignals.push_back(vz);
			vApprox.push_back(&vPursuitX); vApprox.push_back(&vPursuitY); vApprox.push_back(&vPursuitZ);
#if 1
			CStopWatch timer;
			timer.startTimer();
			//ZGeom::SimultaneousMP(vSignals, vBasis, 100, vApprox, 2);
			//ZGeom::GeneralizedSimultaneousMP(vSignals, vBasis, 100, vApprox, innerProdSelected, 2);
			ZGeom::SimultaneousOMP(vSignals, vBasis, 100, vApprox, 2);
			//ZGeom::GeneralizedSimultaneousOMP(vSignals, vBasis, 100, vApprox, innerProdSelected, 2);
			timer.stopTimer("Time to compute Wavelet SOMP: ");
#else
			CStopWatch timer;
			timer.startTimer();
			ZGeom::OMP(vx, vBasis, 100, vPursuitX);
			//ZGeom::GeneralizedOMP(vx, vBasis, 100, vPursuitX, innerProdSelected);
			//ZGeom::MatchingPursuit(vx, vBasis, 100, vPursuitX);
			timer.stopTimer("Time to compute Wavelet OMP: ");
			ZGeom::OMP(vy, vBasis, 100, vPursuitY);
			ZGeom::OMP(vz, vBasis, 100, vPursuitZ);
			//ZGeom::GeneralizedOMP(vy, vBasis, 100, vPursuitY, innerProdSelected);
			//ZGeom::GeneralizedOMP(vz, vBasis, 100, vPursuitZ, innerProdSelected);
			//ZGeom::MatchingPursuit(vy, vBasis, 100, vPursuitY);
			//ZGeom::MatchingPursuit(vz, vBasis, 100, vPursuitZ);
#endif
			std::ofstream ofs3("output/wavelet_mp_approx.txt");
			for (auto t : vPursuitX.getApproxItems()) {
				ofs3 << t.res() << '\t' << t.index() << '\t' << t.coeff() << std::endl;
			}

			writePursuits(waveltPursuitFile, vPursuits);
		}

		//// save reconstruction result
		//
		const int countReconstructAtoms = std::min<int>(50, (int)vPursuitX.size());
		mContReconstructCoords[2].resize(countReconstructAtoms);
		mCoords[3].resize(vertCount);

		for (int i = 0; i < countReconstructAtoms; ++i) {
			mCoords[3].getXCoord() += vPursuitX[i].coeff() * vBasis[vPursuitX[i].index()];
			mCoords[3].getYCoord() += vPursuitY[i].coeff() * vBasis[vPursuitY[i].index()];
			mCoords[3].getZCoord() += vPursuitZ[i].coeff() * vBasis[vPursuitZ[i].index()];
		
			mContReconstructCoords[2][i] = mCoords[3];
		}
		std::cout << "Reconstruct error (3): " << oldCoord.difference(mCoords[3]) << "\n\n";

		mApproxCoeff[2] = vPursuitX;
	} //end of wavelet OMP
	//////////////////////////////////////////////////////////////////////////

#endif	

	std::cout << "Finish reconstructionTest1" << std::endl;
}

void ShapeEditor::reconstructionTest2( bool doWavelet /*= true*/ )
{
	std::cout << "\nStart reconstruction test2\n";

	using ZGeom::VecNd;
	CStopWatch timer;

	/* prepare prerequisite data */
	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);

	const ManifoldHarmonics &graphMHB = mProcessor->getMHB(Umbrella);
	const ManifoldHarmonics &cotMHB = mProcessor->getMHB(CotFormula);
	const ManifoldHarmonics &cotSymMHB = mProcessor->getMHB(SymCot);

	std::vector<double> wDiag;
	mProcessor->getMeshLaplacian(CotFormula).getW().getDiagonal(wDiag);
	ZGeom::InnerProdcutFunc innerProdDiagW = [&](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2) {
		double *y = new double[vertCount];
		vdmul(&vertCount, v1.c_ptr(), &wDiag[0], &y[0]); 
		//for (int i = 0; i < vertCount; ++i) y[i] = v1[i] * wDiag[i];
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

		for (int i = 0; i < nEigTotal; ++i) vAtoms.push_back(cotMHB.getEigVec(i));

		GeneralizedSimultaneousFourierApprox(vSignals, vAtoms, nAtomSel, vApproxCoeff, innerProdDiagW);
		
		computeApproximations(vAtoms, &vApproxCoeff[0], nReconstruct, mContReconstructCoords[0], mCoords[1]);
		
		std::cout << "Reconstruct error (1): " << oldCoord.difference(mCoords[1]) << "\n\n";
	}	
	//////////////////////////////////////////////////////////////////////////


	/* Test 2, Simultaneous Fourier Matching Pursuit */
	{
		for (auto p : vApproxCoeff) p->clear();
		vAtoms.clear();
		
		for (int i = 0; i < nEigTotal; ++i) vAtoms.push_back(cotMHB.getEigVec(i));

		timer.startTimer();
		ZGeom::GeneralizedSimultaneousMP(vSignals, vAtoms, nAtomSel, vApproxCoeff, innerProdDiagW, 2.);
		timer.stopTimer("Time to compute Fourier SMP: ");

		computeApproximations(vAtoms, &vApproxCoeff[0], nReconstruct, mContReconstructCoords[1], mCoords[2]);
		
		std::vector<int> vSelectedAtomIdx = vApproxX.getAllAtomIndex();	
		int selectedFromTop = std::count_if(vSelectedAtomIdx.begin(), vSelectedAtomIdx.end(), [&](int idx) { return idx < vSelectedAtomIdx.size();} );
		std::cout << "Ratio of MP overlapping: " << (double)selectedFromTop / (double)vSelectedAtomIdx.size() << '\n';
		std::cout << "Reconstruct error (2): " << oldCoord.difference(mCoords[2]) << "\n\n";
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
	
	std::cout << "Finish reconstructionTest2" << std::endl;
}

/////// compute various eigenvectors indexed by Fiedler vector ////////////////
////
void ShapeEditor::editTest2()	
{
	const int vertCount = mMesh->vertCount();
	LaplacianType lapType = SymCot;
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW, ZGeom::MAT_UPPER);
	const ManifoldHarmonics &mhb = mProcessor->getMHB(lapType);
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

