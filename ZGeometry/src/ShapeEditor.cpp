#include "ShapeEditor.h"
#include <iostream>
#include <random>
#include <functional>
#include <ppl.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZUtil/timer.h>
#include <ZUtil/zassert.h>
#include <ZUtil/zutil_io.h>
#include "global.h"


Vector3D toVector3D(const ZGeom::Vec3d& v)
{
	return Vector3D(v[0], v[1], v[2]);
}

ZGeom::Vec3d toVec3d(const Vector3D& v)
{
	return ZGeom::Vec3d(v.x, v.y, v.z);
}

void ShapeEditor::init( DifferentialMeshProcessor* processor )
{
	mProcessor = processor; 
	mMesh = processor->getMesh();
	mEngine = processor->getMatlabEngineWrapper();
	
	mCoordSelect = 0;
	mCoords.resize(4);
	processor->getMesh()->getVertCoordinates(mCoords[0]); 

	std::cout << "Shape editor is initialized!" << std::endl;

	//addNoise(0.1);
	//processor->getMesh()->getVertCoordinates(mOldCoord); 
	//editTest1();
// 	ZGeom::DenseMatrixd hkMat1, hkMat2;
// 	processor->computeHeatKernelMat(30, hkMat1);
// 	processor->computeHeatKernelMat_AMP(30, hkMat2);
// 
// 	double maxError(0);
// 	for (int i = 0; i < hkMat1.rowCount(); ++i) {
// 		for (int j = 0; j < hkMat1.colCount(); ++j)
// 			if (std::fabs(hkMat1(i,j) - hkMat2(i,j)) > maxError) maxError = std::fabs(hkMat1(i,j) - hkMat2(i,j));
// 	}
// 
// 	std::cout << "max hkmat error: " << maxError << std::endl;

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

void ShapeEditor::manifoldHarmonicsReconstruct( int nEig )
{
	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	MeshCoordinates newCoord(vertCount);
	MeshLaplacian::LaplacianType lapType = MeshLaplacian::SymCot;
		
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

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(MeshLaplacian::SymCot).getLS();
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

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(MeshLaplacian::SymCot).getLS();
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

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(MeshLaplacian::SymCot).getLS();
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
	if (matSGW.empty()) mProcessor->computeSGW(MeshLaplacian::SymCot);
	const int waveletCount = matSGW.rowCount();
	
	/*std::vector<double> diagW;
	mProcessor->getMeshLaplacian(MeshLaplacian::CotFormula).getW().getDiagonal(diagW);
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
		mProcessor->computeSGW();
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
	const MeshLaplacian& laplacian = mProcessor->getMeshLaplacian(MeshLaplacian::CotFormula);
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

bool readPursuits(const std::string& pursuitFile, std::vector<ZGeom::PursuitApproxItem> *vPursuits[]) {
	if (!ZUtil::fileExist(pursuitFile)) 
		return false;

	std::ifstream ifs(pursuitFile.c_str());
	int atomCount;
	ifs >> atomCount;
	
	for (int b = 0; b < 3; ++b) {
		vPursuits[b]->resize(atomCount);
		for (int i = 0; i < atomCount; ++i) {
			ZGeom::PursuitApproxItem& item = (*vPursuits[b])[i];
			ifs >> std::get<0>(item) >> std::get<1>(item) >> std::get<2>(item);
		}
	}

	return true;
}

void writePursuits(const std::string& pursuitFile, std::vector<ZGeom::PursuitApproxItem> *vPursuits[]) {
	std::ofstream ofs(pursuitFile.c_str());
	int atomCount = vPursuits[0]->size();
	ofs << atomCount << std::endl;
	for (int b = 0; b < 3; ++b) {
		for (int i = 0; i < atomCount; ++i) {
			const ZGeom::PursuitApproxItem& item = (*vPursuits[b])[i];
			ofs << std::get<0>(item) << ' ' << std::get<1>(item) << ' ' << std::get<2>(item) << std::endl;
		}
	}	
}

void ShapeEditor::editTest1()
{
	const int vertCount = mMesh->vertCount();
	MeshLaplacian::LaplacianType lapType = MeshLaplacian::CotFormula;
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
		
	ZGeom::InnerProdcutFunc innerProdDiagW = [&](const ZGeom::VecNd& v1, const ZGeom::VecNd& v2)->double {
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
	ZGeom::VecNd xCoord(vertCount, 0), yCoord(vertCount, 0), zCoord(vertCount, 0);
	std::vector<double> xCoeff(nEig), yCoeff(nEig), zCoeff(nEig);
	std::vector<ZGeom::PursuitApproxItem> vPursuitX;

	for (int i = 0; i < 50; ++i) {
		const ZGeom::VecNd& eigVec = mhb.getEigVec(i);
		xCoeff[i] = innerProdDiagW(vx, eigVec);
		yCoeff[i] = innerProdDiagW(vy, eigVec);
		zCoeff[i] = innerProdDiagW(vz, eigVec);
		xCoord += xCoeff[i] * eigVec;
		yCoord += yCoeff[i] * eigVec;
		zCoord += zCoeff[i] * eigVec;

		vPursuitX.push_back(std::make_tuple((xCoord-vx).norm2(), i, xCoeff[i])); 
	}
	std::ofstream ofs1("output/fourier_approx.txt");
	for (auto t : vPursuitX) {
		ofs1 << std::get<0>(t) << '\t' << std::get<1>(t) << '\t' << std::get<2>(t) << std::endl;
	}
	mCoords[1] = MeshCoordinates(vertCount, &xCoord[0], &yCoord[0], &zCoord[0]);
	//////////////////////////////////////////////////////////////////////////
	
	////////////////// Fourier matching pursuit //////////////////////////////
	////
	{
		std::vector<ZGeom::PursuitApproxItem> vPursuitX, vPursuitY, vPursuitZ;
		std::vector<ZGeom::VecNd> vBasis;
		for (int i = 0; i < nEig; ++i) vBasis.push_back(mhb.getEigVec(i));

		CStopWatch timer;
		timer.startTimer();
		ZGeom::MatchingPursuit(vx, vBasis, innerProdDiagW, nEig, vPursuitX);
		timer.stopTimer("Time to compute Fourier MP: ");

		std::ofstream ofs2("output/fourier_mp_approx.txt");
		for (auto t : vPursuitX) {
			ofs2 << std::get<0>(t) << '\t' << std::get<1>(t) << '\t' << std::get<2>(t) << std::endl;
		}

		ZGeom::MatchingPursuit(vy, vBasis, innerProdDiagW, nEig, vPursuitY);
		ZGeom::MatchingPursuit(vz, vBasis, innerProdDiagW, nEig, vPursuitZ);
		
		mCoords[2].resize(vertCount);
		for (int i = 0; i < 50; ++i) {
			mCoords[2].getXCoord() += std::get<2>(vPursuitX[i]) * vBasis[std::get<1>(vPursuitX[i])];
			mCoords[2].getYCoord() += std::get<2>(vPursuitY[i]) * vBasis[std::get<1>(vPursuitY[i])];
			mCoords[2].getZCoord() += std::get<2>(vPursuitZ[i]) * vBasis[std::get<1>(vPursuitZ[i])];
		}

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
				mProcessor->computeSGW(lapType);
				matSGW.write(sgwFile);
			}			
		}

		ZGeom::InnerProdcutFunc& innerProdSelected = /*innerProdDiagW;*/innerProdRegular;
		std::vector<ZGeom::VecNd>& vBasis = mEditBasis;
		
		//for (int i = 0; i < nEig; ++i) vBasis.push_back(mhb.getEigVec(i));
		for (int i = 0; i < matSGW.rowCount(); ++i) {
			ZGeom::VecNd newBasis = matSGW.getRowVec(i);
			newBasis.normalize(innerProdSelected);
			vBasis.push_back(newBasis);
		}

		std::vector<ZGeom::PursuitApproxItem> vPursuitX, vPursuitY, vPursuitZ;
		std::vector<ZGeom::PursuitApproxItem> *vPursuits[3] = {&vPursuitX, &vPursuitY, &vPursuitZ};

		const std::string waveltPursuitFile = "cache/" + mMesh->getMeshName() + ".omp";
		if (loadCache && readPursuits(waveltPursuitFile, vPursuits)) ;
		else 
		{
			CStopWatch timer;
			timer.startTimer();
			//ZGeom::MatchingPursuit(vx, vBasis, innerProdSelected, 100, vPursuit);
			ZGeom::OrthogonalMatchingPursuit(vx, vBasis, innerProdSelected, 100, vPursuitX);
			//ZGeom::OrthogonalMatchingPursuit(vx, vBasis, innerProdSelected, 100, vPursuitX, *mEngine);
			timer.stopTimer("Time to compute Wavelet MP: ");

			std::ofstream ofs3("output/wavelet_mp_approx.txt");
			for (auto t : vPursuitX) {
				ofs3 << std::get<0>(t) << '\t' << std::get<1>(t) << '\t' << std::get<2>(t) << std::endl;
			}

			ZGeom::OrthogonalMatchingPursuit(vy, vBasis, innerProdSelected, 100, vPursuitY);
			ZGeom::OrthogonalMatchingPursuit(vz, vBasis, innerProdSelected, 100, vPursuitZ);
// 
			writePursuits(waveltPursuitFile, vPursuits);
		}

		mCoords[3].resize(vertCount);
		for (int i = 0; i < 50; ++i) {
			mCoords[3].getXCoord() += std::get<2>(vPursuitX[i]) * vBasis[std::get<1>(vPursuitX[i])];
			mCoords[3].getYCoord() += std::get<2>(vPursuitY[i]) * vBasis[std::get<1>(vPursuitY[i])];
			mCoords[3].getZCoord() += std::get<2>(vPursuitZ[i]) * vBasis[std::get<1>(vPursuitZ[i])];
		}
	} //end of wavelet OMP
	//////////////////////////////////////////////////////////////////////////

	mCoordSelect = 3;
	mMesh->setVertCoordinates(mCoords[mCoordSelect]);
#endif	

	std::cout << "Finish editTest1" << std::endl;
}

/////// compute various eigenvectors indexed by Fiedler vector ////////////////
////
void ShapeEditor::editTest2()	
{
	const int vertCount = mMesh->vertCount();
	MeshLaplacian::LaplacianType lapType = MeshLaplacian::SymCot;
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

void ShapeEditor::deformLaplacian2()
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

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(MeshLaplacian::SymCot).getLS();

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
