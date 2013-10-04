#include "ShapeEditor.h"
#include <iostream>
#include <ppl.h>
#include <ZGeom/SparseMatrix.h>
#include <ZGeom/SparseMatrixCSR.h>
#include <ZGeom/VecN.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZUtil/timer.h>
#include <ZUtil/zassert.h>

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
	processor->getMesh()->getVertCoordinates(mOldCoord); 
	mEngine = processor->getMatlabEngineWrapper();
	std::cout << "Shape editor is initialized!" << std::endl;
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
	ZGeom::addMatMat(matL1, matBiL, kb, matMixedL);
	mEngine->addSparseMat(matMixedL, "matMixedL");

	ZGeom::SparseMatVecMultiplier mulMixedL(matMixedL, true);	
	ZGeom::VecNd diffCoord[3];
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(vertCount);
		mulMixedL.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
	}

	std::set<int> freeVerts;
	for (int i = 0; i < anchorCount; ++i) {
		std::set<int> sFree;
		mMesh->vertRingNeighborVerts(anchorIndex[i], 5, sFree, true);
		for (int v : sFree) freeVerts.insert(v);
	}
	std::vector<int> fixedVerts;
	for (int i = 0; i < vertCount; ++i) {
		if (freeVerts.find(i) == freeVerts.end())
			fixedVerts.push_back(i);
	}
	const int fixedCount = fixedVerts.size();

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount + fixedCount, 0);
		solveRHS[i].copyElements(diffCoord[i], 0);
		for (int l = 0; l < anchorCount; ++l) {
			solveRHS[i][vertCount + l] = anchorWeight * anchorPos[l][i];
		}
		for (int l = 0; l < fixedCount; ++l) {
			solveRHS[i][vertCount + anchorCount + l] = mMesh->getVertexPosition(fixedVerts[l])[i];
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

	MeshCoordinates newCoord(vertCount, lsx, lsy, lsz);
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
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;
	const int vertCount = mMesh->vertCount();

	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	mEngine->addColVec(oldCoord.getXCoord(), "ecx");
	mEngine->addColVec(oldCoord.getYCoord(), "ecy");
	mEngine->addColVec(oldCoord.getZCoord(), "ecz");
	
	const ZGeom::DenseMatrixd& matW = mProcessor->getWaveletMat();

	if (matW.empty()) {
		mProcessor->computeSGW();
		mEngine->addDenseMat(matW, "matSGW");
	}

	timer.stopTimer("Prepare deformation time: ");
#if 0
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
#endif
}
