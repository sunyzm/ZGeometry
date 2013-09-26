#include "ShapeEditor.h"
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

void ShapeEditor::manifoldHarmonicsReconstruct( int nEig )
{
	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	MeshCoordinates newCoord(vertCount);
	MeshLaplacian::LaplacianType lapType = MeshLaplacian::SymCot;
		
	ZGeom::SparseMatrixCSR<double, int> matW;
	mProcessor->getMeshLaplacian(lapType).getW().convertToCSR(matW);
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

void ShapeEditor::differentialDeform()
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
	mEngine->addVariable(oldCoord.getXCoord().c_ptr(), vertCount, 1, false, "ecx");
	mEngine->addVariable(oldCoord.getYCoord().c_ptr(), vertCount, 1, false, "ecy");
	mEngine->addVariable(oldCoord.getZCoord().c_ptr(), vertCount, 1, false, "ecz");

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
	std::vector<int> rowInd, colInd;
	std::vector<double> vals;
	matOptS.copyElements(matLs);
	for (int a = 0; a < anchorCount; ++a) 
		matOptS.insertElem(vertCount + a + 1, anchorIndex[a] + 1, anchorWeight);
	matOptS.convertToCOO(rowInd, colInd, vals, ZGeom::MAT_FULL);
	mEngine->addSparseMat(&rowInd[0], &colInd[0], &vals[0], vertCount + anchorCount, vertCount, matOptS.nonzeroCount(), "matOptS");
	
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

void ShapeEditor::spectralWaveletDeform()
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

	MeshCoordinates newCoord = oldCoord;
	mMesh->setVertCoordinates(newCoord);
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

	std::vector<double> vDist2Handle(vertCount);
	concurrency::parallel_for (0, vertCount, [&](int i) {
		double dist = mMesh->calGeodesic(hIdx, i);
		vDist2Handle[i] = dist;
	});
	double distMax = *std::max_element(vDist2Handle.begin(), vDist2Handle.end());

	std::vector<Vector3D> vDeformedPos(vertCount);
	for (int i = 0; i < vertCount; ++i)
	{
		vDeformedPos[i] = mMesh->getVertex(i)->getPosition() + handleTrans * (1.0 - vDist2Handle[i]/distMax);
	}

	MeshCoordinates newCoord(vertCount);
	ZGeom::VecNd& xNewCoord = newCoord.getXCoord();
	ZGeom::VecNd& yNewCoord = newCoord.getYCoord();
	ZGeom::VecNd& zNewCoord = newCoord.getZCoord();
	
	for (int i = 0; i < vertCount; ++i) {
		xNewCoord[i] = vDeformedPos[i].x;
		yNewCoord[i] = vDeformedPos[i].y;
		zNewCoord[i] = vDeformedPos[i].z;
	}
	mMesh->setVertCoordinates(newCoord);

	timer.stopTimer("Deformation time: ");
}
