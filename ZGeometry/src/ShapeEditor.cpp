#include "ShapeEditor.h"
#include <ZGeom/SparseMatrix.h>
#include <ZGeom/SparseMatrixCSR.h>
#include <ZGeom/VecN.h>
#include <ZGeom/MatVecArithmetic.h>

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
	assert(nEig <= mhb.eigVecCount());

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
	const std::map<int, Vector3D>& anchors = mProcessor->getHandles();
	const int anchorCount = anchors.size();
	if (anchors.size() == 0) {
		std::cout << "At least one anchor need to be picked!";
		return;
	}

	const int anchorWeight = 1.0;

	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord;
	mMesh->getVertCoordinates(oldCoord);
	mEngine->addVariable(oldCoord.getCoordFunc(0).c_ptr(), vertCount, 1, false, "ecx");
	mEngine->addVariable(oldCoord.getCoordFunc(1).c_ptr(), vertCount, 1, false, "ecy");
	mEngine->addVariable(oldCoord.getCoordFunc(2).c_ptr(), vertCount, 1, false, "ecz");

	std::vector<int> anchorIndex;
	std::vector<Vector3D> anchorPos;
	for (auto a : anchors) {
		anchorIndex.push_back(a.first);
		anchorPos.push_back(a.second);
	}

	const ZGeom::SparseMatrix<double>& matLs = mProcessor->getMeshLaplacian(MeshLaplacian::SymCot).getLS();
	ZGeom::SparseMatVecMultiplier mulLs(matLs, true);
	
	ZGeom::DenseMatrixd denseMat(vertCount, vertCount);
	matLs.convertToFull<double>(denseMat.raw_ptr(), ZGeom::MAT_FULL);
	mEngine->addVariable(denseMat.raw_ptr(), vertCount, vertCount, true, "denseLs");
	ZGeom::DenseMatVecMultiplier denseMulLs(denseMat);
	
	ZGeom::VecNd diffCoord[3];
	for (int i = 0; i < 3; ++i) {
		diffCoord[i].resize(vertCount);
		denseMulLs.mul(oldCoord.getCoordFunc(i), diffCoord[i]);		
// 		diffCoord[i].resize(vertCount);
// 		mulLs.mul(oldCoord.getCoordFunc(i), diffCoord[i]);
	}

	ZGeom::VecNd solveRHS[3];
	for (int i = 0; i < 3; ++i ) {
		solveRHS[i].resize(vertCount + anchorCount, 0);
		std::copy_n(diffCoord[i].c_ptr(), vertCount, solveRHS[i].c_ptr());
		for (int l = 0; l < anchorCount; ++l) {
			solveRHS[i][vertCount + l] = anchorWeight * anchorPos[l][i];
		}
	}

	ZGeom::DenseMatrixd matOpt(vertCount + anchorCount, vertCount);
	matOpt.copyRows(denseMat, 0);
	double *matOptData = matOpt.raw_ptr();
	for (int i = 0; i < anchorCount; ++i) {
		for (int j = 0; j < vertCount; ++j) matOptData[(vertCount+i)*vertCount+j] = 0.0;
		matOptData[(vertCount+i)*vertCount + anchorIndex[i]] = anchorWeight;
	}

	mEngine->addVariable(matOptData, vertCount + anchorCount, vertCount, true, "matOpt");
	mEngine->addColVecVariable(solveRHS[0].c_ptr(), vertCount + anchorCount, "dcx");
	mEngine->addColVecVariable(solveRHS[1].c_ptr(), vertCount + anchorCount, "dcy");
	mEngine->addColVecVariable(solveRHS[2].c_ptr(), vertCount + anchorCount, "dcz");

	mEngine->eval("lsx=matOpt\\dcx;");
	mEngine->eval("lsy=matOpt\\dcy;");
	mEngine->eval("lsz=matOpt\\dcz;");

	double *lsx = mEngine->getVariablePtr("lsx");
	double *lsy = mEngine->getVariablePtr("lsy");
	double *lsz = mEngine->getVariablePtr("lsz");

	MeshCoordinates newCoord(vertCount);
	std::copy_n(lsx, vertCount, newCoord.getCoordFunc(0).c_ptr());
	std::copy_n(lsy, vertCount, newCoord.getCoordFunc(1).c_ptr());
	std::copy_n(lsz, vertCount, newCoord.getCoordFunc(2).c_ptr());

	mMesh->setVertCoordinates(newCoord);
}
