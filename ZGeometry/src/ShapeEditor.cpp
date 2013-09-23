#include "ShapeEditor.h"
#include <ZGeom/SparseMatrix.h>
#include <ZGeom/SparseMatrixCSR.h>
#include <ZGeom/VecN.h>
#include <ZGeom/MatVecArithmetic.h>

void ShapeEditor::manifoldHarmonicsReconstruct( int nEig )
{
	const int vertCount = mMesh->vertCount();
	MeshCoordinates oldCoord, newCoord;
	newCoord.resize(vertCount);
	MeshLaplacian::LaplacianType lapType = MeshLaplacian::SymCot;

	mMesh->getVertCoordinates(oldCoord);	
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
