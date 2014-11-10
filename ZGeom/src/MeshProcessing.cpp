#include "MeshProcessing.h"
#include "MatVecArithmetic.h"
#include "arithmetic.h"

namespace ZGeom {

void ConstructMeshMatrix( const CMesh& mesh, MeshMatrixType mmt, SparseMatrix<double>& meshMat )
{
    using namespace std;
	const int vertCount = mesh.vertCount();	

	if (mmt == MM_DEGREE) {
		vector<int> vDiagDegree(vertCount);
		for (int i = 0; i < vertCount; ++i) {
			vDiagDegree[i] = (int)mesh.getVertNeighborVerts(i, 1, false).size();
		}
		meshMat.convertFromDiagonal(vDiagDegree);
	} 
	else if (mmt == MM_INV_DEGREE) {
		vector<double> vDiagDegree(vertCount);
		for (int i = 0; i < vertCount; ++i) {
			vDiagDegree[i] = 1.0 / (int)mesh.getVertNeighborVerts(i, 1, false).size();
		}
		meshMat.convertFromDiagonal(vDiagDegree);
	}
	else if (mmt == MM_ADJACENCY) {
		vector<tuple<int,int,double> > vElem;
		for (int i = 0; i < vertCount; ++i) {
			vector<int> vNeighbors = mesh.getVertNeighborVerts(i, 1, false);
			int valence = (int)vNeighbors.size();
			for (int j = 0; j < valence; ++j)
				vElem.push_back(std::make_tuple(i+1, vNeighbors[j]+1, 1.));
		}
		meshMat.convertFromCOO(vertCount, vertCount, vElem);
	}
    else if (mmt == MM_INV_LENGTH) {
        
    }
	else if (mmt == MM_WALK) {
		// W = D^(-1)*A
		ZGeom::SparseMatrix<double> matInvD, matA;
		ConstructMeshMatrix(mesh, MM_INV_DEGREE, matInvD);
		ConstructMeshMatrix(mesh, MM_ADJACENCY, matA);
		mulMatMat(matInvD, matA, meshMat);
	}
	else if (mmt == MM_GRAPH_LAPLACE) {
		// L = D - A
		ZGeom::SparseMatrix<double> matD, matA;
		ConstructMeshMatrix(mesh, MM_DEGREE, matD);
		ConstructMeshMatrix(mesh, MM_ADJACENCY, matA);
		addMatMat(matD, matA, -1.0, meshMat);
	}	
	else if (mmt == MM_NORMALIZED_GRAPH_LAPLACE) {
		// Ln = D^(-0.5)*L*D^(-0.5) = I - D^(-0.5)*A*D^(-0.5)
		ZGeom::SparseMatrix<double> matD;
		ConstructMeshMatrix(mesh, MM_DEGREE, matD);
		std::vector<double> vInvSqD = matD.getDiagonal();
		for (double& v : vInvSqD) v = std::pow(v, -0.5); // D^(-0.5)

		ConstructMeshMatrix(mesh, MM_GRAPH_LAPLACE, meshMat);
		for (auto& elem : meshMat.allElements()) {
			elem.val() *= vInvSqD[elem.row()-1] * vInvSqD[elem.col()-1];
		}
	}
}

void ComputeHeatKernelMatrix( const EigenSystem& hb, double t, DenseMatrixd& hk )
{
	const int vertCount = hb.eigVecSize();
	const int eigCount = hb.eigVecCount();

	DenseMatrixd matEigVecs(eigCount, vertCount);	
	const double *pEigVals = &(hb.getEigVals()[0]);
	double *pEigVec = matEigVecs.raw_ptr();
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(hb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);
	
	std::vector<double> vDiag(eigCount);
	for (int i = 0; i < eigCount; ++i) 
		vDiag[i] = std::exp(-hb.getEigVal(i) * t);

	hk.resize(vertCount, vertCount);

	quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], hk.raw_ptr());
}

void ComputeKernelMatrix( const EigenSystem& hb, double t, std::function<double(double,double)> gen, DenseMatrixd& hk )
{
	const int vertCount = hb.eigVecSize();
	const int eigCount = hb.eigVecCount();

	DenseMatrixd matEigVecs(eigCount, vertCount);	
	const double *pEigVals = &(hb.getEigVals()[0]);
	double *pEigVec = matEigVecs.raw_ptr();
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(hb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);

	std::vector<double> vDiag(eigCount);
	for (int i = 0; i < eigCount; ++i) 
		vDiag[i] = gen(hb.getEigVal(i), t);

	hk.resize(vertCount, vertCount);

	quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], hk.raw_ptr());
}

} //end of namespace