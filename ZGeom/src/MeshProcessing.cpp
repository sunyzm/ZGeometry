#include "MeshProcessing.h"
#include "MatVecArithmetic.h"

namespace ZGeom {

void ConstructMeshMatrix( const CMesh& mesh, MeshMatrixType mmt, SparseMatrix<double>& meshMat )
{
	const int vertCount = mesh.vertCount();	

	if (mmt == MM_DEGREE) {
		std::vector<int> vDiagDegree(vertCount);
		for (int i = 0; i < vertCount; ++i) {
			vDiagDegree[i] = (int)mesh.getVertNeighborVerts(i, 1, false).size();
			assert(vDiagDegree[i] > 0);
		}
		meshMat.convertFromDiagonal(vDiagDegree);
	} 
	else if (mmt == MM_INV_DEGREE) {
		std::vector<double> vDiagDegree(vertCount);
		for (int i = 0; i < vertCount; ++i) {
			vDiagDegree[i] = 1.0 / (int)mesh.getVertNeighborVerts(i, 1, false).size();
		}
		meshMat.convertFromDiagonal(vDiagDegree);
	}
	else if (mmt == MM_ADJACENCY) {
		std::vector<std::tuple<int,int,double> > vElem;
		for (int i = 0; i < vertCount; ++i) {
			std::vector<int> vNeighbors = mesh.getVertNeighborVerts(i, 1, false);
			int valence = (int)vNeighbors.size();
			for (int j = 0; j < valence; ++j)
				vElem.push_back(std::make_tuple(i+1, vNeighbors[j]+1, 1.));
		}
		meshMat.convertFromCOO(vertCount, vertCount, vElem);
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
		ConstructMeshMatrix(mesh, MM_GRAPH_LAPLACE, meshMat);

		std::vector<double> vInvSqD;	// D^(-0.5)
		matD.getDiagonal(vInvSqD);
		for (double& v : vInvSqD) v = std::pow(v, -0.5);

		for (auto& elem : meshMat.allElements()) {
			elem.val() *= vInvSqD[elem.row()-1] * vInvSqD[elem.col()-1];
		}
	}

}

} //end of namespace