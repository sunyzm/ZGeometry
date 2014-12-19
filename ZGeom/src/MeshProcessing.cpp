#include "MeshProcessing.h"
#include <cstdlib>
#include <random>
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
	const double *pEigVals = &(hb.getAllEigVals()[0]);
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
	const double *pEigVals = &(hb.getAllEigVals()[0]);
	double *pEigVec = matEigVecs.raw_ptr();
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(hb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);

	std::vector<double> vDiag(eigCount);
	for (int i = 0; i < eigCount; ++i) 
		vDiag[i] = gen(hb.getEigVal(i), t);

	hk.resize(vertCount, vertCount);

	quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], hk.raw_ptr());
}

std::vector<int> randomHoleVertex(const CMesh& mesh, int hole_size, int seed /*= -1*/)
{
    using std::vector;
    assert(hole_size >= 1);
    int N = mesh.vertCount();
    std::default_random_engine generator;
    if (seed == -1) {
        std::uniform_int_distribution<int> distribution(0, N-1);
        seed = distribution(generator);
    }

    std::set<int> vertInHole;
    std::vector<int> vertFeasible;

    vertInHole.insert(seed);
    vertFeasible.push_back(seed)
        ;
    while (vertInHole.size() < hole_size && !vertFeasible.empty()) {
        std::uniform_int_distribution<int> distr1(0, (int)vertFeasible.size() - 1);
        int selPos = distr1(generator);
        int selIdx = vertFeasible[selPos];
        vector<int> neighbors = mesh.getVertNeighborVerts(selIdx, 1, false);
        vector<int> feasibleNeighbors;
        for (int vIdx : neighbors) {
            if (vertInHole.find(vIdx) == vertInHole.end())
                feasibleNeighbors.push_back(vIdx);
        }
        if (feasibleNeighbors.empty()) {
            vertFeasible.erase(vertFeasible.begin() + selPos);
        }
        else {
            std::uniform_int_distribution<int> distr2(0, (int)feasibleNeighbors.size() - 1);
            int selNeighbor = feasibleNeighbors[distr2(generator)];
            vertInHole.insert(selNeighbor);
            if (feasibleNeighbors.size() == 1) vertFeasible.erase(vertFeasible.begin() + selPos);
            vertFeasible.push_back(selNeighbor);        
        }
    }

    return vector<int>(vertInHole.begin(), vertInHole.end());
}

} // end of namespace