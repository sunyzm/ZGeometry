#ifndef ZGEOM_MESH_PROCESSING_H
#define ZGEOM_MESH_PROCESSING_H
#include <functional>
#include "Mesh.h"
#include "SparseMatrix.h"
#include "EigenSystem.h"
#include "DenseMatrix.h"

namespace ZGeom
{
	enum MeshMatrixType { 
		MM_ADJACENCY, MM_DEGREE, MM_INV_DEGREE, MM_WALK, 
		MM_GRAPH_LAPLACE, MM_NORMALIZED_GRAPH_LAPLACE
	};

	void ConstructMeshMatrix(const CMesh& mesh, MeshMatrixType mmt, SparseMatrix<double>& meshMat);
	void ComputeHeatKernelMatrix(const EigenSystem& hb, double t, DenseMatrixd& hk);
	void ComputeKernelMatrix(const EigenSystem& hb, double t, std::function<double(double,double)> gen, DenseMatrixd& hk);
}

#endif