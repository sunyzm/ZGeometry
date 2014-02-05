#ifndef ZGEOM_MESH_PROCESSING_H
#define ZGEOM_MESH_PROCESSING_H
#include "Mesh.h"
#include "SparseMatrix.h"

namespace ZGeom
{
	enum MeshMatrixType {MM_ADJACENCY, MM_DEGREE, MM_INV_DEGREE, MM_WALK, MM_GRAPH_LAPLACE, MM_NORMALIZED_GRAPH_LAPLACE};

	void ConstructMeshMatrix(const CMesh& mesh, MeshMatrixType mmt, SparseMatrix<double>& meshMat);
}

#endif