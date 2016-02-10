#include "mesh_partition.h"
#include <metis.h>
#include "SparseMatrix.h"
#include "geometry_processing.h"


std::vector<int> ZGeom::MetisMeshPartition(const CMesh* mesh, int nPart)
{
    int vertCount = mesh->vertCount();
    std::vector<int> vPart(vertCount);
    int ncon = 1;
    std::vector<int> xadj, adjncy;
    getMeshGraphCSR(*mesh, xadj, adjncy);
    int objval;
    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_NUMBERING] = 0;

    METIS_PartGraphKway(&vertCount, &ncon, &xadj[0], &adjncy[0], NULL, NULL, NULL, &nPart, NULL, NULL, NULL, &objval, &vPart[0]);

    return vPart;
}