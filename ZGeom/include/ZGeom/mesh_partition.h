#include <vector>
#include "Mesh.h"

namespace ZGeom {
    std::vector<int> MetisMeshPartition(const CMesh* mesh, int nPart);
    std::vector<int> FiedlerMeshPartition(CMesh* original_mesh, int nParts);
}

