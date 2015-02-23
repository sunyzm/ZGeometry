#ifndef ZMESH_MESH_HOLE_H
#define ZMESH_MESH_HOLE_H

#include <vector>
#include "Mesh.h"

class MeshHole
{
public:
    std::vector<int> mHoleFaces;
    std::vector<int> mHoleVerts;
    std::vector<int> mHoleBoundaryVerts;
};

MeshHole autoGenerateHole(const CMesh& mesh, int seedVert, int holeSize);

MeshHole autoGenerateHole(const CMesh& mesh, const std::vector<int>& seedVerts, int totalSize);

#endif