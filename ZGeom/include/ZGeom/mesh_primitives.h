#ifndef MESH_REGION_H
#define MESH_REGION_H
#include <vector>
#include "Mesh.h"

namespace ZGeom 
{

struct MeshRegion
{
    std::vector<int> vert_on_boundary;
    std::vector<int> he_on_boundary;
    std::vector<int> vert_inside;
    std::vector<int> face_inside;

    bool is_outer_boundary;
    double adjacent_edge_length;

    // constructors
    MeshRegion() : is_outer_boundary(false), adjacent_edge_length(-1) { }
    MeshRegion(const MeshRegion&) = default;
    MeshRegion& operator = (MeshRegion&& hb);
    MeshRegion(MeshRegion&& hb) { *this = std::move(hb); }
    const std::vector<int>& getInsideFaceIdx() const { return face_inside; }

    void determineBoundaryHalfEdges(const CMesh& mesh);
};

std::vector<int> vertSurroundingVerts(const CMesh& mesh, const std::vector<int>& vert_inside, int ring);
std::set<int> meshMultiVertsAdjacentVerts(const CMesh& mesh, const std::vector<int>& vert, int ring, bool inclusive = true);
std::vector<int> getFaceEncompassedByVerts(const CMesh& mesh, const std::vector<int>& verts);

MeshRegion meshRegionFromVerts(CMesh& mesh, const std::vector<int>& inside_verts);
MeshRegion meshRegionFromDistField(CMesh& mesh, const std::vector<double>& dist_field, int seed, std::function<bool(double)> judge_in_region);


}

#endif