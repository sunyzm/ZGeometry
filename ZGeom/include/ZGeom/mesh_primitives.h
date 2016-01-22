#ifndef MESH_REGION_H
#define MESH_REGION_H
#include <vector>
#include "Mesh.h"

namespace ZGeom {

struct MeshRegion
{
    std::vector<int> vert_inside;
    std::vector<int> vert_on_boundary;
    std::vector<int> face_inside;
    std::vector<int> he_on_boundary;
    
    double adjacent_edge_length;

    // constructors
    MeshRegion() : adjacent_edge_length(-1) { }
    MeshRegion(const MeshRegion&) = default;
    MeshRegion& operator = (MeshRegion&& hb);
    MeshRegion(MeshRegion&& hb) { *this = std::move(hb); }
    const std::vector<int>& getInsideFaceIdx() const { return face_inside; }

    void determineBoundaryHalfEdges(const CMesh& mesh);
};

struct BandedMeshRegion {
    void resize(int num)
    {
        band_verts.resize(num);
        concentric_regions.resize(num);
        band_thresholds.resize(num);
    }

    int numOfRegions() const { return (int)band_verts.size(); }
    
    std::vector<std::vector<int>> band_verts;
    std::vector<MeshRegion> concentric_regions;
    std::vector<double> band_thresholds;    //optional thresholds for each band
};

std::vector<int> vertSurroundingVerts(const CMesh& mesh, const std::vector<int>& vert_inside, int ring);
std::set<int> meshMultiVertsAdjacentVerts(const CMesh& mesh, const std::vector<int>& vert, int ring, bool inclusive = true);
std::vector<int> getFaceEncompassedByVerts(const CMesh& mesh, const std::vector<int>& verts);

MeshRegion meshRegionFromInsideVerts(const CMesh& mesh, const std::vector<int>& inside_verts);
MeshRegion meshRegionFromVerts(const CMesh& mesh, const std::vector<int>& region_verts);

MeshRegion meshRegionFromDistField(const CMesh& mesh, const std::vector<double>& dist_field, int seed, std::function<bool(double)> judge_in_region);

BandedMeshRegion meshRegionBandsFromDistField(const CMesh& mesh, int seed_vert, const std::vector<double>& dist_field, const std::vector<double>& thresholds);
BandedMeshRegion meshRegionBandsFromDistField(const CMesh& mesh, int seed_vert, const std::vector<double>& dist_field, double max_threshold, int band_num);

} // end of namespace

#endif