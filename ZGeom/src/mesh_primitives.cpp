#include "mesh_primitives.h"

using namespace std;

namespace ZGeom {
    
void MeshRegion::determineBoundaryHalfEdges(const CMesh& mesh)
{
    std::set<int> faces_in_hole{ face_inside.begin(), face_inside.end() };

    std::set<int> boundary_he;
    for (int fi : faces_in_hole) {
        const CFace* f = mesh.getFace(fi);
        for (CHalfEdge* he : f->getAllHalfEdges()) {
            if (he->twinHalfEdge() == nullptr || !setHas(faces_in_hole, he->twinHalfEdge()->getAttachedFace()->getFaceIndex()))
                boundary_he.insert(he->getIndex());
        }
    }

    he_on_boundary = std::vector<int>(boundary_he.begin(), boundary_he.end());
}

MeshRegion& MeshRegion::operator=(MeshRegion&& hb)
{
    vert_on_boundary = std::move(hb.vert_on_boundary);
    he_on_boundary = std::move(hb.he_on_boundary);
    vert_inside = std::move(hb.vert_inside);
    face_inside = std::move(hb.face_inside);
    is_outer_boundary = hb.is_outer_boundary;
    adjacent_edge_length = hb.is_outer_boundary;
    return *this;
}


std::vector<int> vertSurroundingVerts(const CMesh& mesh, const vector<int>& vert_inside, int ring)
{
    if (ring <= 0) {    // return all remaining vertices 
        set<int> inside_verts(vert_inside.begin(), vert_inside.end());
        vector<int> result;
        for (int i = 0; i < mesh.vertCount(); ++i) {
            if (!setHas(inside_verts, i)) result.push_back(i);
        }
        return result;
    }

    set<int> result;
    set<int> considered_vert(vert_inside.begin(), vert_inside.end());
    set<int> cur_ring = considered_vert;

    for (int level = 1; level <= ring; ++level)
    {
        set<int> new_ring;
        for (int cur_vi : cur_ring) {
            vector<int> cur_neighbor_vert = mesh.getVertNeighborVerts(cur_vi, 1, false);
            for (int vj : cur_neighbor_vert) {
                if (!setHas(considered_vert, vj)) {
                    new_ring.insert(vj);
                    considered_vert.insert(vj);
                }
            }
        }
        for (int new_vi : new_ring) {
            result.insert(new_vi);
        }
        cur_ring = new_ring;
    }

    return vector<int>(result.begin(), result.end());
}


std::set<int> meshMultiVertsAdjacentVerts(const CMesh& mesh, const std::vector<int>& vert, int ring, bool inclusive /*= true*/)
{
    set<int> nbr(vert.begin(), vert.end());

    std::set<int> nbp = nbr;
    for (int r = 1; r <= ring; ++r) {
        std::set<int> nbn;
        for (int vn : nbp) {
            const CVertex* vStart = mesh.vert(vn);
            for (auto he : vStart->m_HalfEdges) {
                int endv = he->getVertIndex(1);
                if (nbr.find(endv) == nbr.end()) {
                    nbr.insert(endv);
                    nbn.insert(endv);
                }
                // to avoid boundary vertex being ignored
                if (he->m_eNext) {
                    int endv2 = he->m_eNext->vert(1)->m_vIndex;
                    if (nbr.find(endv2) == nbr.end()) {
                        nbr.insert(endv2);
                        nbn.insert(endv2);
                    }
                }
            }
        }
        nbp = nbn;
    }

    if (!inclusive) {
        for (int vi : vert) nbr.erase(vi);
    }

    return nbr;
}

std::vector<int> getFaceEncompassedByVerts(const CMesh& mesh, const std::vector<int>& verts)
{
    set<int> face_idx;
    set<int> all_verts(verts.begin(), verts.end());
    for (int vi : verts) {
        for (CFace* fi : mesh.vert(vi)->getAdjacentFaces()) {
            if (setHas(face_idx, fi->getFaceIndex())) continue;
            if (setHasAll(all_verts, fi->getAllVertIdx()))
                face_idx.insert(fi->getFaceIndex());
        }
    }
    return vector<int>(face_idx.begin(), face_idx.end());
}

MeshRegion meshRegionFromVerts(CMesh& mesh, const std::vector<int>& inside_verts)
{
    MeshRegion result;
    result.vert_inside = inside_verts;
    for (int vi : inside_verts) {
        for (CFace* f : mesh.vert(vi)->getAdjacentFaces())
            result.face_inside.push_back(f->getFaceIndex());
    }

    result.vert_on_boundary = vertSurroundingVerts(mesh, result.vert_inside, 1);
    result.determineBoundaryHalfEdges(mesh);

    return result;
}


ZGeom::MeshRegion meshRegionFromDistField(CMesh& mesh, const std::vector<double>& dist_field, int seed, std::function<bool(double)> func_in_region_test)
{
    assert(mesh.vertCount() == dist_field.size());

    set<int> candidate_vert{ seed };
    set<int> vert_visited;
    set<int> vert_inside;
    while (!candidate_vert.empty()) {
        int vi = *candidate_vert.begin();
        candidate_vert.erase(vi);
        vert_visited.insert(vi);
        if (func_in_region_test(dist_field[vi])) {
            vert_inside.insert(vi);
            vector<int> neighbors = mesh.getVertNeighborVerts(vi, 1);
            for (int vj : neighbors) {
                if (vert_visited.find(vj) == vert_visited.end()) {
                    candidate_vert.insert(vj);
                }
            }
        }
    }

    MeshRegion result = meshRegionFromVerts(mesh, vector<int>{vert_inside.begin(), vert_inside.end()});
    return result;
}

ZGeom::BandedMeshRegions meshRegionBandsFromDistField(CMesh& mesh, int seed_vert,
        const std::vector<double>& dist_field, const std::vector<double>& thresholds)
{
    int ring_num = (int)thresholds.size();
    
    // check validity of thresholds
    runtime_assert(thresholds[0] > 0);
    for (int k = 1; k < ring_num; ++k) {
        runtime_assert(thresholds[k] > thresholds[k - 1]);
    }

    set<int> candidate_vert{ seed_vert };
    set<int> vert_visited;
    vector<set<int>> vert_in_ring(ring_num);
    while (!candidate_vert.empty()) {
        int vi = *candidate_vert.begin();
        candidate_vert.erase(vi);
        vert_visited.insert(vi);

        for (int r = 0; r < ring_num; ++r) {
            if (dist_field[vi] <= thresholds[r]) {
                vert_in_ring[r].insert(vi);
                for (int vj : mesh.getVertNeighborVerts(vi, 1)) {
                    if (vert_visited.find(vj) == vert_visited.end()) {
                        candidate_vert.insert(vj);
                    }
                }
                break;
            }
        }
    }
    
    BandedMeshRegions result;
    result.resize(ring_num);
    for (int r = 0; r < ring_num; ++r) {
        result.band_thresholds = thresholds;
        result.band_verts[r] = vector<int>(vert_in_ring[r].begin(), vert_in_ring[r].end());
    }
    return result;
}

ZGeom::BandedMeshRegions meshRegionBandsFromDistField(CMesh& mesh, int seed_vert,
        const std::vector<double>& dist_field, double max_threshold, int band_num)
{
    assert(band_num >= 1 && max_threshold > 0);
    vector<double> band_thresholds(band_num);
    for (int k = 1; k <= band_num; ++k) {
        band_thresholds[k - 1] = max_threshold * double(k) / double(band_num);
    }
    return meshRegionBandsFromDistField(mesh, seed_vert, dist_field, band_thresholds);
}

}