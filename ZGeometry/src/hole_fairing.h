#pragma once
#include <ZGeom/ZGeom.h>
#include <vector>

MeshCoordinates least_square_inpainting(CMesh& mesh, const std::vector<int>& anchor_verts, double anchor_weight);
MeshCoordinates least_square_hole_inpainting(CMesh& mesh, const std::vector<ZGeom::MeshRegion>& hole_region, int anchor_ring, double anchor_weight);


struct ParaL1LsInpainting
{
    int fitting_ring;   // -1 means use whole mesh
    double lambda;      // relative weight of for LS term
    double tol;
    int eigen_count;

    ParaL1LsInpainting() : 
        fitting_ring(0), lambda(1e-3), tol(1e-3), eigen_count(-1) {}
};


MeshCoordinates l1_ls_inpainting(CMesh& mesh, const std::vector<int>& anchors, ParaL1LsInpainting& para);
MeshCoordinates l1_ls_hole_inpainting(CMesh& mesh, const std::vector<ZGeom::MeshRegion>& hole_region, 
        ParaL1LsInpainting& para);
ZGeom::DenseMatrixd matlab_inpaintL1LS(const ZGeom::DenseMatrixd& matCoord, 
        const ZGeom::DenseMatrixd& matDict, const std::vector<int>& vMissingIdx, 
        double lambda, double tol);
MeshCoordinates thin_plate_energy_fairing(CMesh& mesh, int free_vert_count);
MeshCoordinates thin_plate_energy_hole_inpainting(CMesh& mesh, 
        const ZGeom::MeshRegion& hole_region, int nIter, double eps);
MeshCoordinates meshDLRS(CMesh& mesh, double lambda, const std::vector<int>& selected_verts, int ring);
MeshCoordinates meshHoleDLRS(CMesh& mesh, const std::vector<ZGeom::MeshRegion>& hole_regions, 
        double lambda, int ring);