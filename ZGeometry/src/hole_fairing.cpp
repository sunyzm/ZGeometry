#include "hole_fairing.h"
#include <ZGeom/util.h>
#include <map>
#include "MeshLaplacian.h"

using namespace std;
using ZGeom::SparseMatrixd;
using ZGeom::DenseMatrixd;
using ZGeom::VecNd;


MeshCoordinates least_square_fairing_whole_mesh(CMesh& mesh, const std::vector<int>& anchor_verts, double anchor_weight)
{
    const int totalVertCount = mesh.vertCount();
    const MeshCoordinates coordOld = mesh.getVertCoordinates();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();
    MeshLaplacian matLaplacian;
    matLaplacian.constructTutte(&mesh);

    const vector<int> &vControlVerts = anchor_verts;
    set<int> setControlVerts{ anchor_verts.begin(), anchor_verts.end() };
    int controlVertCount = (int)anchor_verts.size();
    double controlWeight = anchor_weight;
    vector<int> rowIdxA, colIdxA;
    vector<double> valsA;
    matLaplacian.getSparseMatrix().convertToCOO(rowIdxA, colIdxA, valsA, ZGeom::MAT_FULL);
    for (int i = 0; i < controlVertCount; ++i) {
        rowIdxA.push_back(totalVertCount + i + 1);
        colIdxA.push_back(vControlVerts[i] + 1);
        valsA.push_back(controlWeight);
    }
    SparseMatrixd matA(totalVertCount + controlVertCount, totalVertCount);
    matA.convertFromCOO(totalVertCount + controlVertCount, totalVertCount, rowIdxA, colIdxA, valsA);

    DenseMatrixd matB(totalVertCount + controlVertCount, 3);
    for (int i = 0; i < controlVertCount; ++i) {
        for (int j = 0; j < 3; ++j)
            matB(totalVertCount + i, j) = controlWeight * vOriginalCoords[j][vControlVerts[i]];
    }

    g_engineWrapper.addSparseMat(matA, "matA");
    g_engineWrapper.addDenseMat(matB, "matB");
    g_engineWrapper.eval("matX=matA\\matB;");
    DenseMatrixd matX = g_engineWrapper.getDenseMat("matX");

    vector<VecNd> vCoordRecovered = vOriginalCoords;
    for (int i = 0; i < totalVertCount; ++i) {
        if (setControlVerts.find(i) == setControlVerts.end()) {
            // replace vertices inside holes with recovered coordinate
            for (int j = 0; j < 3; ++j)  vCoordRecovered[j][i] = matX(i, j);
        }
    }

    return MeshCoordinates(totalVertCount, vCoordRecovered);
}

MeshCoordinates least_square_fairing(CMesh& mesh, const ZGeom::MeshRegion& hole_region, int anchor_ring, double anchor_weight)
{
    ZGeom::logic_assert(anchor_ring >= 1 && anchor_weight >= 0, "Illegal parameter!");
    
    const vector<int> anchor_verts = ZGeom::meshRegionSurroundingVerts(mesh, hole_region, anchor_ring);
    vector<int> submesh_verts;
    int inside_vert_count = (int)hole_region.vert_inside.size();
    int anchor_vert_count = (int)anchor_verts.size();
    vector<int> newVert2oldVert(inside_vert_count);
    int newVertIdx(0);
    for (int vi : hole_region.vert_inside) {
        submesh_verts.push_back(vi);
        newVert2oldVert[newVertIdx++] = vi;
    }
    for (int vi : anchor_verts) submesh_verts.push_back(vi);
    
    CMesh submesh;
    mesh.getSubMesh(submesh_verts, "hole_mesh", submesh);
    vector<int> control_verts;
    for (int i = inside_vert_count; i < (int)submesh_verts.size(); ++i)
        control_verts.push_back(i);
    
    MeshCoordinates faired_sub_coord = least_square_fairing_whole_mesh(submesh, control_verts, anchor_weight);
    MeshCoordinates faired_whole_coord(mesh.getVertCoordinates());
    for (int sub_vIdx = 0; sub_vIdx < inside_vert_count; ++sub_vIdx) {
        faired_whole_coord.setVertCoordinate(newVert2oldVert[sub_vIdx], faired_sub_coord[sub_vIdx]);
    }

    return faired_whole_coord;
}
