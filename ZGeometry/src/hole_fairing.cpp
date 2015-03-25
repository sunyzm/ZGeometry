#include "hole_fairing.h"
#include <ZGeom/util.h>
#include <ZGeom/SparseRepresentation.h>
#include <map>
#include "global.h"
#include "MeshLaplacian.h"
#include "GeometryApproximation.h"

using namespace std;
using ZGeom::SparseMatrixd;
using ZGeom::DenseMatrixd;
using ZGeom::VecNd;


MeshCoordinates least_square_inpainting(CMesh& mesh, const std::vector<int>& anchor_verts, double anchor_weight)
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

    MeshCoordinates result;
    result.fromDenseMatrix(matX);
    return result;
}

MeshCoordinates least_square_hole_inpainting(CMesh& mesh, const ZGeom::MeshRegion& hole_region, int anchor_ring, double anchor_weight)
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
    
    MeshCoordinates faired_sub_coord = least_square_inpainting(submesh, control_verts, anchor_weight);
    MeshCoordinates result(mesh.getVertCoordinates());
    for (int sub_vIdx = 0; sub_vIdx < inside_vert_count; ++sub_vIdx) {
        result.setVertCoordinate(newVert2oldVert[sub_vIdx], faired_sub_coord[sub_vIdx]);
    }

    return result;
}

ZGeom::DenseMatrixd matlab_inpaintL1LS(const DenseMatrixd& matCoord, const DenseMatrixd& matDict, const std::vector<int>& vMissingIdx, double lambda, double tol)
{
    g_engineWrapper.addDenseMat(matCoord, "coord");
    g_engineWrapper.addDenseMat(matDict, "dict");

    ZGeom::VecNd vecMissing(vMissingIdx.size());
    for (int i = 0; i < vecMissing.size(); ++i) vecMissing[i] = double(vMissingIdx[i] + 1);
    g_engineWrapper.addColVec(vecMissing, "missing_idx");
    g_engineWrapper.addDoubleScalar(lambda, "lambda");
    g_engineWrapper.addDoubleScalar(tol, "tol");

    g_engineWrapper.eval("[coord_est,~,err] =  zmesh_inpaint_l1ls(dict, coord, missing_idx, lambda, tol);");
    ZGeom::DenseMatrixd result = g_engineWrapper.getDenseMat("coord_est");

    return result;
}

MeshCoordinates l1_ls_inpainting(CMesh& mesh, const std::vector<int>& anchors, ParaL1LsInpainting& para)
{
    const int totalVertCount = mesh.vertCount();
    MeshLaplacian graphLaplacian;
    graphLaplacian.constructUmbrella(&mesh);
    int eigenCount = para.eigen_count;    // -1 means full decomposition
    ZGeom::EigenSystem es;
    graphLaplacian.meshEigenDecompose(eigenCount, &g_engineWrapper, es);
    ZGeom::Dictionary dictMHB;
    computeDictionary(DT_Fourier, es, dictMHB);

    ZGeom::DenseMatrixd matCoordOld = mesh.getVertCoordinates().toDenseMatrix();
    ZGeom::DenseMatrixd matDict = dictMHB.toDenseMatrix();
    const vector<int>& selected_verts = anchors;

    ZGeom::DenseMatrixd matCoordInpainted = matlab_inpaintL1LS(matCoordOld, matDict, selected_verts, para.lambda, para.tol);
    
    MeshCoordinates coordInpainted;
    coordInpainted.fromDenseMatrix(matCoordInpainted);
    return coordInpainted;
}

MeshCoordinates l1_ls_hole_inpainting(CMesh& mesh, const ZGeom::MeshRegion& hole_region, ParaL1LsInpainting& para)
{
    const vector<int> anchor_verts = ZGeom::meshRegionSurroundingVerts(mesh, hole_region, para.fitting_ring);
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

    MeshCoordinates faired_sub_coord = l1_ls_inpainting(submesh, control_verts, para);

    MeshCoordinates result(mesh.getVertCoordinates());
    for (int sub_vIdx = 0; sub_vIdx < inside_vert_count; ++sub_vIdx) {
        result.setVertCoordinate(newVert2oldVert[sub_vIdx], faired_sub_coord[sub_vIdx]);
    }
    return result;
}
