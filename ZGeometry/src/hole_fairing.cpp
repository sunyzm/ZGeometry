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

MeshCoordinates thin_plate_energy_fairing(CMesh& mesh, int free_vert_count)
{
    const int totalVertCount = mesh.vertCount();

    const MeshCoordinates coordOld = mesh.getVertCoordinates();
    vector<VecNd> vOriginalCoords = coordOld.to3Vec();
    MeshLaplacian mesh_laplacian;
    mesh_laplacian.constructCotFormula(&mesh);
    SparseMatrixd matL = mesh_laplacian.getSparseMatrix();
    SparseMatrixd matL2;
    mulMatMat(matL, matL, matL2);

    vector<int> rowIdxL2, colIdxL2;
    vector<double> valsL2;
    matL2.convertToCOO(rowIdxL2, colIdxL2, valsL2, ZGeom::MAT_FULL);

    vector<int> rowIdxA, colIdxA; vector<double> valsA;
    for (int i = 0; i < (int)rowIdxL2.size(); ++i) {
        if (rowIdxL2[i] < free_vert_count && colIdxL2[i] < free_vert_count) {
            rowIdxA.push_back(rowIdxL2[i]);
            colIdxA.push_back(colIdxL2[i]);
            valsA.push_back(valsL2[i]);
        }
    }
    SparseMatrixd matA;
    matA.convertFromCOO(free_vert_count, free_vert_count, rowIdxA, colIdxA, valsA);

    DenseMatrixd matB(free_vert_count, 3);
    for (int m = 0; m < 3; ++m) {
        for (int i = 0; i < free_vert_count; ++i) {
            for (int j = free_vert_count; j < totalVertCount; ++j)
                matB(i, m) -= matL2(i, j) * vOriginalCoords[m][j];
        }
    }

    g_engineWrapper.addSparseMat(matA, "matA");
    g_engineWrapper.addDenseMat(matB, "matB");
    g_engineWrapper.eval("matX=matA\\matB;");
    DenseMatrixd matX = g_engineWrapper.getDenseMat("matX");

    vector<VecNd> vNewCoord = vOriginalCoords;
    for (int m = 0; m < 3; ++m)
        for (int i = 0; i < free_vert_count; ++i)
            vNewCoord[m][i] = matX(i, m);
    
    MeshCoordinates result(totalVertCount, vNewCoord);
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

MeshCoordinates thin_plate_energy_hole_inpainting(CMesh& mesh, const ZGeom::MeshRegion& hole_region, int nIter, double eps)
{
    const vector<int> anchor_verts = ZGeom::meshRegionSurroundingVerts(mesh, hole_region, 2);
    vector<int> submesh_verts;
    int inside_vert_count = (int)hole_region.vert_inside.size();
    int anchor_vert_count = (int)anchor_verts.size();

    vector<int> newVert2oldVert(inside_vert_count);
    for (int i = 0; i < inside_vert_count; ++i) {
        int vi = hole_region.vert_inside[i];
        submesh_verts.push_back(vi);
        newVert2oldVert[i] = vi;
    }
    for (int vi : anchor_verts) submesh_verts.push_back(vi);

    CMesh submesh;
    mesh.getSubMesh(submesh_verts, "hole_mesh", submesh);

    MeshCoordinates faired_sub_coord = thin_plate_energy_fairing(submesh, inside_vert_count);

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

MeshCoordinates l1_ls_inpainting(CMesh& mesh, const std::vector<int>& missing_idx, ParaL1LsInpainting& para)
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

    ZGeom::DenseMatrixd matCoordInpainted = matlab_inpaintL1LS(matCoordOld, matDict, missing_idx, para.lambda, para.tol);
    
    MeshCoordinates coordInpainted;
    coordInpainted.fromDenseMatrix(matCoordInpainted);
    return coordInpainted;
}

MeshCoordinates l1_ls_hole_inpainting(CMesh& mesh, const ZGeom::MeshRegion& hole_region, ParaL1LsInpainting& para)
{
    if (para.fitting_ring <= 0) {
        MeshCoordinates faired_coord = l1_ls_inpainting(mesh, hole_region.vert_inside, para);
        MeshCoordinates result(mesh.getVertCoordinates());
        for (int vi : hole_region.vert_inside)
            result.setVertCoordinate(vi, faired_coord[vi]);
        return result;
    }

    const vector<int> anchor_verts = ZGeom::meshRegionSurroundingVerts(mesh, hole_region, para.fitting_ring);
    vector<int> submesh_verts;
    int inside_vert_count = (int)hole_region.vert_inside.size();
    int anchor_vert_count = (int)anchor_verts.size();
    vector<int> newVert2oldVert(inside_vert_count);
    vector<int> missing_verts;
    int newVertIdx(0);
    for (int vi : hole_region.vert_inside) {
        submesh_verts.push_back(vi);
        missing_verts.push_back(newVertIdx);
        newVert2oldVert[newVertIdx++] = vi;
    }
    for (int vi : anchor_verts) submesh_verts.push_back(vi);

    CMesh submesh;
    mesh.getSubMesh(submesh_verts, "hole_mesh", submesh);

    MeshCoordinates faired_sub_coord = l1_ls_inpainting(submesh, missing_verts, para);

    MeshCoordinates result(mesh.getVertCoordinates());
    for (int sub_vIdx = 0; sub_vIdx < inside_vert_count; ++sub_vIdx) {
        result.setVertCoordinate(newVert2oldVert[sub_vIdx], faired_sub_coord[sub_vIdx]);
    }
    return result;
}
