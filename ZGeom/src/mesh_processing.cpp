#include "mesh_processing.h"
#include <cstdlib>
#include <random>
#include <ppl.h>
#include "MatVecArithmetic.h"
#include "arithmetic.h"

using namespace std;

namespace ZGeom {

std::function<double(double, double)> heatGen = [](double lambda, double tau){
    return std::exp(-lambda*tau);
};

std::function<double(double, double)> mhwGen = [](double lambda, double tau){
    return lambda * std::exp(-lambda*tau);
};

SparseMatrix<double> constructMeshMatrix(const CMesh& mesh, MeshMatrixType mmt)
{
    SparseMatrix<double> meshMat;
    const int vertCount = mesh.vertCount();

    if (mmt == MM_DEGREE) {
        vector<int> vDiagDegree(vertCount);
        for (int i = 0; i < vertCount; ++i) {
            vDiagDegree[i] = (int)mesh.getVertNeighborVerts(i, 1, false).size();
        }
        meshMat.convertFromDiagonal(vDiagDegree);
    }
    else if (mmt == MM_INV_DEGREE) {
        vector<double> vDiagDegree(vertCount);
        for (int i = 0; i < vertCount; ++i) {
            vDiagDegree[i] = 1.0 / (int)mesh.getVertNeighborVerts(i, 1, false).size();
        }
        meshMat.convertFromDiagonal(vDiagDegree);
    }
    else if (mmt == MM_ADJACENCY) {
        vector<tuple<int, int, double> > vElem;
        for (int i = 0; i < vertCount; ++i) {
            vector<int> vNeighbors = mesh.getVertNeighborVerts(i, 1, false);
            int valence = (int)vNeighbors.size();
            for (int j = 0; j < valence; ++j)
                vElem.push_back(std::make_tuple(i + 1, vNeighbors[j] + 1, 1.));
        }
        meshMat.convertFromCOO(vertCount, vertCount, vElem);
    }
    else if (mmt == MM_INV_LENGTH) {

    }
    else if (mmt == MM_WALK) {
        // W = D^(-1)*A
        SparseMatrixd matInvD = constructMeshMatrix(mesh, MM_INV_DEGREE);
        SparseMatrixd matA = constructMeshMatrix(mesh, MM_ADJACENCY);
        mulMatMat(matInvD, matA, meshMat);
    }
    else if (mmt == MM_GRAPH_LAPLACE) {
        // L = D - A
        SparseMatrixd matD = constructMeshMatrix(mesh, MM_DEGREE);
        SparseMatrixd matA = constructMeshMatrix(mesh, MM_ADJACENCY);
        addMatMat(matD, matA, -1.0, meshMat);
    }
    else if (mmt == MM_NORMALIZED_GRAPH_LAPLACE) {
        // Ln = D^(-0.5)*L*D^(-0.5) = I - D^(-0.5)*A*D^(-0.5)
        SparseMatrixd matD = constructMeshMatrix(mesh, MM_DEGREE);
        vector<double> vInvSqD = matD.getDiagonal();
        for (double& v : vInvSqD) v = std::pow(v, -0.5); // D^(-0.5)
        meshMat = constructMeshMatrix(mesh, MM_GRAPH_LAPLACE);
        for (auto& elem : meshMat.allElements()) {
            elem.val() *= vInvSqD[elem.row() - 1] * vInvSqD[elem.col() - 1];
        }
    }

    return meshMat;
}

DenseMatrixd calSpectralKernelMatrix(const EigenSystem& es, double t, std::function<double(double, double)> genfunc)
{    
    const int vertCount = es.eigVecSize();
    const int eigCount = es.eigVecCount();
    std::vector<double> vDiag(eigCount);
    for (int k = 0; k < eigCount; ++k) vDiag[k] = genfunc(es.getEigVal(k), t);
    DenseMatrixd matEigVecs = es.toDenseMatrix();

    DenseMatrixd hk(vertCount, vertCount);
    quadricFormAMP(vertCount, eigCount, matEigVecs.raw_ptr(), vDiag.data(), hk.raw_ptr());

    return hk;
}

std::vector<double> calSpectralKernelSignature(const EigenSystem& es, double t, std::function<double(double, double)> genfunc)
{
    const int vertCount = es.eigVecSize();
    const int eigCount = es.eigVecCount();
    std::vector<double> vDiag(eigCount);
    for (int k = 0; k < eigCount; ++k) vDiag[k] = genfunc(es.getEigVal(k), t);

    vector<double> result(vertCount, 0);
    for (int i = 0; i < vertCount; ++i) {
        for (int k = 0; k < eigCount; ++k) {
            const VecNd& phi = es.getEigVec(k);
            result[i] += vDiag[k] * phi[i] * phi[i];
        }
    }

    return result;
}

DenseMatrixd calHeatKernelMatrix(const EigenSystem& hb, double t)
{
    return calSpectralKernelMatrix(hb, t, heatGen);
}

std::vector<double> calHeatKernelSignature(const EigenSystem& hb, double t)
{
    return calSpectralKernelSignature(hb, t, heatGen);
}

double calHK(const EigenSystem& es, int i, int j, double t)
{
    double sum = 0;
    for (int k = 0; k < es.eigVecCount(); ++k) {
        auto phi = es.getEigVec(k);
        sum += std::exp(-es.getEigVal(k) * t) * phi[i] * phi[j];
    }
    return sum;
}

double calHeatTrace(const EigenSystem& es, double t)
{
    double sum(0);
    for (int k = 0; k < es.eigVecCount(); ++k)
        sum += std::exp(-es.getEigVal(k) * t);
    return sum;
}

std::vector<Vec3d> computeMeshVertNormals(const CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    int vertCount = mesh.vertCount(), faceCount = mesh.faceCount();
    const vector<Vec3d> faceNormals = computeMeshFaceNormals(mesh);
    vector<double> faceAreas(faceCount);
    for (int i = 0; i < faceCount; ++i) faceAreas[i] = mesh.getFace(i)->calArea();

    vector<Vec3d> vertNormals(vertCount);
    for (int vi = 0; vi < vertCount; ++vi) {
        vector<int> neighborFaceIdx = mesh.getVertexAdjacentFaceIdx(vi, 1);
        Vec3d normalSum(0, 0, 0);
        for (int fj : neighborFaceIdx) {
            double weight(1.0);
            if (vnc == VN_AREA_WEIGHT)
                weight = faceAreas[fj];
            else if (vnc == VN_CENTER_DIST_WEIGHT)
                weight = 1.0 / (mesh.getFace(fj)->calBarycenter() - mesh.vertPos(vi)).length();
            normalSum += weight * faceNormals[fj];
        }
        vertNormals[vi] = normalSum.normalize();
    }

    return vertNormals;
}



std::vector<int> randomHoleVertex(const CMesh& mesh, int hole_size, int seed /*= -1*/)
{
    using std::vector;
    assert(hole_size >= 1);
    int N = mesh.vertCount();
    std::default_random_engine generator;
    if (seed == -1) {
        std::uniform_int_distribution<int> distribution(0, N - 1);
        seed = distribution(generator);
    }

    std::set<int> vertInHole;
    std::vector<int> vertFeasible;

    vertInHole.insert(seed);
    vertFeasible.push_back(seed);


    while (vertInHole.size() < hole_size && !vertFeasible.empty()) {
        std::uniform_int_distribution<int> distr1(0, (int)vertFeasible.size() - 1);
        int selPos = distr1(generator);
        int selIdx = vertFeasible[selPos];
        vector<int> neighbors = mesh.getVertNeighborVerts(selIdx, 1, false);
        vector<int> feasibleNeighbors;
        for (int vIdx : neighbors) {
            if (vertInHole.find(vIdx) == vertInHole.end())
                feasibleNeighbors.push_back(vIdx);
        }
        if (feasibleNeighbors.empty()) {
            vertFeasible.erase(vertFeasible.begin() + selPos);
        }
        else {
            std::uniform_int_distribution<int> distr2(0, (int)feasibleNeighbors.size() - 1);
            int selNeighbor = feasibleNeighbors[distr2(generator)];
            vertInHole.insert(selNeighbor);
            if (feasibleNeighbors.size() == 1) vertFeasible.erase(vertFeasible.begin() + selPos);
            vertFeasible.push_back(selNeighbor);
        }
    }

    return vector<int>(vertInHole.begin(), vertInHole.end());
}

std::vector<int> randomHoleVertex(const CMesh& mesh, int total_size, const std::vector<int>& seeds)
{
    using std::vector;
    assert(total_size >= seeds.size());
    int N = mesh.vertCount();
    std::default_random_engine generator;

    std::set<int> vertInHole{ seeds.begin(), seeds.end() };
    std::vector<int> vertFeasible{ seeds.begin(), seeds.end() };

    while (vertInHole.size() < total_size && !vertFeasible.empty()) {
        std::uniform_int_distribution<int> distr1(0, (int)vertFeasible.size() - 1);
        int selPos = distr1(generator);
        int selIdx = vertFeasible[selPos];
        vector<int> neighbors = mesh.getVertNeighborVerts(selIdx, 1, false);
        vector<int> feasibleNeighbors;
        for (int vIdx : neighbors) {
            if (vertInHole.find(vIdx) == vertInHole.end())
                feasibleNeighbors.push_back(vIdx);
        }
        if (feasibleNeighbors.empty()) {
            vertFeasible.erase(vertFeasible.begin() + selPos);
        }
        else {
            std::uniform_int_distribution<int> distr2(0, (int)feasibleNeighbors.size() - 1);
            int selNeighbor = feasibleNeighbors[distr2(generator)];
            vertInHole.insert(selNeighbor);
            if (feasibleNeighbors.size() == 1) vertFeasible.erase(vertFeasible.begin() + selPos);
            vertFeasible.push_back(selNeighbor);
        }
    }

    return vector<int>(vertInHole.begin(), vertInHole.end());

}

std::vector<Vec3d> computeMeshFaceNormals(const CMesh& mesh)
{
    int faceCount = mesh.faceCount();
    vector<Vec3d> vecNormals(faceCount);

    concurrency::parallel_for(0, faceCount, [&](int i) {
        vecNormals[i] = mesh.getFace(i)->calcNormal();
    });

    return vecNormals;
}

double calBiharmonicDist(const EigenSystem& es, int v1, int v2)
{
    // sqrt(sum_{k=1}^{K}(phi_k(i)-phi_k(j))/lambda_k)^2 )
    double sum = 0;
    for (int k = 1; k < es.eigVecCount(); ++k) {
        sum += pow((es.getEigVec(k)[v1] - es.getEigVec(k)[v2]) / es.getEigVal(k), 2);
    }
    return std::sqrt(sum);
}

void getMeshGraphCSR(const CMesh& mesh, std::vector<int>& xadj, std::vector<int>& adjncy)
{
    int vertNum = mesh.vertCount();
    xadj.resize(vertNum + 1);
    adjncy.clear();
    xadj[0] = 0;
    for (int i = 0; i < vertNum; ++i) {
        std::vector<int> vNbr = mesh.getVertNeighborVerts(i, 1, false);
        for (int j : vNbr) adjncy.push_back(j);
        xadj[i + 1] = (int)adjncy.size();
    }
}

/* conf. "MESH - MEASURING ERRORS BETWEEN SURFACES USING THE HAUSDORFF (2002)" */
double computeMeanHausdorffDistance(const CMesh& mesh1, const CMesh& mesh2)
{
    double result{ 0 };


    return result;
}

/* conf. "Distance between point and triangle in 3D 
/* adapted from GTEngine\Include\GteDistPointTriangleExact */
ResultDistPointTriangle distPointTriangle(Vec3d point, const std::vector<Vec3d>& triangle)
{
    Vec3d diff = point - triangle[0];
    Vec3d edge0 = triangle[1] - triangle[0];
    Vec3d edge1 = triangle[2] - triangle[0];
    double a00 = dot(edge0, edge0);
    double a01 = dot(edge0, edge1);
    double a11 = dot(edge1, edge1);
    double b0 = -dot(diff, edge0);
    double b1 = -dot(diff, edge1);
    const double zero = 0, one = (double)1;
    double det = a00 * a11 - a01 * a01;
    double t0 = a01 * b1 - a11 * b0;
    double t1 = a01 * b0 - a00 * b1;

    if (t0 + t1 <= det)
    {
        if (t0 < zero)
        {
            if (t1 < zero)  // region 4
            {
                if (b0 < zero)
                {
                    t1 = zero;
                    if (-b0 >= a00)  // V0
                    {
                        t0 = one;
                    }
                    else  // E01
                    {
                        t0 = -b0 / a00;
                    }
                }
                else
                {
                    t0 = zero;
                    if (b1 >= zero)  // V0
                    {
                        t1 = zero;
                    }
                    else if (-b1 >= a11)  // V2
                    {
                        t1 = one;
                    }
                    else  // E20
                    {
                        t1 = -b1 / a11;
                    }
                }
            }
            else  // region 3
            {
                t0 = zero;
                if (b1 >= zero)  // V0
                {
                    t1 = zero;
                }
                else if (-b1 >= a11)  // V2
                {
                    t1 = one;
                }
                else  // E20
                {
                    t1 = -b1 / a11;
                }
            }
        }
        else if (t1 < zero)  // region 5
        {
            t1 = zero;
            if (b0 >= zero)  // V0
            {
                t0 = zero;
            }
            else if (-b0 >= a00)  // V1
            {
                t0 = one;
            }
            else  // E01
            {
                t0 = -b0 / a00;
            }
        }
        else  // region 0, interior
        {
            double invDet = one / det;
            t0 *= invDet;
            t1 *= invDet;
        }
    }
    else
    {
        double tmp0, tmp1, numer, denom;

        if (t0 < zero)  // region 2
        {
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)  // V1
                {
                    t0 = one;
                    t1 = zero;
                }
                else  // E12
                {
                    t0 = numer / denom;
                    t1 = one - t0;
                }
            }
            else
            {
                t0 = zero;
                if (tmp1 <= zero)  // V2
                {
                    t1 = one;
                }
                else if (b1 >= zero)  // V0
                {
                    t1 = zero;
                }
                else  // E20
                {
                    t1 = -b1 / a11;
                }
            }
        }
        else if (t1 < zero)  // region 6
        {
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0)
            {
                numer = tmp1 - tmp0;
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)  // V2
                {
                    t1 = one;
                    t0 = zero;
                }
                else  // E12
                {
                    t1 = numer / denom;
                    t0 = one - t1;
                }
            }
            else
            {
                t1 = zero;
                if (tmp1 <= zero)  // V1
                {
                    t0 = one;
                }
                else if (b0 >= zero)  // V0
                {
                    t0 = zero;
                }
                else  // E01
                {
                    t0 = -b0 / a00;
                }
            }
        }
        else  // region 1
        {
            numer = a11 + b1 - a01 - b0;
            if (numer <= zero)  // V2
            {
                t0 = zero;
                t1 = one;
            }
            else
            {
                denom = a00 - ((double)2)*a01 + a11;
                if (numer >= denom)  // V1
                {
                    t0 = one;
                    t1 = zero;
                }
                else  // 12
                {
                    t0 = numer / denom;
                    t1 = one - t0;
                }
            }
        }
    }

    ResultDistPointTriangle result;

    result.closestPoint = triangle[0] + t0 * edge0 + t1 * edge1;
    diff = point - result.closestPoint;
    result.distance = diff.length();
    return result;
}

ZGeom::ResultDistPointPlane distPointPlane(Vec3d point, const Plane3& plane)
{
    ResultDistPointPlane result;
    result.signedDistance = dot(plane.normal, point) - plane.constant;
    result.distance = std::fabs(result.signedDistance);
    result.closestPoint = point - result.signedDistance * plane.normal;
    return result;
}

/* Mixed-Voronoi weighting: "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds (2002)" */
std::vector<double> computeMeshVertArea(const CMesh& mesh, MeshVertAreaScheme scheme /*= VA_VORONOI*/)
{
    int nVerts = mesh.vertCount(), nFaces = mesh.faceCount();
    vector<double> result(nVerts, 0);

    if (scheme == VA_UNIFORM) 
    {
        for (int fi = 0; fi < nFaces; ++fi) {
            double area_fi = mesh.getFace(fi)->calArea();
            vector<int> faceVertIdx = mesh.getFace(fi)->getAllVertIdx();
            for (int i = 0; i < 3; ++i) {
                result[faceVertIdx[i]] += area_fi / 3.;
            }
        }
    }
    else if (scheme == VA_MIXED_VORONOI)
    {
        for (int vIndex = 0; vIndex < nVerts; ++vIndex)
        {
            const CVertex* vi = mesh.vert(vIndex);
            double amix = 0;
            for (CHalfEdge* e0 : vi->m_HalfEdges) {
                CHalfEdge* e1 = e0->m_eNext;
                CHalfEdge* e2 = e1->m_eNext;
                double len0 = e0->length();
                double len1 = e1->length();
                double len2 = e2->length();
                amix += ZGeom::calMixedTriArea(len0, len1, len2);
            }
            result[vIndex] = amix;
        }
    }

    return result;
}

void calMeshAttrMixedVertAreas(CMesh& mesh)
{
    vector<double> vVertAreas = computeMeshVertArea(mesh, VA_MIXED_VORONOI);
    mesh.addAttr<vector<double>>(vVertAreas, StrAttrVertMixedArea, AR_VERTEX, AT_VEC_DBL);
}

const std::vector<double>& getMeshVertMixedAreas(CMesh& mesh)
{
    if (!mesh.hasAttr(StrAttrVertMixedArea)) calMeshAttrMixedVertAreas(mesh);
    return mesh.getAttrValue<std::vector<double>>(StrAttrVertMixedArea);
}

ZGeom::ResultMeshMeanGaussCurvatures computeMeshMeanGaussCurvatures(CMesh &mesh)
{
    const double pi = ZGeom::PI;
    const int vertNum = mesh.vertCount();

    ResultMeshMeanGaussCurvatures result;
    result.mean_curvatures.resize(vertNum, 0);
    result.gauss_curvatures.resize(vertNum, 0);

    mesh.calAttrBoundaryVert();
    const vector<bool>& vertOnBoundary = mesh.getVertsOnBoundary();

    for (int vIndex = 0; vIndex < vertNum; ++vIndex)
    {
        const CVertex* vi = mesh.vert(vIndex);
        double angleSum = 0.0;		// sum of attaching corner's angle
        double amix = 0.0;
        ZGeom::Vec3d kh;

        // boundary vertex default to zero curvature     
        if (vertOnBoundary[vIndex]) continue;

        for (auto he = vi->m_HalfEdges.begin(); he != vi->m_HalfEdges.end(); ++he) {
            CHalfEdge* e0 = *he;
            CHalfEdge* e1 = e0->m_eNext;
            CHalfEdge* e2 = e1->m_eNext;
            double len0 = e0->length();
            double len1 = e1->length();
            double len2 = e2->length();

            // compute corner angle by cosine law 
            double corner = std::acos((len0*len0 + len2*len2 - len1*len1) / (2.0*len0*len2));
            angleSum += corner;
            double cota, cotc;
            amix += ZGeom::calMixedTriArea(len0, len1, len2, cota, cotc);

            const CVertex *pt1 = e1->vert(0), *pt2 = e1->vert(1);
            kh += (vi->pos() - pt1->pos()) * cota + (vi->pos() - pt2->pos()) * cotc;
        }

        kh /= 2.0 * amix;
        result.mean_curvatures[vIndex] = kh.length() / 2.0;	// half magnitude of kh
        result.gauss_curvatures[vIndex] = (2.0 * pi - angleSum) / amix;		// >0: convex; <0: concave
    }

    return result;
}

void calMeshAttrMeanGaussCurvatures(CMesh &mesh)
{
    auto result = computeMeshMeanGaussCurvatures(mesh);
    mesh.addAttr<vector<double>>(result.mean_curvatures, StrAttrVertMeanCurvatures, AR_VERTEX, AT_VEC_DBL);
    mesh.addAttr<vector<double>>(result.gauss_curvatures, StrAttrVertGaussCurvatures, AR_VERTEX, AT_VEC_DBL);
}

const std::vector<double>& getMeshMeanCurvatures(CMesh &mesh)
{
    if (!mesh.hasAttr(StrAttrVertMeanCurvatures)) calMeshAttrMeanGaussCurvatures(mesh);
    return mesh.getAttrValue<vector<double>>(StrAttrVertMeanCurvatures);
}

const std::vector<double>& getMeshGaussCurvatures(CMesh &mesh)
{
    if (!mesh.hasAttr(StrAttrVertGaussCurvatures)) calMeshAttrMeanGaussCurvatures(mesh);
    return mesh.getAttrValue<vector<double>>(StrAttrVertGaussCurvatures);
}

void gatherMeshStatistics(CMesh& mesh)
{
    // collect/compute statistics
    // 1. face normals
    // 2. vertex normals
    // 3. boundary vertex
    // 4. mesh center
    // 5. bounding box
    // 6. average edge length
    // 7. individual vertex curvature value

    mesh.calAttrFaceNormals();
    calMeshAttrVertNormals(mesh, ZGeom::VN_AREA_WEIGHT);
    mesh.calAttrBoundaryVert();
    mesh.calAttrBoundaryLoops();

    ZGeom::Vec3d center = mesh.calMeshCenter();
    ZGeom::Vec3d boundBox = mesh.calBoundingBox(center);
    double edgeLength = 0;
    for (int i = 0; i < mesh.halfEdgeCount(); ++i) {
        edgeLength += mesh.getHalfEdge(i)->length();
    }
    edgeLength /= mesh.halfEdgeCount();

    mesh.addAttr<double>(edgeLength, CMesh::StrAttrAvgEdgeLength, AR_UNIFORM, AT_DBL);
    mesh.addAttr<Vec3d>(center, CMesh::StrAttrMeshCenter, AR_UNIFORM, AT_VEC3);
    mesh.addAttr<Vec3d>(boundBox, CMesh::StrAttrMeshBBox, AR_UNIFORM, AT_VEC3);

    calMeshAttrMeanGaussCurvatures(mesh);
    calMeshAttrMixedVertAreas(mesh);
}

void calMeshAttrVertNormals(CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    vector<Vec3d> vertNormals = computeMeshVertNormals(mesh, vnc);
    mesh.addAttr<std::vector<ZGeom::Vec3d>>(vertNormals, CMesh::StrAttrVertNormal, AR_VERTEX, AT_VEC_VEC3);
}




}   // end of namespace