#include "geometry_processing.h"
#include <cstdlib>
#include <random>
#include <map>
#include <set>
#include <unordered_set>
#include <ppl.h>
#include <concurrent_vector.h>
#include "util.h"
#include "MatVecArithmetic.h"
#include "arithmetic.h"
#include "triBoxOverlap.h"


using namespace std;
using namespace concurrency;

namespace ZGeom {

std::function<double(double, double)> heat_gen_func = [](double lambda, double tau){
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
        vecNormals[i] = mesh.getFace(i)->calNormal();
    });

    return vecNormals;
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


vector<double> computeVertMeshDist(CMesh &mesh1, CMesh &mesh2)
{
    int nVerts1 = mesh1.vertCount();
    int nFaces2 = mesh2.faceCount();
    vector<vector<Vec3d>> mesh2Tri(nFaces2);
    for (int fIdx = 0; fIdx < nFaces2; ++fIdx)
        mesh2Tri[fIdx] = mesh2.getFace(fIdx)->getAllVertPos();

    vector<double> vVertMeshMinDist(nVerts1);
    concurrency::parallel_for(0, nVerts1, [&](int vIdx)
    {
        const Vec3d &vPos = mesh1.vertPos(vIdx);
        double minDistVi = 1e15;
        for (const vector<Vec3d>& tri : mesh2Tri)
            minDistVi = std::min(minDistVi, distPointTriangle(vPos, tri).distance);
        vVertMeshMinDist[vIdx] = minDistVi;
    });

    return vVertMeshMinDist;
}

double computeMeshMeanError(CMesh& mesh1, CMesh& mesh2)
{
    /* see "MESH - MEASURING ERRORS BETWEEN SURFACES USING THE HAUSDORFF (2002)" */
    int nVerts1 = mesh1.vertCount();
    vector<double> vVertAreas = getMeshVertMixedAreas(mesh1);
    vector<double> vVertMeshMinDist = computeVertMeshDist(mesh1, mesh2);

    double result(0);
    for (int vIdx = 0; vIdx < nVerts1; ++vIdx) result += vVertMeshMinDist[vIdx] * vVertAreas[vIdx];
    result /= std::accumulate(vVertAreas.begin(), vVertAreas.end(), (double)0);
    return result;
}

double computeSymMeshMeanError(CMesh &mesh1, CMesh &mesh2)
{
    return 0.5 * (computeMeshMeanError(mesh1, mesh2) + computeMeshMeanError(mesh2, mesh1));
}

double computeMeshRMSE(CMesh &mesh1, CMesh &mesh2)
{
    int nVerts1 = mesh1.vertCount();
    vector<double> vVertAreas = getMeshVertMixedAreas(mesh1);
    vector<double> vVertMeshMinDist = computeVertMeshDist(mesh1, mesh2);

    double result(0);
    for (int vIdx = 0; vIdx < nVerts1; ++vIdx) result += sqr(vVertMeshMinDist[vIdx]) * vVertAreas[vIdx];
    result = sqrt(result / std::accumulate(vVertAreas.begin(), vVertAreas.end(), (double)0));
    return result;
}

double computeSymMeshRMSE(CMesh &mesh1, CMesh &mesh2)
{
    return 0.5 * (computeMeshRMSE(mesh1, mesh2) + computeMeshRMSE(mesh2, mesh1));
}

double computeHausdorffDistance(CMesh &mesh1, CMesh &mesh2)
{
    vector<double> vVertMeshMinDist = computeVertMeshDist(mesh1, mesh2);
    return *std::max_element(vVertMeshMinDist.begin(), vVertMeshMinDist.end());
}

double computeSymHausdorffDistance(CMesh &mesh1, CMesh &mesh2)
{
    double dist12 = computeHausdorffDistance(mesh1, mesh2);
    double dist21 = computeHausdorffDistance(mesh2, mesh1);
    return std::max(dist12, dist21);
}


ZGeom::ResultDistPointPlane distPointPlane(Vec3d point, const Plane3& plane)
{
    ResultDistPointPlane result;
    result.signedDistance = dot(plane.normal, point) - plane.constant;
    result.distance = std::fabs(result.signedDistance);
    result.closestPoint = point - result.signedDistance * plane.normal;
    return result;
}

/* see "Distance between point and triangle in 3D 
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

/************************************************************************/
/* adapted from WildMagic                                               */
/* Much slower than distPointTriangle1                                  */
/************************************************************************/
ZGeom::ResultDistPointTriangle distPointTriangle2(Vec3d point, const std::vector<Vec3d>& triangle)
{
    vector<double> baryCoord(3, 0);

    /**** adapted from WildMagic ****/
    ZGeom::Vec3d V[3] = { triangle[0], triangle[1], triangle[2] };
    ZGeom::Vec3d p = point;

    ZGeom::Vec3d diff = V[0] - p;
    ZGeom::Vec3d edge0 = V[1] - V[0];
    ZGeom::Vec3d edge1 = V[2] - V[0];
    double   a = edge0.length2();
    double   b = edge0.dot(edge1);
    double   c = edge1.length2();
    double   d = diff.dot(edge0);
    double   e = diff.dot(edge1);
    double   f = diff.length2();
    double det = std::abs(a*c - b*b);
    double   s = b*e - c*d;
    double   t = b*d - a*e;
    //double s_bar = s / det, t_bar = t / det;
    double sqrDistance;

    if (s + t <= det)
    {
        if (s < 0.) {
            if (t < 0.) { // region 4
                if (d < 0.) {
                    if (-d >= a) {
                        sqrDistance = a + 2.*d + f;	// on V1
                        s = 1; t = 0;
                    }
                    else {
                        sqrDistance = f - d*d / a; // on E0
                        s = -d / a; t = 0;
                    }
                }
                else {
                    if (e >= 0.) {
                        sqrDistance = f;   // on V0
                        s = 0; t = 0;
                    }
                    else if (-e >= c) {
                        sqrDistance = c + 2.*e + f;	// on V2
                        s = 0; t = 1;
                    }
                    else {
                        sqrDistance = f - e*e / c;	//on E1
                        s = 0; t = -e / c;
                    }
                }
            }
            else {  // region 3
                if (e >= 0.) {
                    sqrDistance = f;	// on V0
                    s = 0; t = 0;
                }
                else if (-e >= c) {
                    sqrDistance = c + 2.*e + f;	// on V2
                    s = 0; t = 1;
                }
                else {
                    sqrDistance = f - e*e / c;	//on E1
                    s = 0; t = -e / c;
                }
            }
        }
        else if (t < 0.)  { // region 5
            if (d >= 0.) {
                sqrDistance = f;	// on V0
                s = 0; t = 0;
            }
            else if (-d >= a) {
                sqrDistance = a + 2.*d + f;	// on V1
                s = 1; t = 0;
            }
            else {
                sqrDistance = d*s + f - d*d / a;	// on E0
                s = -d / a; t = 0;
            }
        }
        else  { // region 0
            // The minimum is at an interior point of the triangle.
            double invDet = 1. / det;
            s *= invDet;
            t *= invDet;
            sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
        }
    } // if (s + t <= det)
    else
    {
        double tmp0, tmp1, numer, denom;

        if (s < 0.)  {// region 2
            tmp0 = b + d;
            tmp1 = c + e;
            if (tmp1 > tmp0) {
                numer = tmp1 - tmp0;
                denom = a - 2.*b + c;
                if (numer >= denom) {
                    sqrDistance = a + 2.*d + f;	// on V1?
                    s = 1; t = 0;
                }
                else {
                    s = numer / denom;
                    t = 1. - s;
                    sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
                }
            }
            else {
                if (tmp1 <= 0.) {
                    sqrDistance = c + 2.*e + f;	//on v2
                    s = 0; t = 1;
                }
                else if (e >= 0.) {
                    sqrDistance = f;	// on v0
                    s = 0; t = 0;
                }
                else {
                    sqrDistance = f - e*e / c;	// on E1?
                    s = 0; t = -e / c;
                }
            }
        }
        else if (t < 0.) { // region 6
            tmp0 = b + e;
            tmp1 = a + d;
            if (tmp1 > tmp0) {
                numer = tmp1 - tmp0;
                denom = a - 2.*b + c;
                if (numer >= denom) {
                    sqrDistance = c + 2.*e + f;	// on V2
                    s = 0.; t = 1.;
                }
                else {
                    t = numer / denom;
                    s = 1. - t;
                    sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
                }
            }
            else {
                if (tmp1 <= 0.) {
                    sqrDistance = a + 2.*d + f;	// on V1
                    s = 1.; t = 0.;
                }
                else if (d >= 0.) {
                    sqrDistance = f;	// on V0
                    s = 0.; t = 0.;
                }
                else {
                    sqrDistance = f - d*d / a;	// on E0
                    s = -d / a; t = 0;
                }
            }
        }
        else { // region 1
            numer = c + e - b - d;
            if (numer <= 0.) {
                sqrDistance = c + 2.*e + f;		// on V2
                s = 0.; t = 1.;
            }
            else {
                denom = a - 2.*b + c;
                if (numer >= denom) {
                    sqrDistance = a + 2.*d + f;	// on V1
                    s = 1.; t = 0.;
                }
                else {
                    s = numer / denom;
                    t = 1. - s;
                    sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
                }
            }
        }
    } // if (s+t > det)

    baryCoord[0] = 1.0 - s - t;
    baryCoord[1] = s;
    baryCoord[2] = t;

    if (baryCoord[0] < -1e-6 || baryCoord[1] < -1e-6 || baryCoord[2] < -1e-6) {
        cout << "Illegal barycentric coordinate: (" << baryCoord[0] << ", " << baryCoord[1] << ", " << baryCoord[2] << ")" << endl;
    }

    ResultDistPointTriangle result;
    result.distance = sqrt(fabs(sqrDistance));
    result.closestPoint = triangle[0] + s * edge0 + t * edge1;
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

std::vector<Vec3d> computeMeshVertNormals(const CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    int vertCount = mesh.vertCount(), faceCount = mesh.faceCount();
    const vector<Vec3d> faceNormals = computeMeshFaceNormals(mesh);
    vector<double> faceAreas(faceCount);
    for (int i = 0; i < faceCount; ++i) faceAreas[i] = mesh.getFace(i)->calArea();

    vector<Vec3d> vertNormals(vertCount);
    for (int vi = 0; vi < vertCount; ++vi) {
        vector<int> neighborFaceIdx = mesh.getVertAdjacentFaceIdx(vi, 1);
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

void calMeshAttrVertNormals(CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    vector<Vec3d> vertNormals = computeMeshVertNormals(mesh, vnc);
    mesh.addAttr<std::vector<ZGeom::Vec3d>>(vertNormals, CMesh::StrAttrVertNormals, AR_VERTEX, AT_VEC_VEC3);
}

const std::vector<Vec3d>& getMeshVertNormals(CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    if (!mesh.hasAttr(CMesh::StrAttrVertNormals)) calMeshAttrVertNormals(mesh, vnc);
    return mesh.getAttrValue<vector<Vec3d>>(CMesh::StrAttrVertNormals);
}

const std::vector<Vec3d>& getMeshVertNormals(const CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    return mesh.getAttrValue<vector<Vec3d>>(CMesh::StrAttrVertNormals);
}

ZGeom::VertCurvature calVertCurvature(const CMesh& mesh, int v_idx, bool compute_principal) {
    VertCurvature result;
    assert(v_idx >= 0 && v_idx < mesh.vertCount());

    // boundary vertex default to zero curvature   
    if (mesh.vert(v_idx)->judgeOnBoundary()) return result;

    const CVertex* vi = mesh.vert(v_idx);
    double angleSum = 0.0;		// sum of attaching corner's angle
    double amix = 0.0;
    ZGeom::Vec3d kh;

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
    result.mean_curv = kh.length() / 2.0;	            // half magnitude of kh
    result.gauss_curv = (2.0 * ZGeom::PI - angleSum) / amix; // >0: ellipse; <0: parabolic; =0: hyperbolic or plane
    
    // determine convex or concave
    const std::vector<Vec3d>& vert_normal = getMeshVertNormals(mesh);
    if (kh.dot(vert_normal[v_idx]) > 0) result.convexity = 1.;
    else result.convexity = -1.;
    
    if (compute_principal) {
        double delta = std::max<double>(0., sqr(result.mean_curv) - result.gauss_curv);
        delta = sqrt(delta);
        result.principal_curv_1 = result.convexity * result.mean_curv + delta;
        result.principal_curv_2 = result.convexity * result.mean_curv - delta;
    }

    return result;
}

void computeMeshCurvatures(CMesh& mesh, bool compute_principal)
{
    int vert_count = mesh.vertCount();
    vector<VertCurvature> all_curvatures(vert_count);
    parallel_for(0, vert_count, [&](int vi) {
        all_curvatures[vi] = calVertCurvature(mesh, vi, compute_principal);
    });

    vector<double> mean_curvs(vert_count), gauss_curvs(vert_count);
    for (int i = 0; i < vert_count; ++i) {
        mean_curvs[i] = all_curvatures[i].mean_curv;
        gauss_curvs[i] = all_curvatures[i].gauss_curv;
    }
    mesh.addAttr<vector<double>>(mean_curvs, StrAttrVertMeanCurvatures, AR_VERTEX, AT_VEC_DBL);
    mesh.addAttr<vector<double>>(gauss_curvs, StrAttrVertGaussCurvatures, AR_VERTEX, AT_VEC_DBL);

    if (compute_principal) {
        vector<double>  pri_curvs_1(vert_count), pri_curvs_2(vert_count);
        for (int i = 0; i < vert_count; ++i) {
            pri_curvs_1[i] = all_curvatures[i].principal_curv_1;
            pri_curvs_2[i] = all_curvatures[i].principal_curv_2;
        }
        mesh.addAttr<vector<double>>(pri_curvs_1, StrAttrVertPrincipalCurvatures1, AR_VERTEX, AT_VEC_DBL);
        mesh.addAttr<vector<double>>(pri_curvs_2, StrAttrVertPrincipalCurvatures2, AR_VERTEX, AT_VEC_DBL);
    }

    mesh.addAttr<vector<VertCurvature>>(all_curvatures, StrAttrVertAllCurvatures, AR_VERTEX, AT_UNKNOWN);
}

void computeShapeIndex(CMesh& mesh)
{
    int vert_count = mesh.vertCount();
    vector<double> result(vert_count, 0);
    if (!mesh.hasAttr(StrAttrVertAllCurvatures)) computeMeshCurvatures(mesh, true);
    const auto& curvatures = mesh.getAttrValue<vector<ZGeom::VertCurvature>>(StrAttrVertAllCurvatures);
    for (int i = 0; i < vert_count; ++i) {
        const VertCurvature& elem = curvatures[i];
        double k1 = elem.principal_curv_1, k2 = elem.principal_curv_2;
        if (fabs(elem.mean_curv) < 1e-5) result[i] = 0;
        else if (fabs(k1 - k2) < 1e-5) {
            result[i] = (k1 >= 0 ? 1 : -1);
        }
        else {
            result[i] = 2 / ZGeom::PI * atan((k1 + k2)/(k1 - k2));
        }
    }

    mesh.addAttr<vector<double>>(result, StrAttrVertShapeIndex, AR_VERTEX, AT_VEC_DBL);
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
    // 7. individual vertex curvatures
    // 8. vertex voronoi area

    mesh.calAttrVertNormals();

    mesh.calAttrBoundaryVert();
    identifyMeshBoundaries(mesh);
    
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

    computeMeshCurvatures(mesh, true);  // mean, gauss, principals
    calMeshAttrMixedVertAreas(mesh);
}

/********************************************************************************/
/* Adapted from AABB-triangle overlap test code                                 */
/* by Tomas Akenine-M�ller                                                      */
/********************************************************************************/
bool testTriBoxOverlap(const std::vector<Vec3d>& triangle, Vec3d boxCenter, Vec3d boxHalfsize)
{
    float boxcenter[3] = { (float)boxCenter[0], (float)boxCenter[1], (float)boxCenter[2] };
    float boxhalfsize[3] = { (float)boxHalfsize[0], (float)boxHalfsize[1], (float)boxHalfsize[2] };
    float triverts[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            triverts[i][j] = (float)triangle[i][j];
    return 1 == triBoxOverlap(boxcenter, boxhalfsize, triverts);
}

int triObtuseEdge(const std::vector<Vec3d>& triVerts)
{
    assert(triVerts.size() == 3);
    double sqlen0 = (triVerts[1] - triVerts[2]).length2();
    double sqlen1 = (triVerts[0] - triVerts[2]).length2();
    double sqlen2 = (triVerts[0] - triVerts[1]).length2();
    if (sqlen0 > sqlen1 + sqlen2) return 0;
    else if (sqlen1 > sqlen0 + sqlen2) return 1;
    else if (sqlen2 > sqlen0 + sqlen1) return 2;
    else return -1; // non-obtuse triangle
}

int triObtuseEdge(const std::vector<double>& triLengths)
{
    if (sqr(triLengths[0]) > sqr(triLengths[1]) + sqr(triLengths[2])) return 0;
    else if (sqr(triLengths[1]) > sqr(triLengths[0]) + sqr(triLengths[2])) return 1;
    else if (sqr(triLengths[2]) > sqr(triLengths[0]) + sqr(triLengths[1])) return 2;
    else return -1; // non-obtuse triangle
}

std::vector<MeshRegion> identifyMeshBoundaries(CMesh& mesh)
{
    typedef vector<vector<int>> VecVecInt;

    int boundaryNum = 0;
    int vertCount = mesh.vertCount();
    const vector<bool>& vVertOnBoundary = mesh.getVertsOnBoundary();

    VecVecInt boundaryEdges;
    vector<bool> vertVisited(vertCount, false);
    for (int i = 0; i < vertCount; i++) {
        // find boundary loop from boundary vertex i if it is not in any loop 
        if (!vVertOnBoundary[i] || vertVisited[i]) continue;
        vector<int> edgeLoop;
        int currentIndex = i;
        do {
            vertVisited[currentIndex] = true;
            int edgeIndex = -1;
            for (CHalfEdge* he : mesh.vert(currentIndex)->getHalfEdges()) {
                if (he->isBoundaryEdge()) {
                    edgeIndex = he->getIndex(); break;
                }
            }
            if (edgeIndex == -1) {
                cerr << "Boundary ended without loop!";
                break;
            }
            currentIndex = mesh.m_vHalfEdges[edgeIndex]->getVertIndex(1);
            edgeLoop.push_back(edgeIndex);
        } while (currentIndex != i);

        boundaryEdges.push_back(edgeLoop);
    }

    // sort by number of edges of boundaries
    std::sort(boundaryEdges.begin(), boundaryEdges.end(),
        [](const vector<int>& v1, const vector<int>& v2) { return v1.size() < v2.size(); });

    vector<MeshRegion> result(boundaryEdges.size());
    for (int i = 0; i < (int)boundaryEdges.size(); ++i) {
        result[i].he_on_boundary = boundaryEdges[i];
        for (int heIdx : boundaryEdges[i])
            result[i].vert_on_boundary.push_back(mesh.getHalfEdge(heIdx)->getVertIndex(0));
    }

    mesh.addAttr<vector<MeshRegion>>(result, StrAttrMeshHoleRegions, AR_UNIFORM);
    return result;
}

const std::vector<MeshRegion>& getMeshBoundaryLoops(CMesh &mesh)
{
    if (!mesh.hasAttr(StrAttrMeshHoleRegions)) identifyMeshBoundaries(mesh);
    return mesh.getAttrValue<vector<MeshRegion>>(StrAttrMeshHoleRegions);
}

std::vector<std::vector<int>> getMeshBoundaryLoopVerts(CMesh &mesh)
{
    const std::vector<MeshRegion>& vBoundaries = getMeshBoundaryLoops(mesh);
    vector<vector<int>> result;
    for (const MeshRegion& hb : vBoundaries)
        result.push_back(hb.vert_on_boundary);
    return result;
}

std::vector<std::vector<int>> getMeshBoundaryLoopHalfEdges(CMesh &mesh)
{
    const std::vector<MeshRegion>& vBoundaries = getMeshBoundaryLoops(mesh);
    vector<vector<int>> result;
    for (const MeshRegion& hb : vBoundaries)
        result.push_back(hb.he_on_boundary);
    return result;
}

int calMeshGenus(CMesh &mesh)
{
    int b = (int)getMeshBoundaryLoops(mesh).size();
    int euler_number = mesh.calEulerNum();
    return (2 - euler_number - b) / 2;
}

std::vector<bool> getMeshVertsOnHoles(CMesh &mesh)
{

    const int MAX_HOLE_SIZE = 200;
    const vector<vector<int>>& boundariesVert = getMeshBoundaryLoopVerts(mesh);
    vector<bool> result(mesh.vertCount(), false);
    for (int i = 0; i < boundariesVert.size(); ++i) {
        const vector<int>& vVertIdx = boundariesVert[i];
        if (vVertIdx.size() > MAX_HOLE_SIZE) continue;
        for (int vIdx : vVertIdx) result[vIdx] = true;
    }
    return result;
}

std::unique_ptr<CMesh> cutFromMesh(CMesh &oldMesh, const std::vector<int>& cutFaceIdx)
{
    std::unique_ptr<CMesh> newMesh = std::make_unique<CMesh>();
    
    vector<int> vRemainingFaceIdx;
    std::unordered_set<int> cutFaces(cutFaceIdx.begin(), cutFaceIdx.end());
    for (int fIdx = 0; fIdx < oldMesh.faceCount(); ++fIdx) {
        if (cutFaces.find(fIdx) == cutFaces.end())
            vRemainingFaceIdx.push_back(fIdx);
    }

    oldMesh.getSubMeshFromFaces(vRemainingFaceIdx, oldMesh.getMeshName() + "_cut", *newMesh);
    gatherMeshStatistics(*newMesh);
    return newMesh;
}


std::unique_ptr<CMesh> cutMeshTo(CMesh &oldMesh, const std::vector<int>& cutFaceIdx)
{
    unique_ptr<CMesh> newMesh = std::make_unique<CMesh>();
    oldMesh.getSubMeshFromFaces(cutFaceIdx, oldMesh.getMeshName() + "_fragment", *newMesh);
    gatherMeshStatistics(*newMesh);
    return newMesh;
}

/* "Filling holes in meshes" (SGP 2003) */
void triangulateMeshHoles(CMesh &mesh)
{ 
    vector<MeshRegion> vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);

    for (int holeIdx = 0; holeIdx < (int)vHoles.size(); ++holeIdx)
    {
        int nOldVerts = mesh.vertCount();
        int nOldEdges = mesh.halfEdgeCount();
        int nOldFaces = mesh.faceCount();
        const vector<int>& boundaryEdgeIdx = vHoles[holeIdx].he_on_boundary;

        const int N = (int)boundaryEdgeIdx.size();    // number of boundary vertices
        vector<int> boundaryVertIdx(N);
        vector<CVertex*> boundaryVertPtr(N);
        vector<Vec3d> boundaryVertPos(N);
        for (int i = 0; i < N; ++i) {
            boundaryVertIdx[i] = mesh.getHalfEdge(boundaryEdgeIdx[i])->getVertIndex(0);
            boundaryVertPtr[i] = mesh.vert(boundaryVertIdx[i]);
            boundaryVertPos[i] = boundaryVertPtr[i]->pos();
        }

        vector<vector<WeightSet>> weights(N);
        vector<vector<int>> midIdx(N);
        for (int i = 0; i < N; ++i) {
            weights[i].resize(N);
            midIdx[i].resize(N, 0);
        }
        for (int i = 0; i < N - 2; ++i) {
            int vIdx1 = boundaryVertIdx[i], vIdx2 = boundaryVertIdx[i + 1], vIdx3 = boundaryVertIdx[i + 2];
            CVertex *v1 = mesh.vert(vIdx1), *v2 = mesh.vert(vIdx2), *v3 = mesh.vert(vIdx3);
            CHalfEdge *e12 = v1->adjacentTo(v2), *e23 = v2->adjacentTo(v3);
            Vec3d vn = triNormal(v1->pos(), v3->pos(), v2->pos());
            Vec3d vn12 = e12->getAttachedFace()->calNormal(), vn23 = e23->getAttachedFace()->calNormal();
            weights[i][i + 2].angle = max(vecAngle(vn, vn12), vecAngle(vn, vn23));
            weights[i][i + 2].area = triArea(boundaryVertPos[i], boundaryVertPos[i + 1], boundaryVertPos[i + 2]);
            midIdx[i][i + 2] = i + 1;
        }

        for (int j = 3; j < N; ++j) {
            for (int i = 0; i < N - j; ++i) {
                int k = i + j;
                CVertex *vi = boundaryVertPtr[i], *vk = boundaryVertPtr[k];
                int minMidIdx = -1;
                WeightSet minWeight_ik(4, 1e30);

                for (int m = i + 1; m < k; ++m) {
                    CVertex *vm = boundaryVertPtr[m];
                    WeightSet triWeight;
                    triWeight.area = ZGeom::triArea(boundaryVertPos[i], boundaryVertPos[m], boundaryVertPos[k]);

                    Vec3d vn_ikm = triNormal(boundaryVertPos[i], boundaryVertPos[k], boundaryVertPos[m]);
                    Vec3d vn_im, vn_mk;
                    if (m == i + 1)
                        vn_im = vi->adjacentTo(vm)->getAttachedFace()->calNormal();
                    else {
                        int l1 = midIdx[i][m];
                        vn_im = triNormal(boundaryVertPos[i], boundaryVertPos[m], boundaryVertPos[l1]);
                    }
                    if (k == m + 1)
                        vn_mk = vm->adjacentTo(vk)->getAttachedFace()->calNormal();
                    else {
                        int l2 = midIdx[m][k];
                        vn_mk = triNormal(boundaryVertPos[m], boundaryVertPos[k], boundaryVertPos[l2]);
                    }
                    triWeight.angle = max(vecAngle(vn_ikm, vn_im), vecAngle(vn_ikm, vn_mk));
                    if (k == i + N - 1) {
                        Vec3d vn_ik = vk->adjacentTo(vi)->getAttachedFace()->calNormal();
                        triWeight.angle = max(triWeight.angle, vecAngle(vn_ikm, vn_ik));
                    }

                    WeightSet weight_imk = weights[i][m] + weights[m][k] + triWeight;
                    if (weight_imk < minWeight_ik) {
                        minWeight_ik = weight_imk;
                        minMidIdx = m;
                    }
                }
                midIdx[i][k] = minMidIdx;
                weights[i][k] = minWeight_ik;
            }
        }

        MeshLineList triangulationLines;
        vector<vector<int>> patchTri;
        queue<pair<int, int>> tracePairs;
        tracePairs.push(make_pair(0, N - 1));
        while (!tracePairs.empty()) {
            pair<int, int> lineIdx = tracePairs.front();
            tracePairs.pop();
            int i = lineIdx.first, k = lineIdx.second;
            if (i + 2 == k) {
                patchTri.push_back(vector < int > { i, i + 1, k });
                triangulationLines.push_back(LineSegment(boundaryVertPos[i], boundaryVertPos[k]));
            }
            else {
                int o = midIdx[i][k];
                if (o > i + 1) tracePairs.push(make_pair(i, o));
                patchTri.push_back(vector < int > { i, o, k });
                triangulationLines.push_back(LineSegment(boundaryVertPos[i], boundaryVertPos[k]));
                if (o < k - 1) tracePairs.push(make_pair(o, k));
            }
        }

        std::cout << "boundary vert count: " << boundaryVertIdx.size() << endl;
        std::cout << "patching face count: " << patchTri.size() << endl;

        vector<CHalfEdge*> patchEdges;
        vector<CFace*> patchFaces;
        for (vector<int> tri : patchTri) {
            CVertex *vi = mesh.vert(boundaryVertIdx[tri[0]]),
                *vj = mesh.vert(boundaryVertIdx[tri[1]]),
                *vk = mesh.vert(boundaryVertIdx[tri[2]]);
            // construct new face (vi->vk->vj->vi)
            CFace *f = new CFace(3);
            CHalfEdge *eik = new CHalfEdge(), *ekj = new CHalfEdge(), *eji = new CHalfEdge();   // be careful with the clockwise
            eik->setVertOrigin(vi);
            ekj->setVertOrigin(vk);
            eji->setVertOrigin(vj);            
            vi->addHalfEdge(eik); vj->addHalfEdge(eji); vk->addHalfEdge(ekj);
            CMesh::makeFace(eik, ekj, eji, f);
            patchEdges.push_back(eik); patchEdges.push_back(ekj); patchEdges.push_back(eji);
            patchFaces.push_back(f);
        }
        for (CHalfEdge* he : patchEdges) {
            if (he->twinHalfEdge() == NULL) {
                CVertex *v1 = he->vert(0), *v2 = he->vert(1);
                CHalfEdge* e21 = v2->adjacentTo(v1);
                CHalfEdge::makeTwins(he, e21);
            }
        }
        for (int i = 0; i < (int)patchEdges.size(); ++i) {
            patchEdges[i]->m_eIndex = nOldEdges + i;
            mesh.m_vHalfEdges.push_back(patchEdges[i]);
        }
        for (int i = 0; i < (int)patchFaces.size(); ++i) {
            patchFaces[i]->m_fIndex = nOldFaces + i;
            mesh.m_vFaces.push_back(patchFaces[i]);
        }

        vHoles[holeIdx].face_inside.clear();
        for (CFace *f : patchFaces) 
            vHoles[holeIdx].face_inside.push_back(f->getFaceIndex());
    } // for each hole

    for (MeshRegion& hole : vHoles) hole.determineBoundaryHalfEdges(mesh);
    mesh.clearNonEssentialAttributes();
    gatherMeshStatistics(mesh);
    mesh.addAttr<vector<MeshRegion>>(vHoles, StrAttrMeshHoleRegions, AR_UNIFORM, AT_UNKNOWN);    
}

/* "Filling holes in meshes" (SGP 2003) */
/* "A multistep approach to restoration of locally undersampled meshes" (GMP 2008) */
void refineMeshHoles(CMesh &mesh, double lambda /*= 0.5*/)
{
    vector<MeshRegion> vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) return;

    /* Perform Delaunay-like refinement */
    for (int holeIdx = 0; holeIdx < (int)vHoles.size(); ++holeIdx)
    {
        MeshRegion& cur_hole = vHoles[holeIdx];
        const double preferred_edge_length = ZGeom::estimateHoleEdgeLength(mesh, cur_hole, 3);

        set<CVertex*> verts_in_hole;
        set<CFace*> faces_in_hole;
        for (int fi : cur_hole.face_inside) faces_in_hole.insert(mesh.getFace(fi));
        cur_hole.determineBoundaryHalfEdges(mesh);
        
        set<CHalfEdge*> boundary_he;
        for (int he_idx : cur_hole.he_on_boundary)
            boundary_he.insert(mesh.getHalfEdge(he_idx));

        while (true)
        {
            // 1. for each triangle inside hole, check if any face is splittable
            bool new_vert_added(false);
            while (true)
            {
                /* keep add new_vert until no new_vert can be selected */
                CVertex* new_vert = nullptr;
                for (CFace* face : faces_in_hole)
                {
                    CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                    int obtuse_edge = ZGeom::triObtuseEdge(face->getAllEdgeLengths());
                    if (obtuse_edge == -1)  // non-obtuse triangle: 3-split face
                    {
                        bool splitTest = true;
#if 0
                        for (int m = 0; m < 3; ++m) {
                            if ((face->vert(m)->pos() - vc).length() <= lambda * preferred_edge_length) {
                                splitTest = false; break;
                            }
                        }
#else
                        Vec3d vc = face->calBarycenter();
                        double totalLength(0);
                        for (int m = 0; m < 3; ++m)
                            totalLength += (face->vert(m)->pos() - vc).length();
                        splitTest = (totalLength >= 3 * lambda * preferred_edge_length);                            
#endif
                        if (splitTest) 
                        {
                            new_vert = mesh.faceSplit3(face);
                            // relax edges opposite to the new added centroid
                            for (int k = 0; k < 3; ++k) {
                                if (!setHas(boundary_he, fe[k]))
                                    mesh.relaxEdge(fe[k]);
                            }

                            break;
                        }
                    }
                    else // obtuse triangle: 2-split obtuse edge
                    {
                        CHalfEdge *obtuse_he = face->getHalfEdge(obtuse_edge);
                        if (setHas(boundary_he, obtuse_he))
                            continue; // don't split boundary edges
                        if (obtuse_he->twinHalfEdge() && !setHas(faces_in_hole, obtuse_he->twinHalfEdge()->getAttachedFace()))
                        {
                            std::cout << "Should not happen!" << std::endl;
                            continue;
                        }

                        bool splitTest = true;
#if 0
                        splitTest = (obtuse_he->length() >= 2 * lambda * preferred_edge_length);
#else
                        Vec3d vc = obtuse_he->midEdge();
                        double totalLength(0);
                        for (int m = 0; m < 3; ++m)
                            totalLength += (face->vert(m)->pos() - vc).length();
                        splitTest = (totalLength >= 3 * lambda * preferred_edge_length);
#endif                        
                        if (splitTest) 
                        {
                            vector<CHalfEdge*> diamond_edges = { obtuse_he->nextHalfEdge(), obtuse_he->prevHalfEdge() };
                            if (obtuse_he->twinHalfEdge() != nullptr) {
                                CHalfEdge* obtuse_he2 = obtuse_he->twinHalfEdge();
                                diamond_edges.push_back(obtuse_he2->nextHalfEdge());
                                diamond_edges.push_back(obtuse_he2->prevHalfEdge());
                            }
                            new_vert = mesh.edgeSplit(obtuse_he);

                            // relax edges opposite to the new added edge_center
                            for (int k = 0; k < (int)diamond_edges.size(); ++k) {
                                if (!setHas(boundary_he, diamond_edges[k]))
                                    mesh.relaxEdge(diamond_edges[k]);
                            }

                            break;
                        }
                    }
                }

                if (new_vert != nullptr) {
                    verts_in_hole.insert(new_vert);
                    for (CFace* newF : new_vert->getAdjacentFaces())
                        faces_in_hole.insert(newF);
                    if(verts_in_hole.size() % 10 == 0) 
                        std::cout << "#verts_in_hole: " << verts_in_hole.size() << '\n';
                    new_vert_added = true;
                }
                else break;
            }
                            
            // 2. If no new points were inserted in step 1, hole filling process completes
            if (!new_vert_added) break;  

            // 3. relax all edges in patch mesh
            vector<CHalfEdge*> inside_he;
            for (CVertex* pv : verts_in_hole) {
                for (CHalfEdge* he : pv->getHalfEdges())
                    inside_he.push_back(he);
            }               

            int swapIterCount = 0;
            while (true) {
                bool noEdgeSwap = true;
                for (int i = 0; i < (int)inside_he.size(); ++i) {
                    if (mesh.relaxEdge(inside_he[i])) {
                        noEdgeSwap = false;
                    }
                }
                if (noEdgeSwap) break;
                if (swapIterCount++ >= 100) break;
            }
            if (swapIterCount >= 100) cout << "Maximum swap reached!\n";

            // 4. go to step 1
        }

        cur_hole.face_inside.clear();
        cur_hole.vert_inside.clear();
        for (CFace* pf : faces_in_hole) cur_hole.face_inside.push_back(pf->getFaceIndex());
        for (CVertex* pv : verts_in_hole) cur_hole.vert_inside.push_back(pv->getIndex());

        std::cout << "[Hole #" << holeIdx
            << "] #vert: " << cur_hole.vert_inside.size() 
            << "; Est. edge length: " << cur_hole.adjacent_edge_length
            << "; Actual edge length: " << calMeshRegionAvgEdgeLen(mesh, cur_hole) 
            << std::endl;

    } // for each hole

    for (MeshRegion& hole : vHoles) hole.determineBoundaryHalfEdges(mesh);

    mesh.clearNonEssentialAttributes();
    gatherMeshStatistics(mesh);
    mesh.addAttr<vector<MeshRegion>>(vHoles, StrAttrMeshHoleRegions, AR_UNIFORM, AT_UNKNOWN);
}

void refineMeshHoles2(CMesh &mesh, double lambda /*= 0.7*/)
{
    vector<MeshRegion> vHoles = getMeshBoundaryLoops(mesh);
    if (vHoles.empty()) return;

    /* Perform Delaunay-like refinement */
    for (int holeIdx = 0; holeIdx < (int)vHoles.size(); ++holeIdx)
    {
        MeshRegion& cur_hole = vHoles[holeIdx];
        const double preferred_edge_length = ZGeom::estimateHoleEdgeLength(mesh, cur_hole, 3);

        set<CVertex*> verts_in_hole;
        set<CFace*> faces_in_hole;
        for (int fi : cur_hole.face_inside) faces_in_hole.insert(mesh.getFace(fi));
        set<CHalfEdge*> boundary_he;
        for (CFace* f : faces_in_hole) {
            for (CHalfEdge* he : f->getAllHalfEdges()) {
                if (he->twinHalfEdge() == nullptr || !setHas(faces_in_hole, he->twinHalfEdge()->getAttachedFace()))
                    boundary_he.insert(he);
            }
        }

        while (true)
        {
            // 1. for each triangle inside hole, check if any face is splittable
            bool new_vert_added(false);
            CVertex* new_vert = nullptr;
            CFace* face_to_split = nullptr;
            double max_len(0);
            for (CFace* face : faces_in_hole)
            {
                bool splitTest = true;
                double center_len(0);
                Vec3d vc = face->calBarycenter();
                for (int m = 0; m < 3; ++m) {
                    double dist = (face->vert(m)->pos() - vc).length();
                    if (dist <= lambda * preferred_edge_length) {
                        splitTest = false;
                        break;
                    }
                    center_len += dist;
                }

                if (splitTest) {
                    if (center_len > max_len) {
                        max_len = center_len;
                        face_to_split = face;
                    }
                }
            }

            if (face_to_split != nullptr)
            {
                vector<CHalfEdge*> fe = face_to_split->getAllHalfEdges();
                new_vert = mesh.faceSplit3(face_to_split);
                
                // relax edges opposite to the new added centroid
                for (int k = 0; k < 3; ++k) {
                    if (!setHas(boundary_he, fe[k]))
                        mesh.relaxEdge(fe[k]);
                }

                verts_in_hole.insert(new_vert);
                for (CFace* newF : new_vert->getAdjacentFaces())
                    faces_in_hole.insert(newF);
            }
            // 2. If no new points were inserted in step 1, hole filling process completes
            else break;

            // 3. relax all edges in patch mesh
            vector<CHalfEdge*> inside_he;
            for (CVertex* pv : verts_in_hole) {
                for (CHalfEdge* he : pv->getHalfEdges())
                    inside_he.push_back(he);
            }

            int swapIterCount = 0;
            while (true) {
                bool noEdgeSwap = true;
                for (int i = 0; i < (int)inside_he.size(); ++i) {
                    if (mesh.relaxEdge(inside_he[i])) {
                        noEdgeSwap = false;
                    }
                }
                if (noEdgeSwap) break;
                if (swapIterCount++ >= 100) break;
            }
            if (swapIterCount >= 100) cout << "!!Maximum swap reached!\n";

            // 4. go to step 1
        }

        cur_hole.face_inside.clear();
        cur_hole.vert_inside.clear();
        for (CFace* pf : faces_in_hole) cur_hole.face_inside.push_back(pf->getFaceIndex());
        for (CVertex* pv : verts_in_hole) cur_hole.vert_inside.push_back(pv->getIndex());

        std::cout << "[Hole #" << holeIdx
            << "] #vert: " << cur_hole.vert_inside.size()
            << "; Est. edge length: " << cur_hole.adjacent_edge_length
            << "; Actual edge length: " << calMeshRegionAvgEdgeLen(mesh, cur_hole)
            << std::endl;

    } // for each hole

    for (MeshRegion& hole : vHoles) hole.determineBoundaryHalfEdges(mesh);

    mesh.clearNonEssentialAttributes();
    gatherMeshStatistics(mesh);
    mesh.addAttr<vector<MeshRegion>>(vHoles, StrAttrMeshHoleRegions, AR_UNIFORM, AT_UNKNOWN);
}

double estimateHoleEdgeLength(CMesh& mesh, MeshRegion& hole, int ring /*= 1*/)
{
    set<int> inside_he;
    for (int fi : hole.face_inside) {
        for (CHalfEdge* he : mesh.getFace(fi)->getAllHalfEdges()) {
            inside_he.insert(he->getIndex());
        }
    }

    set<int> surrounding_verts(hole.vert_on_boundary.begin(), hole.vert_on_boundary.end());    
    if (ring > 1) {
        vector<int> hole_verts(hole.vert_inside);
        for (int vi : hole.vert_on_boundary) hole_verts.push_back(vi);

        vector<int> neighbor_verts = vertSurroundingVerts(mesh, hole_verts, ring - 1);
        for (int vi : neighbor_verts)
            surrounding_verts.insert(vi);
    }

    set<int> neighboring_he;
    for (int vi : surrounding_verts) {
        for (CHalfEdge* he : mesh.vert(vi)->getHalfEdges()) {
            int he_idx = he->getIndex();
            if (!setHas(inside_he, he_idx)) neighboring_he.insert(he_idx);
        }
            
    }

    double lengthSum(0);
    int count(0);
    for (int eIdx : neighboring_he) {
        CHalfEdge *he = mesh.getHalfEdge(eIdx);
        lengthSum += he->length();
        count++;       
    }
 
    double estimated_avg_len = lengthSum / (double)count;
    hole.adjacent_edge_length = estimated_avg_len;
    return estimated_avg_len;
}

double calMeshRegionAvgEdgeLen(CMesh& mesh, MeshRegion& hole)
{
    if (hole.face_inside.empty()) return 0;
    double result(0);
    for (int fi : hole.face_inside) {
        for (CHalfEdge* he : mesh.getFace(fi)->getAllHalfEdges())
            result += he->length();
    }
    return result / double(hole.face_inside.size() * 3);
}

double distSubMesh(CMesh &mesh1, const std::vector<int>& faces1, CMesh &mesh2, const std::vector<int>& faces2, MeshDistMeasure measure /*= RMSE*/)
{
    double result;
    unique_ptr<CMesh> submesh1 = std::move(cutMeshTo(mesh1, faces1)), submesh2 = std::move(cutMeshTo(mesh2, faces2));
    switch (measure)
    {
    case ZGeom::MEAN_ERROR: 
        result = computeMeshMeanError(*submesh1, *submesh2);
        break;
    case ZGeom::SYM_MEAN_ERROR: 
        result = computeSymMeshMeanError(*submesh1, *submesh2);
        break;
    case ZGeom::RMSE:
        result = computeMeshRMSE(*submesh1, *submesh2);
        break;
    case ZGeom::SYM_RMSE:
        result = computeSymMeshRMSE(*submesh1, *submesh2);
        break;
    case ZGeom::HAUSDORFF:
        result = computeHausdorffDistance(*submesh1, *submesh2);
        break;
    case ZGeom::SYM_HAUSDORFF:
        result = computeSymHausdorffDistance(*submesh1, *submesh2);
        break;
    default:
        break;
    }
    return result;
}

std::vector<MeshRegion*> getMeshHoleRegions(CMesh& mesh)
{
    vector<MeshRegion*> result;

    if (mesh.hasAttr(StrAttrMeshHoleRegions)) {
        auto& vHoles = mesh.getAttrValue<vector<MeshRegion>>(StrAttrMeshHoleRegions);
        for (MeshRegion& mr : vHoles)
            result.push_back(&mr);
    }
    if (mesh.hasAttr(StrAttrManualHoleRegions)) {
        auto& vManualHoles = mesh.getAttrValue<vector<MeshRegion>>(StrAttrManualHoleRegions);
        for (MeshRegion& mr : vManualHoles)
            result.push_back(&mr);
    }

    return result;
}

std::vector<int> getMeshRegionsInsideVerts(const std::vector<MeshRegion>& vRegions)
{
    set<int> verts;
    for (auto& mr : vRegions)
        for (int vi : mr.vert_inside)
            verts.insert(vi);

    return vector<int>{verts.begin(), verts.end()};
}

std::vector<int> getMeshRegionsBoundaryVerts(const std::vector<MeshRegion>& vRegions)
{
    set<int> verts;
    for (auto& mr : vRegions)
        for (int vi : mr.vert_on_boundary)
            verts.insert(vi);

    return vector < int > {verts.begin(), verts.end()};
}

void mergeMeshRegions(CMesh& mesh, std::vector<MeshRegion>& vRegions)
{
    vector<set<int>> inside_verts;
    for (MeshRegion& mr : vRegions) 
        inside_verts.push_back(set<int>(mr.vert_inside.begin(), mr.vert_inside.end()));
    
    while(true) 
    {
        int r1(-1), r2(-1);
        for (int i = 0; i < (int)inside_verts.size() - 1; ++i) {            
            set<int> s1 = inside_verts[i];
            vector<int> ring_verts = vertSurroundingVerts(mesh, vector < int > {s1.begin(), s1.end()}, 1);
            s1 = setCombine(s1, set < int > {ring_verts.begin(), ring_verts.end()});
            for (int j = i+1; j < (int)inside_verts.size(); ++j) {
                const set<int> &s2 = inside_verts[j];
                if (setOverlap(s1, s2)) {
                    r1 = i; r2 = j; break;
                }
            }
            if (r1 != -1) break;
        }

        if (r1 != -1) {
            inside_verts[r1] = setCombine(inside_verts[r1], inside_verts[r2]);
            inside_verts.erase(inside_verts.begin() + r2);
        }
        else break;
    } 
    
    vRegions.clear();
    for (set<int>& region_vert : inside_verts) {
        vRegions.push_back(meshRegionFromInsideVerts(mesh, vector<int>(region_vert.begin(), region_vert.end())));
    }
}

ZGeom::MeshRegion generateRingMeshRegion(const CMesh& mesh, int seedVert, int ring)
{
    vector<int> vert_inside = mesh.getVertNeighborVerts(seedVert, ring - 1, true);
    vector<int> vert_on_boundary = mesh.getVertIsoNeighborVerts(seedVert, ring);
    set<int> face_inside;
    for (int vi : vert_inside) {
        for (CFace* f : mesh.vert(vi)->getAdjacentFaces())
            face_inside.insert(f->getFaceIndex());
    }

    MeshRegion result;
    result.face_inside = vector<int>(face_inside.begin(), face_inside.end());
    result.vert_inside = vert_inside;
    result.vert_on_boundary = vert_on_boundary;
    result.determineBoundaryHalfEdges(mesh);
    return result;
}


MeshRegion generateRandomMeshRegion(const CMesh& mesh, int seedVert, int holeSize)
{
    MeshRegion result;

    set<int> vertInHole;
    set<int> faceInHole;
    set<int> boundaryVerts;
    default_random_engine generator((unsigned int)time(NULL));

    vertInHole.insert(seedVert);
    for (auto f : mesh.vert(seedVert)->getAdjacentFaces()) {
        faceInHole.insert(f->getFaceIndex());
        for (int k = 0; k < 3; ++k) {
            int faceVertIdx = f->vertIdx(k);
            if (vertInHole.find(faceVertIdx) == vertInHole.end())
                boundaryVerts.insert(faceVertIdx);
        }
    }

    while (vertInHole.size() < holeSize) {
        uniform_int_distribution<int> distr1(0, (int)boundaryVerts.size() - 1);
        auto iter_bv = boundaryVerts.begin();
        std::advance(iter_bv, distr1(generator));             
        int newHoleVert = *iter_bv;
        vertInHole.insert(newHoleVert);
        for (auto f : mesh.vert(newHoleVert)->getAdjacentFaces()) {
            if (faceInHole.find(f->getFaceIndex()) != faceInHole.end())
                continue;       // face already considered     
            faceInHole.insert(f->getFaceIndex());
            for (int k = 0; k < 3; ++k) {
                int faceVertIdx = f->vertIdx(k);
                if (vertInHole.find(faceVertIdx) == vertInHole.end())
                    boundaryVerts.insert(faceVertIdx);
            }
        }
        boundaryVerts.erase(newHoleVert);
        for (auto iter = boundaryVerts.begin(); iter != boundaryVerts.end();) {
            int vIdx = *iter;
            bool encompassed = true;
            for (auto f : mesh.vert(vIdx)->getAdjacentFaces()) {
                if (faceInHole.find(f->getFaceIndex()) == faceInHole.end()) {
                    // found a face not in hole, so vIdx is not in hole yet
                    encompassed = false;
                    break;
                }
            }
            if (encompassed) {
                iter = boundaryVerts.erase(iter);
                vertInHole.insert(vIdx);
            }
            else iter++;
        }
    }

    result.face_inside = vector < int > {faceInHole.begin(), faceInHole.end()};
    result.vert_inside = vector < int > {vertInHole.begin(), vertInHole.end()};
    result.vert_on_boundary = vector < int > {boundaryVerts.begin(), boundaryVerts.end()};
    return result;
}

MeshRegion generateRandomMeshRegion(const CMesh& mesh, const std::vector<int>& seedVerts, int targetSize)
{
    set<int> inside_verts{ seedVerts.begin(), seedVerts.end() };
    set<int> inside_faces;
    set<int> boundary_verts;
    
    default_random_engine generator((unsigned int)time(NULL));
    for (int vIdx : inside_verts) {
        for (auto f : mesh.vert(vIdx)->getAdjacentFaces()) {
            inside_faces.insert(f->getFaceIndex());
            for (int k = 0; k < 3; ++k) {
                int faceVertIdx = f->vertIdx(k);
                if (!setHas(inside_verts, faceVertIdx))
                    boundary_verts.insert(faceVertIdx);
            }
        }
    }
    set<int> hole_verts = setCombine<int>(inside_verts, boundary_verts);

    while (inside_verts.size() < targetSize) 
    {
        uniform_int_distribution<int> distr1(0, (int)boundary_verts.size() - 1);
        std::set<int>::const_iterator it(boundary_verts.begin());
        std::advance(it, distr1(generator));
        int newHoleVert = *it;
        inside_verts.insert(newHoleVert);
        for (auto f : mesh.vert(newHoleVert)->getAdjacentFaces()) {
            if (inside_faces.find(f->getFaceIndex()) != inside_faces.end())
                continue;       // face already considered     
            inside_faces.insert(f->getFaceIndex());
            for (int k = 0; k < 3; ++k) {
                int faceVertIdx = f->vertIdx(k);
                if (inside_verts.find(faceVertIdx) == inside_verts.end()) {
                    boundary_verts.insert(faceVertIdx);
                    hole_verts.insert(faceVertIdx);
                }
            }
        }
        boundary_verts.erase(newHoleVert);

        /* consider faces whose all verts are hole_verts; include them in inside_faces */
        for (int vi : boundary_verts) {
            for (CFace* f : mesh.vert(vi)->getAdjacentFaces()) {
                if (setHas<int>(inside_faces, f->getFaceIndex())) continue;
                if (setHasAll<int>(hole_verts, f->getAllVertIdx())) {
                    inside_faces.insert(f->getFaceIndex());
                }
            }
        }

        /* consider boundary verts encompassed by inside_faces; include them in inside_verts */
        for (auto iter = boundary_verts.begin(); iter != boundary_verts.end();) {
            int vIdx = *iter;
            bool encompassed = true;
            for (auto f : mesh.vert(vIdx)->getAdjacentFaces()) {
                if (!setHas<int>(inside_faces, f->getFaceIndex())) {
                    // found a face not in hole, so vIdx is not in hole yet
                    encompassed = false;
                    break;
                }
            }
            if (encompassed) {
                iter = boundary_verts.erase(iter);
                inside_verts.insert(vIdx);
            }
            else iter++;
        }
    }

    MeshRegion result;
    result.face_inside = vector < int > {inside_faces.begin(), inside_faces.end()};
    result.vert_inside = vector < int > {inside_verts.begin(), inside_verts.end()};
    result.vert_on_boundary = vector < int > {boundary_verts.begin(), boundary_verts.end()};
    result.determineBoundaryHalfEdges(mesh);
    return result;
}

double compareCoordRMSE(const MeshCoordinates& coord1, const MeshCoordinates& coord2, const std::vector<int>& selctedVerts)
{
    double errorSum = 0;
    for (int vIdx : selctedVerts) {
        errorSum += (coord1.getVertCoord(vIdx) - coord2.getVertCoord(vIdx)).length2();
    }
    return std::sqrt(errorSum / selctedVerts.size());
}

MeshCoordinates addMeshNoise(CMesh& mesh, double phi, std::vector<int>& selectedVerts)
{
    assert(phi > 0 && phi < 1);
    const int vertCount = mesh.vertCount();
    const double avgLen = mesh.getAvgEdgeLength();
    double bbDiag = mesh.getBoundingBox().length() * 2;
    MeshCoordinates newCoord = mesh.getVertCoordinates();

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, phi);

    for (int vIdx : selectedVerts) {
        for (int a = 0; a < 3; ++a) {
            double noise = bbDiag * distribution(generator);
            newCoord(vIdx, a) += noise;
        }
    }

    return newCoord;
}

std::vector<int> getMeshRegionsFaces(const std::vector<MeshRegion>& vRegions)
{
    set<int> faces;
    for (auto& mr : vRegions)
        for (int vi : mr.face_inside)
            faces.insert(vi);

    return vector < int > {faces.begin(), faces.end()};
}

void refineMeshHoles3(CMesh &mesh, double lambda /*= 0.7*/)
{   
    if (!mesh.hasAttr(ZGeom::StrAttrMeshHoleRegions)) { std::cout << "Hole not found!\n"; return; }
    vector<MeshRegion>& vHoles = mesh.getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
    if (vHoles.empty()) { std::cout << "Hole not found!\n"; return; }

    vector<set<CVertex*>> vHoleInsideVerts;
    vector<set<CFace*>> vHoleInsideFaces;

    vector<double> est_hole_edge_lengths(vHoles.size());
    for (int holeIdx = 0; holeIdx < (int)vHoles.size(); ++holeIdx) {
        vHoles[holeIdx].determineBoundaryHalfEdges(mesh);
        est_hole_edge_lengths[holeIdx] = estimateHoleEdgeLength(mesh, vHoles[holeIdx], 3);
    }

    /* Perform Delaunay-like refinement */
    for (int holeIdx = 0; holeIdx < (int)vHoles.size(); ++holeIdx)
    {
        const MeshRegion& cur_hole = vHoles[holeIdx];
        const double preferred_edge_length = est_hole_edge_lengths[holeIdx];

        set<CVertex*> verts_in_hole;
        set<CFace*> faces_in_hole;
        for (int fi : cur_hole.face_inside) faces_in_hole.insert(mesh.getFace(fi));

        auto isHalfEdgeInside = [&faces_in_hole](CHalfEdge* he) {
            if (he->twinHalfEdge() == nullptr) return false;
            return faces_in_hole.find(he->twinHalfEdge()->getAttachedFace()) != faces_in_hole.end();
        };

        set<CHalfEdge*> boundary_he;
        for (int he_idx : cur_hole.he_on_boundary)
            boundary_he.insert(mesh.getHalfEdge(he_idx));

        while (true)
        {
            // 1. for each triangle inside hole, check if any face is splittable

            CVertex* new_vert = nullptr;
            for (CFace* face : faces_in_hole)
            {
                int obtuse_edge = ZGeom::triObtuseEdge(face->getAllEdgeLengths());
                if (obtuse_edge == -1)  // non-obtuse triangle: 3-split face
                {
                    Vec3d vc = face->calBarycenter();
                    bool splitTest = false;
#if 1
                    splitTest = true;
                    for (int m = 0; m < 3; ++m) {
                        if ((face->vert(m)->pos() - vc).length() <= lambda * preferred_edge_length) {
                            splitTest = false; break;
                        }
                    }                    
#else
                    double totalLength(0);
                    for (int m = 0; m < 3; ++m)
                        totalLength += (face->vert(m)->pos() - vc).length();
                    splitTest = (totalLength >= 3 * lambda * preferred_edge_length);
#endif
                    if (splitTest)
                    {
                        CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                        new_vert = mesh.faceSplit3(face);
                        
                        verts_in_hole.insert(new_vert);
                        for (CFace* newF : new_vert->getAdjacentFaces())
                            faces_in_hole.insert(newF);
                        
                        // relax edges opposite to the new added centroid
                        for (int k = 0; k < 3; ++k) {
                            if (isHalfEdgeInside(fe[k]))
                                mesh.relaxEdge(fe[k]);
                        }
                        break;
                    }
                }
                else // obtuse triangle: 2-split obtuse edge
                {
                    CHalfEdge *obtuse_he = face->getHalfEdge(obtuse_edge);
                    if (!isHalfEdgeInside(obtuse_he))
                        continue; // don't split boundary edges

                    bool splitTest = false;
#if 1
                    splitTest = (obtuse_he->length() >= 2 * lambda * preferred_edge_length);
#else
                    Vec3d vc = obtuse_he->midEdge();
                    double totalLength(0);
                    for (int m = 0; m < 3; ++m)
                        totalLength += (face->vert(m)->pos() - vc).length();
                    splitTest = (totalLength >= 3 * lambda * preferred_edge_length);
#endif                        
                    if (splitTest)
                    {
                        vector<CHalfEdge*> diamond_edges = { obtuse_he->nextHalfEdge(), obtuse_he->prevHalfEdge() };
                        if (obtuse_he->twinHalfEdge() != nullptr) {
                            CHalfEdge* obtuse_he2 = obtuse_he->twinHalfEdge();
                            diamond_edges.push_back(obtuse_he2->nextHalfEdge());
                            diamond_edges.push_back(obtuse_he2->prevHalfEdge());
                        }
                        new_vert = mesh.edgeSplit(obtuse_he);

                        verts_in_hole.insert(new_vert);
                        for (CFace* newF : new_vert->getAdjacentFaces())
                            faces_in_hole.insert(newF);

                        // relax edges opposite to the new added edge_center
                        for (int k = 0; k < (int)diamond_edges.size(); ++k) {
                            if (isHalfEdgeInside(diamond_edges[k]))
                                mesh.relaxEdge(diamond_edges[k]);
                        }
                        break;
                    }
                }
            }

            // 2. If no new points were inserted in step 1, hole filling process completes
            if (new_vert == nullptr) break;
            else if (verts_in_hole.size() % 10 == 0) {
                 std::cout << "#verts_in_hole: " << verts_in_hole.size() << '\n';
            }
            
            // 3. relax all edges in patch mesh
            vector<CHalfEdge*> inside_he;
            for (CVertex* pv : verts_in_hole) {
                for (CHalfEdge* he : pv->getHalfEdges())
                    inside_he.push_back(he);
            }

            int swapIterCount = 0;
            while (true) {
                bool noEdgeSwap = true;
                for (int i = 0; i < (int)inside_he.size(); ++i) {
                    if (mesh.relaxEdge(inside_he[i])) {
                        noEdgeSwap = false;
                    }
                }
                if (noEdgeSwap) break;
                if (swapIterCount++ >= 100) break;
            }
            if (swapIterCount >= 100) cout << "!!! Maximum swap reached!\n";

            // 4. go to step 1
        }

        vHoleInsideVerts.push_back(verts_in_hole);
        vHoleInsideFaces.push_back(faces_in_hole);
    } // for each hole

    vector<MeshRegion> newHoles(vHoles.size());
    for (int holeIdx = 0; holeIdx < (int)vHoles.size(); ++holeIdx)
    {
        vector<int> inside_verts;
        for (CVertex* pv : vHoleInsideVerts[holeIdx])
            inside_verts.push_back(pv->getIndex());
        vector<int> inside_faces;
        for (CFace* pf : vHoleInsideFaces[holeIdx])
            inside_faces.push_back(pf->getFaceIndex());

        newHoles[holeIdx].vert_inside = inside_verts;
        newHoles[holeIdx].face_inside = inside_faces;
        newHoles[holeIdx].vert_on_boundary = vHoles[holeIdx].vert_on_boundary;
        newHoles[holeIdx].determineBoundaryHalfEdges(mesh);

        std::cout << "[Hole #" << holeIdx
            << "] #vert: " << newHoles[holeIdx].vert_inside.size()
            << "; Est. edge length: " << est_hole_edge_lengths[holeIdx]
            << "; Actual edge length: " << calMeshRegionAvgEdgeLen(mesh, newHoles[holeIdx])
            << std::endl;
    }

    mesh.clearNonEssentialAttributes();
    gatherMeshStatistics(mesh);
    mesh.addAttr<vector<MeshRegion>>(newHoles, StrAttrMeshHoleRegions, AR_UNIFORM, AT_UNKNOWN);
}

bool refineMeshHoleByNum(CMesh &mesh, MeshRegion& hole, int new_vert_count)
{
    if (new_vert_count <= 0) return false;
    if (hole.face_inside.empty()) return false;
    if (!hole.vert_inside.empty()) return false;
    
    MeshRegion& cur_hole = hole;
    int num_new_vert = 0;

    set<CVertex*> verts_in_hole;
    set<CFace*> faces_in_hole;
    for (int fi : cur_hole.face_inside) 
        faces_in_hole.insert(mesh.getFace(fi));
    set<CHalfEdge*> boundary_he;
    for (CFace* f : faces_in_hole) {
        for (CHalfEdge* he : f->getAllHalfEdges()) {
            if (he->twinHalfEdge() == nullptr || !setHas(faces_in_hole, he->twinHalfEdge()->getAttachedFace()))
                boundary_he.insert(he);
        }
    }

    while (num_new_vert++ < new_vert_count)
    {
        // 1. for each triangle inside hole, check if any face is splittable
        concurrency::concurrent_vector<std::pair<CFace*, double>> face_priority;
        concurrency::parallel_for_each(faces_in_hole.begin(), faces_in_hole.end(), [&face_priority](CFace *f) {
            double length_sum(0);
            Vec3d vc = f->calBarycenter();
            for (int m = 0; m < 3; ++m)
                length_sum += (f->vert(m)->pos() - vc).length();
            face_priority.push_back(make_pair(f, length_sum));
        });

        CFace* face_to_split = nullptr;
        double max_len(0);
        for (int i = 0; i < (int)face_priority.size(); ++i) {
            if (face_priority[i].second > max_len) {
                face_to_split = face_priority[i].first;
                max_len = face_priority[i].second;
            }
        }
            
        assert(face_to_split);
        vector<CHalfEdge*> fe = face_to_split->getAllHalfEdges();
        CVertex* new_vert = mesh.faceSplit3(face_to_split);

        verts_in_hole.insert(new_vert);
        for (CFace* newF : new_vert->getAdjacentFaces())
            faces_in_hole.insert(newF);

        // relax edges opposite to the new added centroid
        for (int k = 0; k < 3; ++k) {
            if (!setHas(boundary_he, fe[k]))
                mesh.relaxEdge(fe[k]);
        }
    
        // 3. relax all edges in patch mesh
        vector<CHalfEdge*> inside_he;
        for (CVertex* pv : verts_in_hole) {
            for (CHalfEdge* he : pv->getHalfEdges())
                inside_he.push_back(he);
        }

        int swapIterCount = 0;
        while (true) {
            bool noEdgeSwap = true;
            for (int i = 0; i < (int)inside_he.size(); ++i) {
                if (mesh.relaxEdge(inside_he[i])) {
                    noEdgeSwap = false;
                }
            }
            if (noEdgeSwap) break;
            if (swapIterCount++ >= 100) break;
        }
        if (swapIterCount >= 100) cout << "!!Maximum swap reached!\n";    
    }

    cur_hole.face_inside.clear();
    cur_hole.vert_inside.clear();
    for (CFace* pf : faces_in_hole) cur_hole.face_inside.push_back(pf->getFaceIndex());
    for (CVertex* pv : verts_in_hole) cur_hole.vert_inside.push_back(pv->getIndex());

    std::cout << "[Hole #" << 0
        << "] #vert: " << cur_hole.vert_inside.size()
        << "; Est. edge length: " << cur_hole.adjacent_edge_length
        << "; Actual edge length: " << calMeshRegionAvgEdgeLen(mesh, cur_hole)
        << std::endl;


    hole.determineBoundaryHalfEdges(mesh);

    mesh.clearNonEssentialAttributes();
    gatherMeshStatistics(mesh);
    mesh.addAttr<vector<MeshRegion>>(vector < MeshRegion > {hole}, StrAttrMeshHoleRegions, AR_UNIFORM, AT_UNKNOWN);
    
    return true;
}

const std::vector<double>& getMeshCurvatures(CMesh& mesh, VertCurvature::CurvatureType curvature_type)
{
    string curv_type_name;
    switch (curvature_type)
    {
    case ZGeom::VertCurvature::MEAN:
        curv_type_name = StrAttrVertMeanCurvatures;
        break;
    case ZGeom::VertCurvature::GAUSS:
        curv_type_name = StrAttrVertGaussCurvatures;
        break;
    case ZGeom::VertCurvature::PRINCIPAL_1:
        curv_type_name = StrAttrVertPrincipalCurvatures1;
        break;
    case ZGeom::VertCurvature::PRINCIPAL_2:
        curv_type_name = StrAttrVertPrincipalCurvatures2;
        break;
    default:
        throw runtime_error("Unrecognized curvature type!");
        break;
    }

    if (!mesh.hasAttr(curv_type_name)) computeMeshCurvatures(mesh, true);
    
    return mesh.getAttrValue<vector<double>>(curv_type_name);
}

}   // end of namespace