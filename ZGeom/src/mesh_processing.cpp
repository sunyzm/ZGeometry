#include "mesh_processing.h"
#include <cstdlib>
#include <random>
#include <map>
#include <unordered_set>
#include <ppl.h>
#include <concurrent_vector.h>
#include "MatVecArithmetic.h"
#include "arithmetic.h"
#include "triBoxOverlap.h"

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
    mesh.addAttr<std::vector<ZGeom::Vec3d>>(vertNormals, CMesh::StrAttrVertNormal, AR_VERTEX, AT_VEC_VEC3);
}

std::vector<Vec3d> getMeshVertNormals(CMesh& mesh, VertNormalCalcMethod vnc /*= VN_AREA_WEIGHT*/)
{
    if (!mesh.hasAttr(CMesh::StrAttrVertNormal)) calMeshAttrVertNormals(mesh, vnc);

    return mesh.getAttrValue<vector<Vec3d>>(CMesh::StrAttrVertNormal);
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
    // 8. vertex voronoi area

    mesh.calAttrFaceNormals();
    calMeshAttrVertNormals(mesh, ZGeom::VN_AREA_WEIGHT);
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

    calMeshAttrMeanGaussCurvatures(mesh);
    calMeshAttrMixedVertAreas(mesh);
}

/********************************************************************************/
/* Adapted from AABB-triangle overlap test code                                 */
/* by Tomas Akenine-Möller                                                      */
/********************************************************************************/
bool testTriBoxOverlap(const std::vector<Vec3d>& triangle, Vec3d boxCenter, Vec3d boxHalfsize)
{
    float boxcenter[3] = { (float)boxCenter[0], (float)boxCenter[1], (float)boxCenter[2] };
    float boxhalfsize[3] = { (float)boxHalfsize[0], (float)boxHalfsize[1], (float)boxHalfsize[2] };
    float triverts[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            triverts[i][j] = (float)triangle[i][j];
    return ((1 == triBoxOverlap(boxcenter, boxhalfsize, triverts)) ? true : false);
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

std::vector<HoleBoundary> identifyMeshBoundaries(CMesh& mesh)
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
            for (CHalfEdge* he : mesh.m_vVertices[currentIndex]->m_HalfEdges) {
                if (he->isBoundaryEdge()) {
                    edgeIndex = he->getIndex(); break;
                }
            }
            currentIndex = mesh.m_vHalfEdges[edgeIndex]->getVertIndex(1);
            edgeLoop.push_back(edgeIndex);
        } while (currentIndex != i);

        boundaryEdges.push_back(edgeLoop);
    }

    // sort by number of edges of boundaries
    std::sort(boundaryEdges.begin(), boundaryEdges.end(),
        [](const vector<int>& v1, const vector<int>& v2) { return v1.size() < v2.size(); });

    vector<HoleBoundary> result(boundaryEdges.size());
    for (int i = 0; i < (int)boundaryEdges.size(); ++i) {
        result[i].he_on_boundary = boundaryEdges[i];
        for (int heIdx : boundaryEdges[i])
            result[i].vert_on_boundary.push_back(mesh.getHalfEdge(heIdx)->getVertIndex(0));
    }

    mesh.addAttr<vector<HoleBoundary>>(result, StrAttrMeshHoleBoundaries, AR_UNIFORM);
    return result;
}

const std::vector<HoleBoundary>& getMeshBoundaryLoops(CMesh &mesh)
{
    if (!mesh.hasAttr(StrAttrMeshHoleBoundaries)) identifyMeshBoundaries(mesh);
    return mesh.getAttrValue<vector<HoleBoundary>>(StrAttrMeshHoleBoundaries);
}

std::vector<std::vector<int>> getMeshBoundaryLoopVerts(CMesh &mesh)
{
    const std::vector<HoleBoundary>& vBoundaries = getMeshBoundaryLoops(mesh);
    vector<vector<int>> result;
    for (const HoleBoundary& hb : vBoundaries)
        result.push_back(hb.vert_on_boundary);
    return result;
}

std::vector<std::vector<int>> getMeshBoundaryLoopHalfEdges(CMesh &mesh)
{
    const std::vector<HoleBoundary>& vBoundaries = getMeshBoundaryLoops(mesh);
    vector<vector<int>> result;
    for (const HoleBoundary& hb : vBoundaries)
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
    const vector<vector<int>>& boundariesVert = getMeshBoundaryLoopVerts(mesh);
    vector<bool> result(mesh.vertCount(), false);
    for (int i = 0; i < boundariesVert.size(); ++i) {
        const vector<int>& vVertIdx = boundariesVert[i];
        if (vVertIdx.size() > MAX_HOLE_SIZE) continue;
        for (int vIdx : vVertIdx) result[vIdx] = true;
    }
    return result;
}

std::unique_ptr<CMesh> cutMesh(CMesh &oldMesh, const std::vector<int>& cutFaceIdx)
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

void triangulateMeshHoles(CMesh &mesh)
{ 
    vector<HoleBoundary> old_holes = getMeshBoundaryLoops(mesh);

    for (int holeIdx = 0; holeIdx < (int)old_holes.size(); ++holeIdx)
    {
        int nOldVerts = mesh.vertCount();
        int nOldEdges = mesh.halfEdgeCount();
        int nOldFaces = mesh.faceCount();
        const vector<int>& boundaryEdgeIdx = old_holes[holeIdx].he_on_boundary;

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

        cout << "boundary vert count: " << boundaryVertIdx.size() << endl;
        cout << "patching face count: " << patchTri.size() << endl;

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

        for (CFace *f : patchFaces) 
            old_holes[holeIdx].face_inside.push_back(f->getFaceIndex());
    }

    mesh.clearAttributes();
    gatherMeshStatistics(mesh);
    mesh.initNamedCoordinates();
    mesh.addAttr<vector<HoleBoundary>>(old_holes, StrAttrMeshHoleBoundaries, AR_UNIFORM, AT_UNKNOWN);    
}

/* "A multistep approach to restoration of locally undersampled meshes" (GMP 2008) */
void refineMeshHoles(CMesh &mesh, double lambda /*= 0.5*/)
{
    vector<HoleBoundary> oldHoles = getMeshBoundaryLoops(mesh);

    /* estimate edge length from hold surroundings */
    for (HoleBoundary& hb : oldHoles) ZGeom::estimateHoleEdgeLength(mesh, hb, 1);

    /* Perform Delaunay-like refinement */
    for (int holeIdx = 0; holeIdx < (int)oldHoles.size(); ++holeIdx)
    {
        int nOldVerts = mesh.vertCount();
        int nOldEdges = mesh.halfEdgeCount();
        int nOldFaces = mesh.faceCount();
        const vector<int>& boundaryEdgeIdx = oldHoles[holeIdx].he_on_boundary;
        const int N = (int)boundaryEdgeIdx.size();    // number of boundary vertices
        vector<int> boundaryVertIdx(N);
        vector<CVertex*> boundaryVertPtr(N);
        vector<Vec3d> boundaryVertPos(N);
        for (int i = 0; i < N; ++i) {
            boundaryVertIdx[i] = mesh.getHalfEdge(boundaryEdgeIdx[i])->getVertIndex(0);
            boundaryVertPtr[i] = mesh.vert(boundaryVertIdx[i]);
            boundaryVertPos[i] = boundaryVertPtr[i]->pos();
        }
        set<CHalfEdge*> boundaryHalfEdges;
        for (int ei : boundaryEdgeIdx) boundaryHalfEdges.insert(mesh.getHalfEdge(ei));
        set<int> vertsInHole;
        set<int> facesInHole = set<int>(oldHoles[holeIdx].face_inside.begin(), oldHoles[holeIdx].face_inside.end());
        const double preferred_edge_length = oldHoles[holeIdx].adjacent_edge_length;
        assert(preferred_edge_length > 0);

        while (true)
        {
            // 1. for each triangle inside hole, check if any face is splittable
            int selectedFIdx = -1;
            for (int fIdx : facesInHole) {
                CFace* face = mesh.getFace(fIdx);
                Vec3d vc = face->calBarycenter();
                CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                bool splitTest = true;
                for (int m = 0; m < 3; ++m) {
                    if ((fv[m]->pos() - vc).length() <= lambda * preferred_edge_length) {
                        splitTest = false; break;
                    }
                }
                if (splitTest) { selectedFIdx = fIdx; break; }
            }

            if (selectedFIdx != -1) // split the selected face
            {
                CFace *face = mesh.getFace(selectedFIdx);
                CHalfEdge *he[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                CVertex* centroid = mesh.faceSplit3(face->getFaceIndex());

                vertsInHole.insert(centroid->getIndex());
                facesInHole.erase(selectedFIdx);
                for (CFace* newF : centroid->getAdjacentFaces())
                    facesInHole.insert(newF->getFaceIndex());

                // relax edges opposite to the new added centroid
                for (int k = 0; k < 3; ++k) {
                    if (boundaryHalfEdges.find(he[k]) == boundaryHalfEdges.end()) { // not boundary half-edge
                        mesh.relaxEdge(he[k]);
                    }
                }
            }
            else break;  // 2. no new points were inserted in step 1, hole filling process completes

            // 3. relax all edges in patch mesh
            int swapIterCount = 0;
            while (true && swapIterCount++ < 100) {
                bool noEdgeSwap = true;
                for (int i = nOldEdges; i < mesh.m_vHalfEdges.size(); ++i) {
                    if (mesh.relaxEdge(mesh.getHalfEdge(i))) {
                        noEdgeSwap = false;
                        // break;
                    }
                }
                if (noEdgeSwap) break;
            }            

            // 4. go to step 1
        }

        oldHoles[holeIdx].face_inside = vector<int>(facesInHole.begin(), facesInHole.end());
        oldHoles[holeIdx].vert_inside = vector<int>(vertsInHole.begin(), vertsInHole.end());
    }

    mesh.clearAttributes();
    gatherMeshStatistics(mesh);
    mesh.initNamedCoordinates();
    mesh.addAttr<vector<HoleBoundary>>(oldHoles, StrAttrMeshHoleBoundaries, AR_UNIFORM, AT_UNKNOWN);
}

/* "Filling holes in meshes" (SGP 2003) */
void refineMeshHoles2(CMesh &mesh, double lambda /*= 2.0*/)
{
    vector<HoleBoundary> oldHoles = getMeshBoundaryLoops(mesh);

    /* Perform Delaunay-like refinement */
    for (int holeIdx = 0; holeIdx < (int)oldHoles.size(); ++holeIdx)
    {
        int nOldVerts = mesh.vertCount();
        int nOldEdges = mesh.halfEdgeCount();
        int nOldFaces = mesh.faceCount();
        const vector<int>& boundaryEdgeIdx = oldHoles[holeIdx].he_on_boundary;
        const int N = (int)boundaryEdgeIdx.size();    // number of boundary vertices
        vector<int> boundaryVertIdx(N);
        vector<CVertex*> boundaryVertPtr(N);
        vector<Vec3d> boundaryVertPos(N);
        for (int i = 0; i < N; ++i) {
            boundaryVertIdx[i] = mesh.getHalfEdge(boundaryEdgeIdx[i])->getVertIndex(0);
            boundaryVertPtr[i] = mesh.vert(boundaryVertIdx[i]);
            boundaryVertPos[i] = boundaryVertPtr[i]->pos();
        }
        set<CHalfEdge*> boundaryHalfEdges;
        for (int ei : boundaryEdgeIdx) boundaryHalfEdges.insert(mesh.getHalfEdge(ei));
        set<int> vertsInHole;
        set<int> facesInHole = set<int>(oldHoles[holeIdx].face_inside.begin(), oldHoles[holeIdx].face_inside.end());
        map<int, double> lengthAttr;
        for (int i = 0; i < N; ++i) {
            CVertex* pv = boundaryVertPtr[i];
            double length_sum(0);
            double valid_valence = 0;
            for (CHalfEdge *he : pv->m_HalfEdges) {
                if (facesInHole.find(he->getAttachedFace()->getFaceIndex()) != facesInHole.end()) {
                    length_sum += he->length();
                    valid_valence++;
                }
            }            
            lengthAttr[pv->getIndex()] = length_sum / valid_valence;
        }

        while (true)
        {
            // 1. for each triangle inside hole, check if any face is splittable
            int selectedFIdx = -1;
            for (int fIdx : facesInHole) {
                CFace* face = mesh.getFace(fIdx);
                Vec3d vc = face->calBarycenter();
                CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                double vcLenAttr = (lengthAttr[fv[0]->getIndex()] + lengthAttr[fv[1]->getIndex()] + lengthAttr[fv[2]->getIndex()]) / 3.0;
                bool splitTest = true;
                for (int m = 0; m < 3; ++m) {
                    if ((fv[m]->pos() - vc).length() <= lambda * vcLenAttr) {
                        splitTest = false; break;
                    }
                }
                if (splitTest) { selectedFIdx = fIdx; break; }
            }

            if (selectedFIdx != -1) // split the selected face
            {
                CFace *face = mesh.getFace(selectedFIdx);
                CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                double vcLenAttr = (lengthAttr[fv[0]->getIndex()] + lengthAttr[fv[1]->getIndex()] + lengthAttr[fv[2]->getIndex()]) / 3.0;
                CVertex* centroid = mesh.faceSplit3(face->getFaceIndex());
                lengthAttr[centroid->getIndex()] = vcLenAttr;

                vertsInHole.insert(centroid->getIndex());
                facesInHole.erase(selectedFIdx);
                for (CFace* newF : centroid->getAdjacentFaces())
                    facesInHole.insert(newF->getFaceIndex());

                // relax edges opposite to the new added centroid
                for (int k = 0; k < 3; ++k) {
                    if (boundaryHalfEdges.find(fe[k]) == boundaryHalfEdges.end()) {
                        mesh.relaxEdge(fe[k]);
                    }
                }                
            }             
            else break;  // 2. no new points were inserted in step 1, hole filling process completes

            // 3. relax all edges in patch mesh
            int swapIterCount = 0;
            while (true && swapIterCount++ < 100) {
                bool noEdgeSwap = true;
                for (int i = nOldEdges; i < mesh.m_vHalfEdges.size(); ++i) {
                    if (mesh.relaxEdge(mesh.getHalfEdge(i))) {
                        noEdgeSwap = false; 
                        // break;
                    }
                }
                if (noEdgeSwap) break;
            }
            std::cout << "Swap iterations count: " << swapIterCount << std::endl;

            // 4. go to step 1
        }

        oldHoles[holeIdx].face_inside = vector<int>(facesInHole.begin(), facesInHole.end());
        oldHoles[holeIdx].vert_inside = vector<int>(vertsInHole.begin(), vertsInHole.end());
    }

    mesh.clearAttributes();
    gatherMeshStatistics(mesh);
    mesh.initNamedCoordinates();
    mesh.addAttr<vector<HoleBoundary>>(oldHoles, StrAttrMeshHoleBoundaries, AR_UNIFORM, AT_UNKNOWN);
}

void refineMeshHoles3(CMesh &mesh, double lambda /*= 2.0*/)
{
    vector<HoleBoundary> oldHoles = getMeshBoundaryLoops(mesh);

    /* Perform Delaunay-like refinement */
    for (int holeIdx = 0; holeIdx < (int)oldHoles.size(); ++holeIdx)
    {
        int nOldVerts = mesh.vertCount();
        int nOldEdges = mesh.halfEdgeCount();
        int nOldFaces = mesh.faceCount();
        const vector<int>& boundaryEdgeIdx = oldHoles[holeIdx].he_on_boundary;
        const int N = (int)boundaryEdgeIdx.size();    // number of boundary vertices
        vector<int> boundaryVertIdx(N);
        vector<CVertex*> boundaryVertPtr(N);
        vector<Vec3d> boundaryVertPos(N);
        for (int i = 0; i < N; ++i) {
            boundaryVertIdx[i] = mesh.getHalfEdge(boundaryEdgeIdx[i])->getVertIndex(0);
            boundaryVertPtr[i] = mesh.vert(boundaryVertIdx[i]);
            boundaryVertPos[i] = boundaryVertPtr[i]->pos();
        }
        set<CHalfEdge*> boundaryHalfEdges;
        for (int ei : boundaryEdgeIdx) boundaryHalfEdges.insert(mesh.getHalfEdge(ei));

        set<int> facesInHole = set<int>(oldHoles[holeIdx].face_inside.begin(), oldHoles[holeIdx].face_inside.end());
        map<int, double> lengthAttr;
        for (int i = 0; i < N; ++i) {
            CVertex* pv = boundaryVertPtr[i];
            double avgLen(0);
            double valid_valence = 0;
            for (CHalfEdge *he : pv->m_HalfEdges) {
                if (facesInHole.find(he->getAttachedFace()->getFaceIndex()) != facesInHole.end()) {
                    avgLen += he->length();
                    valid_valence++;
                }
            }
            avgLen /= valid_valence;
            lengthAttr[pv->getIndex()] = avgLen;
        }

        bool noFaceSplit(true), noEdgeSwap(true);
        while (true)
        {
            vector<CFace*> splitCandidates;
            for (int faceIdx : facesInHole)
            {
                CFace* face = mesh.getFace(faceIdx);
                Vec3d vc = mesh.getFace(faceIdx)->calBarycenter();
                CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                double vcLenAttr = (lengthAttr[fv[0]->getIndex()] + lengthAttr[fv[1]->getIndex()] + lengthAttr[fv[2]->getIndex()]) / 3.0;
                bool splitTest = true;
                for (int m = 0; m < 3; ++m) {
                    double a = lambda * (fv[m]->pos() - vc).length();
                    if (a <= max(lengthAttr[fv[m]->getIndex()], vcLenAttr)) {
                        splitTest = false; break;
                    }
                }
                if (splitTest) splitCandidates.push_back(face);
            }

            if (splitCandidates.empty()) break;
            else {
                for (CFace* face : splitCandidates) {
                    facesInHole.erase(face->getFaceIndex());
                    CHalfEdge *fe[3] = { face->getHalfEdge(0), face->getHalfEdge(1), face->getHalfEdge(2) };
                    CVertex *fv[3] = { face->vert(0), face->vert(1), face->vert(2) };
                    double vcLenAttr = (lengthAttr[fv[0]->getIndex()] + lengthAttr[fv[1]->getIndex()] + lengthAttr[fv[2]->getIndex()]) / 3.0;
                    CVertex* centroid = mesh.faceSplit3(face->getFaceIndex());
                    lengthAttr[centroid->getIndex()] = vcLenAttr;
                    for (CFace* newF : centroid->getAdjacentFaces())
                        facesInHole.insert(newF->getFaceIndex());
                }
            }

            // relaxing all interior half-edges
            int swapCount = 0;
            while (true && swapCount++ < 100) {
                noEdgeSwap = true;
                for (int i = nOldEdges; i < mesh.m_vHalfEdges.size(); ++i) {
                    if (boundaryHalfEdges.find(mesh.getHalfEdge(i)->twinHalfEdge()) != boundaryHalfEdges.end()) continue;
                    if (mesh.relaxEdge(mesh.getHalfEdge(i))) {
                        noEdgeSwap = false; // break;
                    }
                }
                if (noEdgeSwap) break;
            }
        }

        oldHoles[holeIdx].face_inside = vector<int>(facesInHole.begin(), facesInHole.end());
        oldHoles[holeIdx].vert_inside.clear();
        for (int vIdx = nOldVerts; vIdx < (int)mesh.vertCount(); ++vIdx)
            oldHoles[holeIdx].vert_inside.push_back(vIdx);
    }

    mesh.clearAttributes();
    gatherMeshStatistics(mesh);
    mesh.initNamedCoordinates();
    mesh.addAttr<vector<HoleBoundary>>(oldHoles, StrAttrMeshHoleBoundaries, AR_UNIFORM, AT_UNKNOWN);
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

void estimateHoleEdgeLength(CMesh& mesh, HoleBoundary& hole, int ring /* = 1 */)
{
    set<int> vNeighborVerts = meshMultiVertsAdjacentVerts(mesh, hole.vert_on_boundary, ring, true);
    set<int> vNeighborHalfEdges;
    for (int vi : vNeighborVerts) {
        for (CHalfEdge* he : mesh.vert(vi)->getHalfEdges())
            vNeighborHalfEdges.insert(he->getIndex());
    }
    double lengthSum(0);
    for (int eIdx : vNeighborHalfEdges)
        lengthSum += mesh.getHalfEdge(eIdx)->length();
 
    hole.adjacent_edge_length = lengthSum / (double)vNeighborHalfEdges.size();
}

}   // end of namespace