#include "MeshProcessing.h"
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
    for (int k = 0; k < es.eigVecCount(); ++k)	{
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
        double weightSum(0);
        for (int fj : neighborFaceIdx) {
            normalSum += faceAreas[fj] * faceNormals[fj];
            weightSum += faceAreas[fj];
        }
        vertNormals[vi] = normalSum / weightSum;
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

}   // end of namespace