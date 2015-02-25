#ifndef ZGEOM_MESH_PROCESSING_H
#define ZGEOM_MESH_PROCESSING_H
#include <functional>
#include "Mesh.h"
#include "SparseMatrix.h"
#include "EigenSystem.h"
#include "DenseMatrix.h"

namespace ZGeom {

double distPointTriangle(Vec3d p, const std::vector<ZGeom::Vec3d>& T);
double distPoint2Mesh(Vec3d p, CMesh& mesh);
double computeMeanHausdorffDistance(const CMesh& mesh1, const CMesh& mesh2);

enum MeshMatrixType 
{ 
	MM_ADJACENCY, MM_DEGREE, MM_INV_DEGREE, MM_WALK, 
	MM_GRAPH_LAPLACE, MM_NORMALIZED_GRAPH_LAPLACE,
    MM_INV_LENGTH
};
SparseMatrix<double> constructMeshMatrix(const CMesh& mesh, MeshMatrixType mmt);

void getMeshGraphCSR(const CMesh& mesh, std::vector<int>& xadj, std::vector<int>& adjncy);

DenseMatrixd calSpectralKernelMatrix(const EigenSystem& hb, double t, std::function<double(double, double)> gen);
std::vector<double> calSpectralKernelSignature(const EigenSystem& es, double t, std::function<double(double, double)> gen);

double calHK(const EigenSystem& es, int i, int j, double t);
DenseMatrixd calHeatKernelMatrix(const EigenSystem& es, double t);
std::vector<double> calHeatKernelSignature(const EigenSystem& es, double t);
double calHeatTrace(const EigenSystem& es, double t);

double calBiharmonicDist(const EigenSystem& es, int v1, int v2);

std::vector<Vec3d> computeMeshFaceNormals(const CMesh& mesh);

enum VertNormalCalcMethod {VN_UNIFORM_WEIGHT, VN_AREA_WEIGHT};
std::vector<Vec3d> computeMeshVertNormals(const CMesh& mesh, VertNormalCalcMethod vnc = VN_AREA_WEIGHT);

struct VertPrincipalCurvatures
{
public:
    double curvatures[2];
    Vec3d directions[2];

    VertPrincipalCurvatures(const VertPrincipalCurvatures& vc2)
    {
        for (int k = 0; k < 2; ++k) {
            this->curvatures[k] = vc2.curvatures[k];
            this->directions[k] = vc2.directions[k];
        }
    }

    const VertPrincipalCurvatures& operator=(const VertPrincipalCurvatures& vc2)
    {
        for (int k = 0; k < 2; ++k) {
            this->curvatures[k] = vc2.curvatures[k];
            this->directions[k] = vc2.directions[k];
        }
        return *this;
    }
};

std::vector<VertPrincipalCurvatures> estimateMeshPrincipalCurvatures(const CMesh& mesh);

std::vector<int> randomHoleVertex(const CMesh& mesh, int hole_size, int seed = -1);
std::vector<int> randomHoleVertex(const CMesh& mesh, int total_size, const std::vector<int>& seeds);

}   // end of namespace

#endif