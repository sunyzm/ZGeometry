#ifndef ZGEOM_MESH_PROCESSING_H
#define ZGEOM_MESH_PROCESSING_H
#include <functional>
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Plane.h"
#include "EigenSystem.h"
#include "Mesh.h"

namespace ZGeom {

const std::string StrAttrVertMixedArea = "vert_scalar_mixed_area";
const std::string StrAttrVertGaussCurvatures = "vert_gauss_curvature";
const std::string StrAttrVertMeanCurvatures = "vert_mean_curvature";

struct ResultDistPointTriangle { Vec3d closestPoint; double distance; };
ResultDistPointTriangle distPointTriangle(Vec3d point, const std::vector<ZGeom::Vec3d>& triangle);

struct ResultDistPointPlane { Vec3d closestPoint; double distance, signedDistance; };
ResultDistPointPlane distPointPlane(Vec3d point, const Plane3& plane);

double distPoint2Mesh(Vec3d p, CMesh& mesh);
double computeMeanHausdorffDistance(const CMesh& mesh1, const CMesh& mesh2);


/* spectral geometry*/
DenseMatrixd calSpectralKernelMatrix(const EigenSystem& hb, double t, std::function<double(double, double)> gen);
std::vector<double> calSpectralKernelSignature(const EigenSystem& es, double t, std::function<double(double, double)> gen);
double calHK(const EigenSystem& es, int i, int j, double t);
DenseMatrixd calHeatKernelMatrix(const EigenSystem& es, double t);
std::vector<double> calHeatKernelSignature(const EigenSystem& es, double t);
double calHeatTrace(const EigenSystem& es, double t);
double calBiharmonicDist(const EigenSystem& es, int v1, int v2);


/* differential geometry */
enum MeshMatrixType
{
    MM_ADJACENCY, MM_DEGREE, MM_INV_DEGREE, MM_WALK,
    MM_GRAPH_LAPLACE, MM_NORMALIZED_GRAPH_LAPLACE,
    MM_INV_LENGTH
};
SparseMatrix<double> constructMeshMatrix(const CMesh& mesh, MeshMatrixType mmt);

void getMeshGraphCSR(const CMesh& mesh, std::vector<int>& xadj, std::vector<int>& adjncy);

std::vector<Vec3d> computeMeshFaceNormals(const CMesh& mesh);
enum VertNormalCalcMethod {VN_UNIFORM_WEIGHT, VN_AREA_WEIGHT, VN_CENTER_DIST_WEIGHT};
std::vector<Vec3d> computeMeshVertNormals(const CMesh& mesh, VertNormalCalcMethod vnc = VN_AREA_WEIGHT);
void calMeshAttrVertNormals(CMesh& mesh, VertNormalCalcMethod vnc = VN_AREA_WEIGHT);

enum MeshVertAreaScheme {VA_UNIFORM, VA_MIXED_VORONOI};
std::vector<double> computeMeshVertArea(const CMesh& mesh, MeshVertAreaScheme scheme = VA_MIXED_VORONOI);
void calMeshAttrMixedVertAreas(CMesh& mesh);
const std::vector<double>& getMeshVertMixedAreas(CMesh& mesh);

struct ResultMeshMeanGaussCurvatures 
{ 
    std::vector<double> mean_curvatures, gauss_curvatures; 
    
    ResultMeshMeanGaussCurvatures() = default;
    ResultMeshMeanGaussCurvatures(const ResultMeshMeanGaussCurvatures&) = default;
    ResultMeshMeanGaussCurvatures(ResultMeshMeanGaussCurvatures&& c2) { 
        mean_curvatures = std::move(c2.mean_curvatures); 
        gauss_curvatures = std::move(c2.gauss_curvatures); 
    }
};
ResultMeshMeanGaussCurvatures computeMeshMeanGaussCurvatures(CMesh &mesh);
void calMeshAttrMeanGaussCurvatures(CMesh &mesh);
const std::vector<double>& getMeshMeanCurvatures(CMesh &mesh);
const std::vector<double>& getMeshGaussCurvatures(CMesh &mesh);

struct VertPrincipalCurvatures
{
public:
    double curvatures[2];
    Vec3d directions[2];

    VertPrincipalCurvatures(const VertPrincipalCurvatures& vc2) { *this = vc2; }

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

void gatherMeshStatistics(CMesh& mesh);

}   // end of namespace

#endif