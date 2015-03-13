#ifndef ZGEOM_MESH_PROCESSING_H
#define ZGEOM_MESH_PROCESSING_H
#include <functional>
#include <array>
#include <memory>
#include <algorithm>
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Plane.h"
#include "EigenSystem.h"
#include "Mesh.h"

namespace ZGeom {

const int MAX_HOLE_SIZE = 200;

const std::string StrAttrVertMixedArea = "vert_scalar_mixed_area";
const std::string StrAttrVertGaussCurvatures = "vert_gauss_curvature";
const std::string StrAttrVertMeanCurvatures = "vert_mean_curvature";
const std::string StrAttrVertPrincipalCurvatures1 = "vert_principal_curvature_1";
const std::string StrAttrVertPrincipalCurvatures2 = "vert_principal_curvature_2";
const std::string StrAttrMeshHoleBoundaries = "mesh_hole_boundaries";
const std::string StrAttrMeshGeneratedHols = "mesh_generated_holes";


struct ResultDistPointTriangle { Vec3d closestPoint; double distance; };
ResultDistPointTriangle distPointTriangle(Vec3d point, const std::vector<Vec3d>& triangle);
ResultDistPointTriangle distPointTriangle2(Vec3d point, const std::vector<Vec3d>& triangle);

struct ResultDistPointPlane { Vec3d closestPoint; double distance, signedDistance; };
ResultDistPointPlane distPointPlane(Vec3d point, const Plane3& plane);

bool testTriBoxOverlap(const std::vector<Vec3d>& triangle, Vec3d boxCenter, Vec3d boxHalfsize);

int triObtuseEdge(const std::vector<Vec3d>& triVerts);

std::set<int> meshMultiVertsAdjacentVerts(const CMesh& mesh, const std::vector<int>& vert, int ring, bool inclusive = true);

double distPointMesh(Vec3d p, CMesh& mesh);
std::vector<double> computeVertMeshDist(CMesh &mesh1, CMesh &mesh2);
double computeMeshMeanError(CMesh &mesh1, CMesh &mesh2);
double computeSymMeshMeanError(CMesh &mesh1, CMesh &mesh2);
double computeMeshRMSE(CMesh &mesh1, CMesh &mesh2);
double computeSymMeshRMSE(CMesh &mesh1, CMesh &mesh2);
double computeHausdorffDistance(CMesh &mesh1, CMesh &mesh2);
double computeSymHausdorffDistance(CMesh &mesh1, CMesh &mesh2);

enum MeshDistMeasure {MEAN_ERROR, SYM_MEAN_ERROR, RMSE, SYM_RMSE, HAUSDORFF, SYM_HAUSDORFF};
double distSubMesh(CMesh &mesh1, const std::vector<int>& faces1, CMesh &mesh2, const std::vector<int>& faces2, MeshDistMeasure measure = RMSE);

struct HoleBoundary
{
    std::vector<int> vert_on_boundary;
    std::vector<int> he_on_boundary;
    std::vector<int> vert_inside;
    std::vector<int> face_inside;  

    bool is_outer_boundary;
    double adjacent_edge_length;

    // constructors
    HoleBoundary() : is_outer_boundary(false), adjacent_edge_length(-1) { }
    HoleBoundary(const HoleBoundary&) = default;
    HoleBoundary& operator = (HoleBoundary&& hb) {
        vert_on_boundary = std::move(hb.vert_on_boundary);
        he_on_boundary = std::move(hb.he_on_boundary);
        vert_inside = std::move(hb.vert_inside);
        face_inside = std::move(hb.face_inside);
        is_outer_boundary = hb.is_outer_boundary;
        adjacent_edge_length = hb.is_outer_boundary;
        return *this;
    }
    HoleBoundary(HoleBoundary&& hb) { *this = std::move(hb); }
    const std::vector<int>& getInsideFaceIdx() const { return face_inside; }
};
ZGeom::HoleBoundary autoGenerateHole(const CMesh& mesh, int seedVert, int holeSize);
ZGeom::HoleBoundary autoGenerateHole(const CMesh& mesh, const std::vector<int>& seedVerts, int totalSize);



std::vector<HoleBoundary> identifyMeshBoundaries(CMesh& mesh); // compute number of (connective) boundaries
void estimateHoleEdgeLength(CMesh& mesh, HoleBoundary& hole, int ring = 1);
double calAvgHoleEdgeLength(CMesh& mesh, HoleBoundary& hole);

const std::vector<HoleBoundary>& getMeshBoundaryLoops(CMesh &mesh);
std::vector<std::vector<int>> getMeshBoundaryLoopVerts(CMesh &mesh);
std::vector<std::vector<int>> getMeshBoundaryLoopHalfEdges(CMesh &mesh);
int calMeshGenus(CMesh &mesh);
std::vector<bool> getMeshVertsOnHoles(CMesh &mesh);

std::unique_ptr<CMesh> cutFromMesh(CMesh &oldMesh, const std::vector<int>& cutFaceIdx);
std::unique_ptr<CMesh> cutMeshTo(CMesh &oldMesh, const std::vector<int>& cutFaceIdx);

struct WeightSet
{
    double angle, area;

    WeightSet() : angle(0), area(0) {}
    WeightSet(double a, double s) : angle(a), area(s) {}
    bool operator < (const WeightSet &w2) const {
        return angle < w2.angle || (angle == w2.angle && area < w2.area);
    }
    WeightSet operator+ (const WeightSet &w2) const {
        WeightSet result;
        result.angle = std::max(this->angle, w2.angle);
        result.area = this->area + w2.area;
        return result;
    }
};
void triangulateMeshHoles(CMesh &oldMesh);
void refineMeshHoles(CMesh &oldMesh, double lambda = 2.0);  // lambda: density control factor
void refineMeshHoles2(CMesh &oldMesh, double lambda = 2.0);
void refineMeshHoles3(CMesh &oldMesh, double lamdba = 2.0);

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
std::vector<Vec3d> getMeshVertNormals(CMesh& mesh, VertNormalCalcMethod vnc = VN_AREA_WEIGHT);

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