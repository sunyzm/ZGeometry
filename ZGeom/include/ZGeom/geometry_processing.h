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
#include "curvature.h"
#include "spectral_analysis.h"
#include "mesh_primitives.h"

namespace ZGeom {

const std::string StrAttrVertMixedArea = "vert_scalar_mixed_area";
const std::string StrAttrVertGaussCurvatures = "vert_gauss_curvature";
const std::string StrAttrVertMeanCurvatures = "vert_mean_curvature";
const std::string StrAttrVertPrincipalCurvatures1 = "vert_principal_curvature_1";
const std::string StrAttrVertPrincipalCurvatures2 = "vert_principal_curvature_2";
const std::string StrAttrVertAllCurvatures = "vert_all_curvatures";
const std::string StrAttrVertShapeIndex = "vert_shape_index";
const std::string StrAttrMeshHoleRegions = "mesh_hole_regions";
const std::string StrAttrManualHoleRegions = "mesh_manually_generated_hole_regions";

struct ResultDistPointTriangle { Vec3d closestPoint; double distance; };
ResultDistPointTriangle distPointTriangle(Vec3d point, const std::vector<Vec3d>& triangle);
ResultDistPointTriangle distPointTriangle2(Vec3d point, const std::vector<Vec3d>& triangle);

struct ResultDistPointPlane { Vec3d closestPoint; double distance, signedDistance; };
ResultDistPointPlane distPointPlane(Vec3d point, const Plane3& plane);

bool testTriBoxOverlap(const std::vector<Vec3d>& triangle, Vec3d boxCenter, Vec3d boxHalfsize);

int triObtuseEdge(const std::vector<Vec3d>& triVerts);

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

MeshCoordinates addMeshNoise(CMesh& mesh, double phi, std::vector<int>& selectedVerts);

std::vector<int> getMeshRegionsInsideVerts(const std::vector<MeshRegion>& vRegions);
std::vector<int> getMeshRegionsBoundaryVerts(const std::vector<MeshRegion>& vRegions);
std::vector<int> getMeshRegionsFaces(const std::vector<MeshRegion>& vRegions);
void mergeMeshRegions(CMesh& mesh, std::vector<MeshRegion>& vRegions);

std::vector<MeshRegion*> getMeshHoleRegions(CMesh& mesh);
ZGeom::MeshRegion generateRandomMeshRegion(const CMesh& mesh, int seedVert, int holeSize);
ZGeom::MeshRegion generateRandomMeshRegion(const CMesh& mesh, const std::vector<int>& seedVerts, int totalSize);
ZGeom::MeshRegion generateRingMeshRegion(const CMesh& mesh, int seedVert, int ring);

std::vector<MeshRegion> identifyMeshBoundaries(CMesh& mesh); // compute number of (connective) boundaries
double estimateHoleEdgeLength(CMesh& mesh, MeshRegion& hole, int ring = 1);
double calMeshRegionAvgEdgeLen(CMesh& mesh, MeshRegion& hole);

const std::vector<MeshRegion>& getMeshBoundaryLoops(CMesh &mesh);
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

void refineMeshHoles(CMesh &oldMesh, double lambda = 0.7);  // lambda: density control factor
void refineMeshHoles2(CMesh &oldMesh, double lambda = 0.7);
void refineMeshHoles3(CMesh &oldMesh, double lambda = 0.7);
bool refineMeshHoleByNum(CMesh &oldMesh, MeshRegion& hole, int nNewVerts);


/* differential geometry */
//
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
const std::vector<Vec3d>& getMeshVertNormals(const CMesh& mesh, VertNormalCalcMethod vnc = VN_AREA_WEIGHT);
const std::vector<Vec3d>& getMeshVertNormals(CMesh& mesh, VertNormalCalcMethod vnc = VN_AREA_WEIGHT);

enum MeshVertAreaScheme {VA_UNIFORM, VA_MIXED_VORONOI};
std::vector<double> computeMeshVertArea(const CMesh& mesh, MeshVertAreaScheme scheme = VA_MIXED_VORONOI);
void calMeshAttrMixedVertAreas(CMesh& mesh);
const std::vector<double>& getMeshVertMixedAreas(CMesh& mesh);

//* curvatures *//
//
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

ZGeom::VertCurvature calVertCurvature(const CMesh& mesh, int v_idx, bool compute_principal = true);
void computeMeshCurvatures(CMesh& mesh, bool compute_principal = true);
const std::vector<double>& getMeshCurvatures(CMesh& mesh, VertCurvature::CurvatureType curvature_type);
void computeShapeIndex(CMesh& mesh);

std::vector<int> randomHoleVertex(const CMesh& mesh, int hole_size, int seed = -1);
std::vector<int> randomHoleVertex(const CMesh& mesh, int total_size, const std::vector<int>& seeds);

void gatherMeshStatistics(CMesh& mesh);

double compareCoordRMSE(const MeshCoordinates& coord1, const MeshCoordinates& coord2, const std::vector<int>& selctedVerts);

}   // end of namespace

#endif
