#pragma once
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <engine.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/Mesh.h>
#include "MeshLaplacian.h"
#include "global.h"

enum KernelType { HEAT_KERNEL, MHW_KERNEL, SGW_KERNEL, BIHARMONIC_KERNEL };
enum DistanceType { DISTANCE_GEODESIC, DISTANCE_BIHARMONIC, DISTANCE_HK, DISTANCE_MHW };
enum PointSimilarityType {SIM_TYPE_1, SIM_TYPE_2, SIM_TYPE_3, SIM_TYPE_COUNT};

enum FeatureID { FEATURE_ID	= 0x0200, FEATURE_NEIGHBORS, FEATURE_HKS,
				 FEATURE_MHWS, FEATURE_MULTI_HKS, FEATURE_DEMO, FEATURE_DEMO2,
				 FEATURE_SGW_SOMP,
				 FEATURE_ID_COUNT};

double transferFunc2(double lambda, double t);	// Mexican-hat wavelet
double heatKernelTransferFunc(double lambda, double t);  // heat kernel
double mhwTransferFunc1(double lambda, double t);	// Mexican-hat wavelet (Tingbo's)

typedef double (*TransferFunc)(double, double);
typedef double (*ScalelessTransferFunc)(double);

class DifferentialMeshProcessor
{
public:
	DifferentialMeshProcessor();
	DifferentialMeshProcessor(CMesh* tm, CMesh* originalMesh);
	~DifferentialMeshProcessor();
	void init(CMesh* tm);
	void init_lite(CMesh* tm, CMesh* originalMesh);	
	const CMesh* getMesh_const() const { return mMesh; }
	const CMesh* getOriMesh_const() const { return mOriMesh; }
	CMesh* getMesh() const { return mMesh; }

	/* Laplacian and MHB related */
	void	constructLaplacian(LaplacianType laplacianType = CotFormula);
	void	decomposeLaplacian(int nEigFunc, LaplacianType laplacianType = CotFormula);
	const MeshLaplacian& getMeshLaplacian(LaplacianType laplacianType) const { return mMeshLaplacians[laplacianType]; }
	bool hasLaplacian(LaplacianType laplacianType) { return mMeshLaplacians[laplacianType].isLaplacianConstructed(); }
	bool isLaplacianDecomposed(LaplacianType laplacianType) { return !mMHBs[laplacianType].empty(); }
	void loadMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	void saveMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	const ZGeom::EigenSystem& getMHB(LaplacianType laplacianType) const { return mMHBs[laplacianType]; }
	bool isMHBCacheValid(const std::string& pathMHB, int eigenCount);
	std::string generateMHBPath(const std::string& prefix, LaplacianType laplacianType);
	const ZGeom::EigenSystem& prepareEigenSystem(const MeshLaplacian& laplaceMat, int eigenCount);

	double	calHK(int v1, int v2, double timescale) const;
	void	calHeat(int vSrc, double tMultiplier, std::vector<double>& vHeat);
	double	calHeatTrace(double timescale) const;
	double	calBiharmonic(int v1, int v2) const;
	double	calMHW(int v1, int v2, double timescale) const;
	void	computeHeatDiffuseMat(double tMultiplier);
	void	calKernelSignature(double scale, KernelType kernelType, std::vector<double>& values) const;
	void	computeKernelSignatureFeatures(const std::vector<double>& timescales, KernelType kernelType);
	ZGeom::SparseSymMatVecSolver& getHeatSolver() { return mHeatDiffuseSolver; }

	/* editing */
	void addNewHandle(int hIdx);
	int getActiveHandle() const { return mActiveHandle; }
	void setActiveHandle(int h) { mActiveHandle = h; }
	std::map<int, ZGeom::Vec3d>& getHandles() { return mHandles; }
	const std::map<int, ZGeom::Vec3d>& getHandles() const { return mHandles; }
	void clearAllHandles() { mHandles.clear(); }
	int  getRefPointIndex() const { return mRefVert; }
	void setRefPointIndex(int i) { mRefVert = i; }
	void setRefPointPosition(int x, int y, int z) { mRefPos = ZGeom::Vec3d(x, y, z); }

private:
	CMesh* mMesh;
	CMesh* mOriMesh;

	int mRefVert;
	ZGeom::Vec3d mRefPos;
	std::map<int, ZGeom::Vec3d> mHandles;
	int mActiveHandle;

	MeshLaplacian mMeshLaplacians[LaplacianTypeCount];
	ZGeom::EigenSystem mMHBs[LaplacianTypeCount];
	ZGeom::SparseMatrix<double> mHeatDiffuseMat;
	ZGeom::SparseSymMatVecSolver mHeatDiffuseSolver;
};

