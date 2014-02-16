#pragma once
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <engine.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/Mesh.h>
#include <ZGeom/MeshProcessor.h>
#include "MeshLaplacian.h"
#include "global.h"

enum KernelType {HEAT_KERNEL, MHW_KERNEL, SGW_KERNEL, BIHARMONIC_KERNEL};
enum DistanceType {DISTANCE_GEODESIC, DISTANCE_BIHARMONIC, DISTANCE_HK, DISTANCE_MHW};
enum PointSimilarityType {SIM_TYPE_1, SIM_TYPE_2, SIM_TYPE_3, SIM_TYPE_COUNT};

enum FeatureID { FEATURE_ID	= 0x0200, FEATURE_NEIGHBORS, FEATURE_HKS,
				 FEATURE_MHWS, FEATURE_SGWS, FEATURE_MULTI_HKS, FEATURE_DEMO, FEATURE_DEMO2,
				 FEATURE_SGW_SOMP,
				 FEATURE_ID_COUNT};

double transferScalingFunc1(double lambda);

double transferFunc1(double lambda, double t);	// Mexican-hat square wavelet
double transferFunc2(double lambda, double t);	// Mexican-hat wavelet
double heatKernelTransferFunc(double lambda, double t);  // heat kernel
double mhwTransferFunc1(double lambda, double t);	// Mexican-hat wavelet (Tingbo)

typedef double (*TransferFunc)(double, double);
typedef double (*ScalelessTransferFunc)(double);

class DifferentialMeshProcessor : public MeshProcessor
{
public:
	DifferentialMeshProcessor();
	DifferentialMeshProcessor(CMesh* tm, CMesh* originalMesh);
	~DifferentialMeshProcessor();

	void init(CMesh* tm);
	void init_lite(CMesh* tm, CMesh* originalMesh);
	
	// ---- computation ---- //
	void	constructLaplacian(LaplacianType laplacianType = CotFormula);
	void	decomposeLaplacian(int nEigFunc, LaplacianType laplacianType = CotFormula);
	void	computeCurvature(std::vector<double>& vCurvature, int curvatureType = 0); //0: mean; 1: Gauss
	double	calHK(int v1, int v2, double timescale) const;
	void	calHeat(int vSrc, double tMultiplier, std::vector<double>& vHeat);
	double	calHeatTrace(double timescale) const;
	double	calBiharmonic(int v1, int v2) const;
	double	calMHW(int v1, int v2, double timescale) const;
	double	calSGW(int v1, int v2, double timescale, LaplacianType laplacianType) const;	
	void	computeSGW1(LaplacianType laplacianType = CotFormula);
	void	computeSGW2(LaplacianType laplacianType = CotFormula);
	void    computeMixedAtoms1(LaplacianType laplacianType = CotFormula);
	void	computeHeatKernelMat(double t, ZGeom::DenseMatrix<double>& hkmat);
	void	computeHeatKernelMat_AMP(double t, ZGeom::DenseMatrix<double>& hkmat);
	void	computeHeatDiffuseMat(double tMultiplier);

	void	calKernelSignature(double scale, KernelType kernelType, std::vector<double>& values) const;
	void	computeKernelSignatureFeatures(const std::vector<double>& timescales, KernelType kernelType);
	void	computeSimilarityMap1(int refPoint);
	void	computeSimilarityMap2(int refPoint);
	void	computeSimilarityMap3(int refPoint);

	// ---- editing ---- //
	void addNewHandle(int hIdx);
	int getActiveHandle() const { return mActiveHandle; }
	void setActiveHandle(int h) { mActiveHandle = h; }
	std::map<int, Vector3D>& getHandles() { return mHandles; }
	const std::map<int, Vector3D>& getHandles() const { return mHandles; }
	void clearAllHandles() { mHandles.clear(); }

	// ---- attribute access --- //
	const MeshLaplacian& getMeshLaplacian(LaplacianType laplacianType) const { return mMeshLaplacians[laplacianType]; }
	bool hasLaplacian(LaplacianType laplacianType) { return mMeshLaplacians[laplacianType].isLaplacianConstructed(); }
	bool isLaplacianDecomposed(LaplacianType laplacianType) { return !mMHBs[laplacianType].empty(); }

	std::string generateMHBPath(const std::string& prefix, LaplacianType laplacianType);
	void loadMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	void saveMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	const ZGeom::EigenSystem& getMHB(LaplacianType laplacianType) const { return mMHBs[laplacianType]; }
	const ZGeom::DenseMatrixd& getWaveletMat() const { return mMatAtoms; }
	ZGeom::DenseMatrixd& getWaveletMat() { return mMatAtoms; }
	ZGeom::SparseSymMatVecSolver& getHeatSolver() { return mHeatDiffuseSolver; }

	int  getRefPointIndex() const { return mRefVert; }
	void setRefPointIndex(int i) { mRefVert = i; }
	void setRefPointPosition(int x, int y, int z) { mRefPos = Vector3D(x, y, z); }
	const MeshFeatureList* getActiveFeatures() const;
	void setActiveFeaturesByID(int feature_id) { mActiveFeature = feature_id; }
	bool isMHBCacheValid( const std::string& pathMHB, int eigenCount );

	static void computeGeometricLaplacianCoordinate(const CMesh& mesh, const MeshCoordinates& eCoord, MeshCoordinates& lCoord);
	static void computeSGWMat1(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW);
	static void computeSGWMat2(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW);

private:
	int mRefVert;
	Vector3D mRefPos;
	std::map<int, Vector3D> mHandles;
	int mActiveHandle;
	int mActiveFeature;

	MeshLaplacian mMeshLaplacians[LaplacianTypeCount];
	ZGeom::EigenSystem mMHBs[LaplacianTypeCount];

	ZGeom::DenseMatrixd mMatAtoms;
	ZGeom::SparseMatrix<double> mHeatDiffuseMat;
	ZGeom::SparseSymMatVecSolver mHeatDiffuseSolver;
};

