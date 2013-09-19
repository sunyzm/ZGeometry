#pragma once

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <engine.h>
#include <ZMesh/ZMesh.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include "Laplacian.h"
#include "MatlabWrapper.h"

enum DeformType {Simple, Shell, Laplace, SGW};
enum KernelType {HEAT_KERNEL, MHW_KERNEL, SGW_KERNEL, BIHARMONIC_KERNEL};
enum DistanceType {DISTANCE_GEODESIC, DISTANCE_BIHARMONIC, DISTANCE_HK, DISTANCE_MHW};
enum PointSimilarityType {SIM_TYPE_1, SIM_TYPE_2, SIM_TYPE_3, SIM_TYPE_COUNT};

enum SignatureID {	SIGNATURE_ID = 0x0100, SIGNATURE_EIG_FUNC, SIGNATURE_HKS, SIGNATURE_HK, 
					SIGNATURE_MEAN_CURVATURE, SIGNATURE_GAUSS_CURVATURE, SIGNATURE_WKS, 
					SIGNATURE_MHWS, SIGNATURE_MHW, SIGNATURE_SGWS, SIGNATURE_SGW, 
					SIGNATURE_BIHARMONIC_DISTANCE, SIGNATURE_SIMILARITY_MAP, 
					SIGNATURE_ID_COUNT};

enum FeatureID {	FEATURE_ID	= 0x0200, FEATURE_NEIGHBORS, FEATURE_HKS,
					FEATURE_MHWS, FEATURE_SGWS, FEATURE_MULTI_HKS, FEATURE_DEMO, FEATURE_DEMO2,
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

	void init(CMesh* tm, ZGeom::MatlabEngineWrapper* e);
	void init_lite(CMesh* tm, CMesh* originalMesh);
	void loadMHB(const std::string& path, MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	void saveMHB(const std::string& path, MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	
	// ---- computation ---- //
	void constructLaplacian(MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	void decomposeLaplacian(int nEigFunc, MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula);
	void computeCurvature(std::vector<double>& vCurvature, int curvatureType = 0); //0: mean; 1: Gauss
	void calKernelSignature(double scale, KernelType kernelType, std::vector<double>& values) const;
	void calNormalizedKernelSignature(double scale, KernelType kernelType, std::vector<double>& normalized_values) const;
	void computeKernelSignature(double timescale, KernelType kernelType);
	void computeKernelSignatureFeatures(const std::vector<double>& timescales, KernelType kernelType);
	void computeKernelDistanceSignature(double timescale, KernelType kernelType, int refPoint);
	void computeBiharmonicDistanceSignature(int refPoint);
	void computeSimilarityMap1(int refPoint);
	void computeSimilarityMap2(int refPoint);
	void computeSimilarityMap3(int refPoint);
	double calHK(int v1, int v2, double timescale) const;
	double calHeatTrace(double timescale) const;
	double calBiharmonic(int v1, int v2) const;
	double getVertexHKS(int index, double timescale) const { return calHK(index, index, timescale); }

	void computeSGW(const std::vector<double>& timescales, TransferFunc transferWavelet = &transferFunc1, bool withScaling = false, double (*transferScaling)(double) = &transferScalingFunc1);
	void computeMexicanHatWavelet(std::vector<double>& vMHW, double scale, int wtype = 1);
	void computeExperimentalWavelet(std::vector<double>& vExp, double scale);
	void computeDWTCoefficient(std::vector<double>& vCoeff, const std::vector<double>& vScales, const std::vector<double>& vfunc);
	
	// ---- editing ---- //
	void addNewHandle(int hIdx);
	void deform(const std::vector<int>& vHandleIdx, const std::vector<Vector3D>& vHanldelPos, const std::vector<int>& vFreeIdx, std::vector<Vector3D>& vDeformedPos, DeformType dfType);
	void calGeometryDWT();
	void reconstructExperimental1(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false);
	void filterBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz);
	void reconstructByDifferential(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructByMHB(int approxN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
	
	// ---- boolean query ---- //
	bool isSGWComputed() const { return m_bSGWComputed; }
	bool isLaplacianConstructed(MeshLaplacian::LaplacianType laplacianType) { return vMeshLaplacian[laplacianType].isLaplacianConstructed(); }
	bool isLaplacianDecomposed(MeshLaplacian::LaplacianType laplacianType) { return !vMHB[laplacianType].empty(); }

	// ---- attributes access --- //
	const MeshLaplacian& getMeshLaplacian(MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula) const { return vMeshLaplacian[laplacianType]; }
	const ManifoldHarmonics& getMHB(MeshLaplacian::LaplacianType laplacianType = MeshLaplacian::CotFormula) const { return vMHB[laplacianType]; }
	void setActiveLaplacian(MeshLaplacian::LaplacianType laplacianType) { mActiveLaplacianType = laplacianType;}
	int  getRefPointIndex() const { return pRef; }
	void setRefPointIndex(int i) { pRef = i; }
	void setRefPointPosition(int x, int y, int z) { posRef = Vector3D(x, y, z); }
	const MeshFeatureList* getActiveFeatures() const;
	void setActiveFeaturesByID(int feature_id);
	int getActiveHandle() const { return active_handle; }
	void setActiveHandle(int h) { active_handle = h; }
	int getConstrainWeight() const { return constrain_weight; }
	void setConstrainWeight(double w) { constrain_weight = w; }
	std::map<int, Vector3D>& getHandles() { return mHandles; }
	const std::map<int, Vector3D>& getHandles() const { return mHandles; }
		
private:
	ZGeom::MatlabEngineWrapper *mpEngineWrapper;
	MatlabWrapper matlabWrapper;
	MeshLaplacian vMeshLaplacian[MeshLaplacian::LaplacianTypeCount];
	ManifoldHarmonics vMHB[MeshLaplacian::LaplacianTypeCount];
	MeshLaplacian::LaplacianType mActiveLaplacianType;
	bool m_bSGWComputed;
	int pRef;
	Vector3D posRef;
	int active_feature_id;
	std::map<int, Vector3D> mHandles;
	int active_handle;
	double constrain_weight;
	std::vector<double> m_vTimescales;
	std::vector<std::vector<double> > m_vSGW;
};

