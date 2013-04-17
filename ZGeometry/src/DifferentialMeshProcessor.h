#pragma once
#include <engine.h>
#include <ZMesh.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "Laplacian.h"
#include "MatlabWrapper.h"

enum DeformType {Simple, Shell, Laplace, SGW};
enum KernelType {HEAT_KERNEL, MHW_KERNEL, SGW_KERNEL, BIHARMONIC_KERNEL};
enum DistanceType {DISTANCE_GEODESIC, DISTANCE_BIHARMONIC, DISTANCE_HK, DISTANCE_MHW};

#define SIGNATURE_HKS					0x0101
#define SIGNATURE_HK					0x0102
#define SIGNATURE_MEAN_CURVATURE		0x0103
#define SIGNATURE_GAUSS_CURVATURE		0x0104
#define SIGNATURE_WKS					0x0105
#define SIGNATURE_MHWS					0x0106
#define SIGNATURE_MHW					0x0106
#define SIGNATURE_SGWS					0x0107
#define SIGNATURE_SGW					0x0108
#define SIGNATURE_BIHARMONIC_DISTANCE   0x0109
#define SIGNATURE_SIMILARITY_MAP		0x0110

#define FEATURE_NEIGHBORS				0x0201
#define FEATURE_HKS						0x0202
#define FEATURE_MHWS					0x0203
#define FEATURE_SGWS					0x0204
#define FEATURE_MULTI_HKS				0x0205

double transferScalingFunc1(double lambda);

double transferFunc1(double lambda, double t);	// Mexican-hat square wavelet
double transferFunc2(double lambda, double t);	// Mexican-hat wavelet
double transferFunc3(double lambda, double t);  // heat kernel
double transferFunc4(double lambda, double t);	// Mexican-hat wavelet (Tingbo)

typedef double (*TransferFunc)(double, double);
typedef double (*ScalelessTransferFunc)(double);

class DifferentialMeshProcessor : public MeshProcessor
{
public:
	DifferentialMeshProcessor(void);
	DifferentialMeshProcessor(CMesh* tm);
	~DifferentialMeshProcessor(void);
	void init(CMesh* tm, Engine* e);
	void init_lite(CMesh* tm);
	void readMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	void writeMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	void addNewHandle(int hIdx);

	// ---- computation ---- //
	void computeLaplacian(LaplacianType laplacianType = CotFormula);
	void decomposeLaplacian(int nEigFunc, LaplacianType laplacianType = CotFormula);
	void selectMHB(LaplacianType laplacianType) { if (!vMHB[laplacianType].empty()) mhb = vMHB[laplacianType]; };
	void computeCurvature(std::vector<double>& vCurvature, int curvatureType = 0); //0: mean; 1: Gauss
	void calKernelSignature(double timescale, KernelType kernelType, std::vector<double>& values) const;
	void computeKernelSignature(double timescale, KernelType kernelType);
	void computeKernelSignatureFeatures(const std::vector<double>& timescales, KernelType kernelType);
	void computeKernelDistanceSignature(double timescale, KernelType kernelType, int refPoint);
	void computeBiharmonicDistanceSignature(int refPoint);
	void computeSimilarityMap(int refPoint);
	void computeSimilarityMap2(int refPoint);
	double calHK(int v1, int v2, double timescale) const;
	double calBiharmonic(int v1, int v2) const;
	double getVertexHKS(int index, double timescale) const { return calHK(index, index, timescale); }

	void computeSGW(const std::vector<double>& timescales, TransferFunc transferWavelet = &transferFunc1, bool withScaling = false, double (*transferScaling)(double) = &transferScalingFunc1);
	void computeMexicanHatWavelet(std::vector<double>& vMHW, double scale, int wtype = 1);
	void computeExperimentalWavelet(std::vector<double>& vExp, double scale);
	void computeDWTCoefficient(std::vector<double>& vCoeff, const std::vector<double>& vScales, const std::vector<double>& vfunc);
	
	// ---- editing ---- //
	void deform(const std::vector<int>& vHandleIdx, const std::vector<Vector3D>& vHanldelPos, const std::vector<int>& vFreeIdx, std::vector<Vector3D>& vDeformedPos, DeformType dfType);
	void calGeometryDWT();
	void reconstructExperimental1(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false);
	void filterBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz);
	void reconstructByDifferential(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructByMHB(int approxN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
	
	// ---- boolean query ---- //
	bool isSGWComputed() const { return m_bSGWComputed; }
	bool isLaplacianDecomposed() const { return  m_bLaplacianDecomposed; }

	// ---- attribute access --- //
	const ManifoldHarmonics& getMHB() const { return mhb; }
	int  getRefPointIndex() const { return pRef; }
	void setRefPointIndex(int i) { pRef = i; }
	void setRefPointPosition(int x, int y, int z) { posRef = Vector3D(x, y, z); }
	const std::vector<MeshFeature*>* getActiveFeatures() const { return pvActiveFeatures; }
	void setActiveFeatures(std::vector<MeshFeature*>* pvMF) { pvActiveFeatures = pvMF; }
	void setActiveFeaturesByID(int feature_id);
	void setActiveMHB(LaplacianType laplacianType) { mhb = vMHB[laplacianType]; }
public:
	int active_handle;
	std::map<int, Vector3D> mHandles;
	double constrain_weight;
	
private:
	Engine *m_ep;
	MatlabWrapper matlabWrapper;

	bool m_bSGWComputed;
	bool m_bLaplacianDecomposed;		// mhb available

//	Laplacian mLaplacian;
//	MeshLaplacian meshKernel;
	ManifoldHarmonics mhb;	

	int	m_currentLaplacian;
public:
	MeshLaplacian vMeshLaplacian[LaplacianEnd];
	ManifoldHarmonics vMHB[LaplacianEnd];
private:	
	int pRef;
	Vector3D posRef;
	std::vector<MeshFeature*>* pvActiveFeatures;	

	std::vector<double> m_vTimescales;
	std::vector<std::vector<double> > m_vSGW;
};

