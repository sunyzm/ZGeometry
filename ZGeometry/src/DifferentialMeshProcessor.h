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
enum KernelType {HEAT_KERNEL, MEXICAN_HAT_KERNEL};

#define SIGNATURE_HKS					0x0101
#define SIGNATURE_HK					0x0102
#define SIGNATURE_MEAN_CURVATURE		0x0103
#define SIGNATURE_GAUSS_CURVATURE		0x0104
#define SIGANTURE_WKS					0x0105

#define FEATURE_NEIGHBORS				0x0201

double transferScalingFunc1(double lambda);

double transferFunc1(double lambda, double t);	// Mexican-hat square
double transferFunc2(double lambda, double t);	// Mexican-hat
double transferFunc3(double lambda, double t);  // heat kernel
double transferFunc4(double lambda, double t);	// Mexican-hat (v2)

class DifferentialMeshProcessor : public MeshProcessor
{
public:
	DifferentialMeshProcessor(void);
	~DifferentialMeshProcessor(void);
	void init(CMesh* tm, Engine* e);
	void decomposeLaplacian(int nEigFunc);
	void readMHB(const std::string& path);
	void writeMHB(std::string path);
	void addNewHandle(int hIdx);
	void computeCurvature(std::vector<double>& vCurvature, int curvatureType = 0); //0: mean; 1: Gauss
	void deform(const std::vector<int>& vHandleIdx, const std::vector<Vector3D>& vHanldelPos, const std::vector<int>& vFreeIdx, std::vector<Vector3D>& vDeformedPos, DeformType dfType);
	void computeKernelSignature(double timescale, KernelType kernelType);
	void computeKernelDistanceSignature(double timescale, KernelType kernelType, int refPoint);
	void computeSGW(const std::vector<double>& timescales, double (*transferWavelet)(double, double) = &transferFunc1, bool withScaling = false, double (*transferScaling)(double) = &transferScalingFunc1);
	void computeMexicanHatWavelet(std::vector<double>& vMHW, double scale, int wtype = 1);
	void computeExperimentalWavelet(std::vector<double>& vExp, double scale);
	void computeDWTCoefficient(std::vector<double>& vCoeff, const std::vector<double>& vScales, const std::vector<double>& vfunc);
	void calGeometryDWT();
	void reconstructExperimental1(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false);
	void filterBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz);
	void reconstructByDifferential(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructByMHB(int approxN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
	void getSGWSignature(double timescale, std::vector<double>& values) const;
	
	bool isSGWComputed() const { return m_bSGWComputed; }
	bool isLaplacianDecomposed() const { return  m_bLaplacianDecomposed; }

	// ---- attribute access --- //
	const ManifoldHarmonics& getMHB() const { return mhb; }
	int getRefPointIndex() const { return pRef; }
	void setRefPointIndex(int i) { pRef = i; }
	void setRefPointPosition(int x, int y, int z) { posRef = Vector3D(x, y, z); }

public:
	int active_handle;
	std::map<int, Vector3D> mHandles;
	std::vector<MeshFeature> vFeatures;	
	double constrain_weight;
	
private:
	Engine *m_ep;
	MatlabWrapper matlabWrapper;

	bool m_bSGWComputed;
	bool m_bLaplacianDecomposed;		// mhb available

	Laplacian mLaplacian;
	ManifoldHarmonics mhb;
	
	int pRef;
	Vector3D posRef;

	std::vector<double> m_vTimescales;
	std::vector<std::vector<double> > m_vSGW;
};

