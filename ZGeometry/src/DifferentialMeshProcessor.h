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

double transferScalingFunc1(double lambda);

double transferFunc1(double lambda, double t);	// Mexican-hat square
double transferFunc2(double lambda, double t);	// Mexican-hat

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
	void normalizeSignatureFrom(const std::vector<double>& vFrom);
	void logNormalizeSignatureFrom(const std::vector<double>& vFrom);
	void bandCurveSignatureFrom(const std::vector<double>& vFrom, double lowend, double highend);
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

public:
	MatlabWrapper matlabWrapper;
	Engine *m_ep;
	
	Laplacian mLaplacian;
	ManifoldHarmonics mhb;

	std::vector<double> vDisplaySignature;
	double sigMin, sigMax;

	int pRef;
	Vector3D posRef;

	int active_handle;
	std::map<int, Vector3D> mHandles;
	std::vector<MeshFeature> vFeatures;	
	double constrain_weight;
	
private:
	bool m_bSGWComputed;
	bool m_bLaplacianDecomposed;		// mhb available
	std::vector<double> m_vTimescales;
	std::vector<std::vector<double> > m_vSGW;
};

