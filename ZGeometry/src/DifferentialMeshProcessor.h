#pragma once
#include <engine.h>
#include <ZMesh.h>
#include <string>
#include <vector>
#include <map>
#include "Laplacian.h"

enum DeformType {Simple, Shell, Differential};

class DifferentialMeshProcessor : public MeshProcessor
{
public:
	DifferentialMeshProcessor(void);
	~DifferentialMeshProcessor(void);
	void init(CMesh* tm, Engine* e);
	void decomposeLaplacian(int nEigFunc);
	void readMHB(const std::string& path);
	void writeMHB(std::string path);
	void computeCurvature(std::vector<double>& vCurvature, int curvatureType = 0); //0: mean; 1: Gauss

	void normalizeSignatureFrom(const std::vector<double>& vFrom);
	void logNormalizeSignatureFrom(const std::vector<double>& vFrom);
	void bandCurveSignatureFrom(const std::vector<double>& vFrom, double lowend, double highend);

public:
	void deform(const std::vector<int>& vHandleIdx, const std::vector<Vector3D>& vHanldelPos, const std::vector<int>& vFreeIdx, std::vector<Vector3D>& vDeformedPos, DeformType dfType);

	Engine *m_ep;
	
	Laplacian mLaplacian;
	ManifoldHarmonics mhb;
	bool isMHBBuilt;

	std::vector<double> vDisplaySignature;
	double sigMin, sigMax;

	int pRef;
	Vector3D posRef;

	std::map<int, Vector3D> mHandles;
	int active_handle;
	std::vector<MeshFeature> vFeatures;	
	void addNewHandle(int hIdx);
	
	/// spectral wavelet transform related
private:
	std::vector<double> m_vTimescales;
	std::vector<std::vector<double> > m_vSGW;
public:
	void computeMexicanHatWavelet(std::vector<double>& vMHW, double scale, int wtype = 1);
	void computeExperimentalWavelet(std::vector<double>& vExp, double scale);
	void computeDWTCoefficient(std::vector<double>& vCoeff, const std::vector<double>& vScales, const std::vector<double>& vfunc);
	void calGeometryDWT();
	void reconstructExperimental1(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false);
	void filterBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz);
	void reconstructByDifferential(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructByMHB(int approxN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
	void computeSGW(const std::vector<double>& timescales);
	void getSGWSignature(double timescale, std::vector<double>& values) const;
};

