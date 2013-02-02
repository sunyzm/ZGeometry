#pragma once
#include <engine.h>
#include <ZMesh.h>
#include <string>
#include <vector>
#include <map>
#include "Laplacian.h"

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
};

