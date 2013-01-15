#pragma once
#include <engine.h>
#include <ZMesh.h>
#include <string>
#include <vector>

class ManifoldFunction
{
public:
	int m_size;
	std::vector<double> m_function;
	
	ManifoldFunction() : m_size(0) {}
	ManifoldFunction(int s) : m_size(s) {m_function.resize(s); }
	void setSize(int s) { m_size = s; m_function.resize(s); }
};

class SparseManifoldFunction
{
public:
	int size;
	std::vector<int> m_funcIdx;
	std::vector<double> m_funcVal;
};

class ManifoldFeature
{
public:
	int index;
	int scale;
	double value;
	int featureType;
	ManifoldFeature(int i, int s) : index(i), scale(s), featureType(0){}
};

class MeshProcessor
{
public:
	MeshProcessor(void);
	~MeshProcessor(void);
	void init(CMesh* tm, Engine* e);
	void decomposeLaplacian(int nEigFunc);
	void readMHB(const std::string& path);
	void writeMHB(std::string path);
	void computeCurvature(std::vector<double>& vCurvature, int curvatureType = 0); //0: mean; 1: Gauss

	void normalizeSignatureFrom(const std::vector<double>& vFrom);
	void logNormalizeSignatureFrom(const std::vector<double>& vFrom);
	void bandCurveSignatureFrom(const std::vector<double>& vFrom, double lowend, double highend);
		
	CMesh* mesh;
	Engine *m_ep;
	Laplacian mLaplacian;
	ManifoldHarmonics mhb;
	std::vector<double> vDisplaySignature;
	double sigMin, sigMax;
	std::vector<ManifoldFeature> vFeatures;
	int pRef;
	Vector3D posRef;
	bool isMHBuilt;
	int m_size;	
};

