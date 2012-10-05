#pragma once
#include <mesh/Mesh.h>
#include <engine.h>
#include <mesh/Laplacian.h>
#include <string>
#include <vector>


class ManifoldFunction
{
public:
	int m_size;
	std::vector<double> m_function;
	
	ManifoldFunction() : m_size(0) {}
	ManifoldFunction(int s) : m_size(s) {}
};

class SparseManifoldFunction
{
public:
	int size;
	std::vector<int> m_funcIdx;
	std::vector<double> m_funcVal;
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

	void normalizeFrom(const std::vector<double>& vFrom);
	void logNormalizeFrom(const std::vector<double>& vFrom);
	void bandCurveFrom(const std::vector<double>& vFrom, double lowend, double highend);
	
	
	CMesh* mesh;
	Engine *ep;
	Laplacian mLaplacian;
	ManifoldHarmonics mhb;
	std::vector<double> vDisplaySignature;
	int pRef;
	bool isMHBuilt;
	int m_size;
};

