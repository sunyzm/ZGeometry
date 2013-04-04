#pragma once
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include <ZMesh.h>

#define LBO_EPS 0
#define LBO_COT 1
#define LBO_UMBRELLA 2

// const double PZERO  = 0.00001;
// const double NZERO  = -0.00001;
// const int HOLE_SIZE = 10;

class ManifoldHarmonics;

class MeshMatrix
{
public:
	int m_size;
	std::vector<double> vWeights;
	std::vector<int> vII, vJJ;
	std::vector<double> vSS;

public:
	MeshMatrix() : m_size(0) {}
	double innerProduct(const std::vector<double>& vf, const std::vector<double>& vg) const;
	const std::vector<double>& getVerticesWeight() const { return vWeights; };
	void multiply(Engine *ep, const std::vector<double>& func, std::vector<double>& result) const;
	int getNonzeroNum() const { return vSS.size(); }
};

class Laplacian : public MeshMatrix
{
public:
	enum LaplacianType {Umbrella, CotFormula};
	bool isBuilt;

public:
	Laplacian() : isBuilt(false) {}
	void computeLaplacian(const CMesh* tmesh, LaplacianType laplacianType = CotFormula);
	void decompose(ManifoldHarmonics& mhb, int nEig, Engine *ep) const;
	void getSparseLaplacian(std::vector<int>& II, std::vector<int>& JJ, std::vector<double>& SS) const;
	void dumpLaplacian(const std::string& path) const;
};

class ManifoldBasis
{
public:
	std::vector<double> m_vec;	// eigenvector
	double m_val;			    // eigenvalue
};

class ManifoldHarmonics
{
public:
	ManifoldHarmonics() : m_size(0), m_nEigFunc(0) {}
	bool decompLaplacian(Engine *ep, const CMesh *tmesh, int nEigFunc, Laplacian::LaplacianType lbo_type = Laplacian::CotFormula);
	void decompLaplacian(Engine *ep, const Laplacian& matLaplacian, int nEigFunc) { matLaplacian.decompose(*this, nEigFunc, ep); }
	void write(const std::string& meshPath) const;
	void read(const std::string& meshPath);
	MeshFunction getManifoldHarmonic(int k) const;
	void dumpEigenValues(const std::string& evlPath) const;
public:
	std::vector<ManifoldBasis> m_func;	// manifold harmonic basis
	int m_size;	    // shape size
	int m_nEigFunc; // number of basis < shape size
};
