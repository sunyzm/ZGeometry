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

class Laplacian
{
public:
	enum LaplacianType {Umbrella, CotFormula};
	void computeLaplacian(const CMesh* tmesh, LaplacianType laplacianType = CotFormula);
	void decompose(ManifoldHarmonics& mhb, int nEig, Engine *ep) const;
	double innerProduct(const std::vector<double>& vf, const std::vector<double>& vg) const;
	const std::vector<double>& getVerticesWeight() const { return vWeights; };
	void multiply(Engine *ep, const std::vector<double>& func, std::vector<double>& result) const;
	Laplacian() : isBuilt(false), size(0) {}
	void getSparseLaplacian(std::vector<int>& II, std::vector<int>& JJ, std::vector<double>& SS) const;
private:
	int size;
	std::vector<double> vWeights;
	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	bool isBuilt;
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
	void write(const std::string& meshPath) const;
	void read(const std::string& meshPath);
	MeshFunction getManifoldHarmonic(int k) const;
public:
	std::vector<ManifoldBasis> m_func;	// manifold harmonic basis
	int m_size;	    // shape size
	int m_nEigFunc; // number of basis < shape size
};
