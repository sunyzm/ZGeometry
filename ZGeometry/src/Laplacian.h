#pragma once
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include <ZMesh.h>

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
	void write(const std::string& meshPath) const;
	void read(const std::string& meshPath);
	MeshFunction getManifoldHarmonic(int k) const;
	void dumpEigenValues(const std::string& evlPath) const;
public:
	std::vector<ManifoldBasis> m_func;	// manifold harmonic basis
	int m_size;	    // shape size
	int m_nEigFunc; // number of basis < shape size
};

class SparseMeshMatrix
{
public:
	int m_size;
	bool m_bMatrixBuilt;
	std::vector<double> vWeights;	// weight associated with each vertex
	std::vector<int> vII, vJJ;		// 1-based index of non-zero elements
	std::vector<double> vSS;		// value of non-zero elements

public:
	SparseMeshMatrix() : m_size(0), m_bMatrixBuilt(false) {}
	double innerProduct(const std::vector<double>& vf, const std::vector<double>& vg) const;
	const std::vector<double>& getVerticesWeight() const { return vWeights; };
	void multiply(Engine *ep, const std::vector<double>& func, std::vector<double>& result) const;
	int getNonzeroNum() const { return vSS.size(); }
	void getSparseMatrix(std::vector<int>& II, std::vector<int>& JJ, std::vector<double>& SS) const;
	void dumpMatrix(const std::string& path) const;

	virtual void constructFromMesh(const CMesh* tmesh) {}
	virtual void decompose(ManifoldHarmonics& mhb, int nEig, Engine *ep) const;

};

class Laplacian : public SparseMeshMatrix
{
public:
	enum LaplacianType {Umbrella, CotFormula} m_laplacianType;

public:
	Laplacian() : m_laplacianType(CotFormula) {}
	void setLaplacianType(LaplacianType lt) { m_laplacianType = lt; }
	void constructFromMesh(const CMesh* tmesh);
};

class ManifoldLaplaceHarmonics : public ManifoldHarmonics
{
public:
	bool decompLaplacian(Engine *ep, const CMesh *tmesh, int nEigFunc, Laplacian::LaplacianType lbo_type = Laplacian::CotFormula);
};

class AnisotropicLaplacian : public SparseMeshMatrix
{
public:
	AnisotropicLaplacian() : anisotropicType(0) {}
	void constructFromMesh1(const CMesh* tmesh);
	void constructFromMesh2(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	int anisotropicType; 
	void setAnisotropicType(int a) { anisotropicType = a; }
};