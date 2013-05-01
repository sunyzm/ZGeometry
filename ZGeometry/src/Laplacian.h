#pragma once
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include <ZMesh.h>

enum LaplacianType {Umbrella = 0, CotFormula, Anisotropic1,
	                Anisotropic2, Anisotropic3, IsoApproximate, LaplacianTypeCount};

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
	void write(const std::string& meshPath, bool binaryMode = true) const;
	void read(const std::string& meshPath, bool binaryMode = true);
	MeshFunction getManifoldHarmonic(int k) const;
	void dumpEigenValues(const std::string& evlPath) const;
	bool empty() const { return m_func.empty(); }
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
	Laplacian(){}
	void constructFromMesh(const CMesh* tmesh);
};

class ManifoldLaplaceHarmonics : public ManifoldHarmonics
{
public:
	bool decompLaplacian(Engine *ep, const CMesh *tmesh, int nEigFunc);
};

class MeshLaplacian : public SparseMeshMatrix
{
public:
	MeshLaplacian() {}
	void constructFromMesh1(const CMesh* tmesh);	// graph Laplacian
	void constructFromMesh2(const CMesh* tmesh);	// Cot formula
	void constructFromMesh3(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh5(const CMesh* tmesh);	

	LaplacianType m_laplacianType;
};