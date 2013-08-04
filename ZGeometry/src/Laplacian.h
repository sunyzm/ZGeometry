#pragma once
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include <ZMesh/ZMesh.h>
#include <ZGeom/SparseMatrix.h>

enum LaplacianType {Umbrella = 0, CotFormula = 1, Anisotropic1 = 2,
	                Anisotropic2 = 3, IsoApproximate = 4, LaplacianTypeCount};

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

class ManifoldLaplaceHarmonics : public ManifoldHarmonics
{
public:
    bool decompLaplacian(Engine *ep, const CMesh *tmesh, int nEigFunc);
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
	const std::vector<double>& getVerticesWeight() const { return vWeights; };
	int getNonzeroNum() const { return vSS.size(); }
	void getSparseMatrix(std::vector<int>& II, std::vector<int>& JJ, std::vector<double>& SS) const;
	void dumpMatrix(const std::string& path) const;

	virtual void constructFromMesh(const CMesh* tmesh) {}
	virtual void decompose(ManifoldHarmonics& mhb, int nEig, Engine *ep) const;
};



class MeshLaplacian : public SparseMeshMatrix
{
public:
    static const std::string LaplacianTypeNames[];
    LaplacianType m_laplacianType;

	MeshLaplacian() {}
	void constructFromMesh1(const CMesh* tmesh);	// graph Laplacian
	void constructFromMesh2(const CMesh* tmesh);	// Cot formula
	void constructFromMesh3(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh5(const CMesh* tmesh);		
};