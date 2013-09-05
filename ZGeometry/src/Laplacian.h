#pragma once
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include <ZMesh/ZMesh.h>
#include <ZGeom/SparseMatrix.h>
#include <ZGeom/EigenSystem.h>
#include <ZGeom/MatlabEngineWrapper.h>

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
    int eigVecSize() const { return m_size; }
    int eigVecCount() const { return m_nEigFunc; } 

public:
	std::vector<ManifoldBasis> m_func;	// manifold harmonic basis
	int m_size;	    // shape size
	int m_nEigFunc; // number of basis < shape size
};

class MeshLaplacian
{
public:
    static const std::string LaplacianTypeNames[];
    enum LaplacianType {Umbrella = 0, CotFormula = 1, Anisotropic1 = 2,
                        Anisotropic2 = 3, IsoApproximate = 4, 
                        LaplacianTypeCount} m_laplacianType;
        

	MeshLaplacian() : mLaplacianConstructed(false) {}
	void constructFromMesh1(const CMesh* tmesh);	// graph Laplacian
	void constructFromMesh2(const CMesh* tmesh);	// Cot formula
	void constructFromMesh3(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh5(const CMesh* tmesh);		

    bool isLaplacianConstructed() const { return mLaplacianConstructed; }
    const ZGeom::SparseMatrix<double>& getLS() const { return mLS; }
    const ZGeom::SparseMatrix<double>& getW() const { return mW; }
    void decompose(ManifoldHarmonics& mhb, int nEig, Engine *ep) const;
    void decompose(int nEig, const ZGeom::MatlabEngineWrapper* ep, ZGeom::EigenSystem& eigSys);

private:
    int mOrder;
    ZGeom::SparseMatrix<double> mLS, mW;
    bool mLaplacianConstructed;
};