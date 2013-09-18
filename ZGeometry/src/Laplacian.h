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

class ManifoldHarmonics : public ZGeom::EigenSystem
{
};

class MeshLaplacian
{
public:
	static const std::string LaplacianTypeNames[];
	enum LaplacianType {Tutte = 0, Umbrella,  CotFormula, SymCot, Anisotropic1,
						Anisotropic2, IsoApproximate, 
						LaplacianTypeCount} m_laplacianType;
		

	MeshLaplacian() : mLaplacianConstructed(false) {}
	void constructUmbrella(const CMesh* tmesh);		// symmetric graph Laplacian 
	void constructTutte(const CMesh* tmesh);		// asymmetric graph Laplacian
	void constructCotFormula(const CMesh* tmesh);	// Cotangent formula
	void constructSymCot(const CMesh* tmesh);		// symmetric cotangent formula
	void constructFromMesh3(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2);
	void constructFromMesh5(const CMesh* tmesh);		

	bool isLaplacianConstructed() const { return mLaplacianConstructed; }
	const ZGeom::SparseMatrix<double>& getLS() const { return mLS; }
	const ZGeom::SparseMatrix<double>& getW() const { return mW; }
	void decompose(int nEig, ZGeom::MatlabEngineWrapper* ep, ZGeom::EigenSystem& eigSys);

private:
	int mOrder;
	ZGeom::SparseMatrix<double> mLS, mW;
	bool mLaplacianConstructed;
};