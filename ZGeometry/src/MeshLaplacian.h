#pragma once
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include <ZGeom/Laplacian.h>
#include <ZGeom/SparseMatrix.h>
#include <ZGeom/EigenSystem.h>
#include <ZGeom/Mesh.h>
#include "global.h"

class MeshLaplacian : public ZGeom::Laplacian
{
public:
	typedef void (MeshLaplacian::*MeshLaplacianConstruct)(CMesh* tmesh);	

	MeshLaplacian() : ZGeom::Laplacian()
	{ 
		mConstructFunc[Tutte] = &MeshLaplacian::constructTutte; 
		mConstructFunc[Umbrella] = &MeshLaplacian::constructUmbrella;
		mConstructFunc[NormalizedUmbrella] = &MeshLaplacian::constructNormalizedUmbrella;
		mConstructFunc[CotFormula] = &MeshLaplacian::constructCotFormula;
		mConstructFunc[SymCot] = &MeshLaplacian::constructSymCot;
		mConstructFunc[Anisotropic1] = &MeshLaplacian::constructAnisotropic1;
		mConstructFunc[Anisotropic2] = &MeshLaplacian::constructAnisotropic2;
	}
	virtual ~MeshLaplacian() {}

    void constructUmbrella(CMesh* tmesh);				// negative symmetric graph Laplacian. L = A - D
    void constructGeometricUmbrella(CMesh *tmesh);    // negative symmetric graph Laplacian. w_ij is weighted by the inverse of length between i and j
    void constructNormalizedUmbrella(CMesh* tmesh);	// negative symmetric, normalized graph Laplacian; L = D^(-1/2) * (A-D) * D^(-1/2) = D^(-1/2) * A * D^(-1/2) - I
    void constructTutte(CMesh* tmesh);				// negative asymmetric graph Laplacian; random walk. L = D^(-1) * (A-D) = D^(-1)*A - I
    void constructCotFormula(CMesh* tmesh);			// negative cotangent formula
    void constructSymCot(CMesh* tmesh);				// negative symmetric cotangent formula
    void constructAnisotropic1(CMesh* tmesh);
    void constructAnisotropic2(CMesh* tmesh);		// combine distance and curvature difference
    void constructAnisotropic3(CMesh* tmesh, int nRing, double hPara1, double hPara2);
    void constructAnisotropic4(CMesh* tmesh, int nRing, double hPara1, double hPara2);
	MeshLaplacianConstruct getConstructFunc(LaplacianType laplacianType) { return mConstructFunc[laplacianType]; }
	void meshEigenDecompose(int nEig, ZGeom::MatlabEngineWrapper *eng, ZGeom::EigenSystem& es) const;
	LaplacianType laplacianType() const { return mLaplacianType;  }

private:
	MeshLaplacianConstruct mConstructFunc[LaplacianTypeCount];
    LaplacianType mLaplacianType;
};