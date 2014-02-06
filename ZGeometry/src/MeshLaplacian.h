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
	LaplacianType m_laplacianType;

	typedef void (MeshLaplacian::*MeshLaplacianConstruct)(const CMesh* tmesh);	

	MeshLaplacian() : Laplacian() 
	{ 
		mConstructFunc[Tutte] = &MeshLaplacian::constructTutte; 
		mConstructFunc[Umbrella] = &MeshLaplacian::constructUmbrella;
		mConstructFunc[NormalizedUmbrella] = &MeshLaplacian::constructNormalizedUmbrella;
		mConstructFunc[CotFormula] = &MeshLaplacian::constructCotFormula;
		mConstructFunc[SymCot] = &MeshLaplacian::constructSymCot;
	}
	virtual ~MeshLaplacian() {}

	void constructUmbrella(const CMesh* tmesh);				// negative symmetric graph Laplacian. L = A - D
	void constructNormalizedUmbrella(const CMesh* tmesh);	// negative symmetric, normalized graph Laplacian; L = D^(-1/2) * (A-D) * D^(-1/2) = D^(-1/2) * A * D^(-1/2) - I
	void constructTutte(const CMesh* tmesh);				// negative asymmetric graph Laplacian; random walk. L = D^(-1) * (A-D) = D^(-1)*A - I
	void constructCotFormula(const CMesh* tmesh);			// negative cotangent formula
	void constructSymCot(const CMesh* tmesh);				// negative symmetric cotangent formula

	void constructAnisotropic1(const CMesh* tmesh, int nRing, double hPara1, double hPara2);
	void constructFromMesh4(const CMesh* tmesh, int nRing, double hPara1, double hPara2);
	void constructFromMesh5(const CMesh* tmesh);		

	MeshLaplacianConstruct getConstructFunc(LaplacianType laplacianType) { return mConstructFunc[laplacianType]; }

private:
	MeshLaplacianConstruct mConstructFunc[LaplacianTypeCount];
};