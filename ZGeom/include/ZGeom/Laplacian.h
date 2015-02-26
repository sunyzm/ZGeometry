#ifndef ZGEOM_LAPLACIAN_H
#define ZGEOM_LAPLACIAN_H
#include <vector>
#include "SparseMatrix.h"
#include "MatlabEngineWrapper.h"
#include "EigenSystem.h"
#include "Mesh.h"

namespace ZGeom {

class Laplacian
{
public:
	Laplacian() : mOrder(0), mSymmetric(true) {}
	virtual ~Laplacian(){}

	const ZGeom::SparseMatrix<double>& getLS() const { return mLS; }
	const ZGeom::SparseMatrix<double>& getW() const { return mW; }
	void decompose(int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, bool generalized = true) const;
	void decomposeGeneralized(int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, const SparseMatrix<double>& matB) const;
	bool isLaplacianConstructed() const { return !mLS.empty(); }
	void computeSubLaplacian(const std::vector<int>& vSelected, Laplacian& subLaplacian) const;
    ZGeom::SparseMatrix<double> getSparseMatrix() const;

    void constructUmbrella(CMesh* tmesh);				// negative symmetric graph Laplacian. L = A - D
    void constructGeometricUmbrella(CMesh *tmesh);    // negative symmetric graph Laplacian. w_ij is weighted by the inverse of length between i and j
    void constructNormalizedUmbrella(CMesh* tmesh);	// negative symmetric, normalized graph Laplacian; L = D^(-1/2) * (A-D) * D^(-1/2) = D^(-1/2) * A * D^(-1/2) - I
    void constructTutte(CMesh* tmesh);				// negative asymmetric graph Laplacian; random walk. L = D^(-1) * (A-D) = D^(-1)*A - I
    void constructCotFormula(CMesh* tmesh);			// negative cotangent formula
    void constructSymCot(CMesh* tmesh);				// negative symmetric cotangent formula

    void constructAniso(CMesh *tmesh);

public:
	int mOrder;
	SparseMatrix<double> mLS, mW;
    bool mSymmetric;
};

} // end of namespace

#endif