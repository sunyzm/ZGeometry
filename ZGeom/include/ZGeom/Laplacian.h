#ifndef ZGEOM_LAPLACIAN_H
#define ZGEOM_LAPLACIAN_H
#include <vector>
#include "SparseMatrix.h"
#include "MatlabEngineWrapper.h"
#include "EigenSystem.h"

namespace ZGeom
{

class Laplacian
{
public:
	Laplacian() : mOrder(0), mConstructed(false) {}
	virtual ~Laplacian(){}

	const ZGeom::SparseMatrix<double>& getLS() const { return mLS; }
	const ZGeom::SparseMatrix<double>& getW() const { return mW; }
	void decompose(int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys);
	bool isLaplacianConstructed() const { return mConstructed; }
	void computeSubLaplacian(const std::vector<int>& vSelected, Laplacian& subLaplacian) const;

protected:
	int mOrder;
	SparseMatrix<double> mLS, mW;
	bool mConstructed;
};

} // end of namespace

#endif