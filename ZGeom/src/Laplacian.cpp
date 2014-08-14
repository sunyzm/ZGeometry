#include "Laplacian.h"
#include "EigenCompute.h"
#include "zassert.h"

namespace ZGeom
{

void Laplacian::decompose( int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, bool generalized /*= true*/ )
{
	runtime_assert(ep->isOpened(), "Matlab engine not opened!");
	EigenCompute eigenCompute(ep);
	if (generalized) {
		std::cout << "Do generalized eigendecomposition!\n";
		eigenCompute.solveGenSym(mLS, mW, nEig, eigSys);
	} else {
		std::cout << "Do standard eigendecomposition!\n";
		eigenCompute.solveStdSym(mLS, nEig, eigSys);
	}
}

void Laplacian::decomposeGeneralized(int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, const SparseMatrix<double>& matB)
{
	runtime_assert(ep->isOpened(), "Matlab engine not opened!");
	EigenCompute eigenCompute(ep);
	std::cout << "Do generalized eigendecomposition!\n";
	eigenCompute.solveGenSym(mLS, matB, nEig, eigSys);
}

void Laplacian::computeSubLaplacian( const std::vector<int>& vSelected, Laplacian& subLaplacian ) const
{
	const int subSize = (int)vSelected.size();
	subLaplacian.mOrder = subSize;
	mW.computeSubMatrix(vSelected, subLaplacian.mW);
	mLS.computeSubMatrix(vSelected, subLaplacian.mLS);

	subLaplacian.mLS.makeLaplacian();
}

}// end of namespace