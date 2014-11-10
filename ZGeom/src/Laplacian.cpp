#include "Laplacian.h"
#include "EigenCompute.h"
#include "zassert.h"

namespace ZGeom {

using namespace std;

void Laplacian::decompose( int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, bool generalized /*= true*/ ) const
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

void Laplacian::decomposeGeneralized(int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, const SparseMatrix<double>& matB) const
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

    //TODO: 
	//subLaplacian.mLS.makeLaplacian();
}

ZGeom::SparseMatrix<double> Laplacian::getSparseMatrix() const
{
    SparseMatrix<double> result = getLS();
    vector<double> vWeight = getW().getDiagonal();

    for (MatElem<double> &elem : result.allElements()) {
        elem.mVal /= vWeight[elem.row() - 1];
    }
    return result;
}

}   // end of namespace