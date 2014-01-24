#include "Laplacian.h"
#include "EigenCompute.h"

namespace ZGeom
{
	void Laplacian::decompose( int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, bool generalized /*= true*/ )
	{
		EigenCompute eigenCompute(ep);
		if (generalized) {
			std::cout << "Do generalized eigendecomposition!\n";
			eigenCompute.solveGenSym(mLS, mW, nEig, eigSys);
		} else {
			std::cout << "Do standard eigendecomposition!\n";
			eigenCompute.solveStdSym(mLS, nEig, eigSys);
		}
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