#include "Laplacian.h"
#include "EigenCompute.h"

namespace ZGeom
{
	void Laplacian::decompose( int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys )
	{
		EigenCompute eigenCompute(ep);
		eigenCompute.solveGenSym(mLS, mW, nEig, eigSys);
	}

	void Laplacian::computeSubLaplacian( const std::vector<int>& vSelected, Laplacian& subLaplacian ) const
	{
		const int subSize = (int)vSelected.size();
		subLaplacian.mOrder = subSize;


	}

}// end of namespace