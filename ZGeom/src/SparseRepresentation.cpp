#include "SparseRepresentation.h"

namespace ZGeom{

ZGeom::VecNd ReconstructApproximationSingleChannel(const Dictionary& dict, const SparseCoding& approx)
{
	ZGeom::VecNd vApprox(dict.atomDim(), 0);
	for (auto item : approx.getApproxItems())
	{
		vApprox += dict[item.index()] * item.coeff();
	}

	return vApprox;
}

ZGeom::VecNd ReconstructApproximationSingleChannel(const std::vector<ZGeom::VecNd>& vAtoms, const SparseCoding& approx)
{
	ZGeom::VecNd vApprox(vAtoms[0].size());
	for (auto item : approx.getApproxItems())
	{
		vApprox += vApprox[item.index()] * item.coeff();
	}

	return vApprox;
}

void ReconstructApproximationMultiChannel(const std::vector<VecNd>& vAtoms, const std::vector<SparseCoding>& vApprox, std::vector<VecNd>& vReconstructed)
{
	int nChannels = (int)vApprox.size();
	vReconstructed.resize(nChannels);
	for (int i = 0; i < nChannels; ++i)
		vReconstructed[i] = ReconstructApproximationSingleChannel(vAtoms, vApprox[i]);
}

} // end of namespace