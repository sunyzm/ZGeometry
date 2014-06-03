#ifndef ZGEOM_MCA_H
#define ZGEOM_MCA_H
#include <vector>
#include "VecN.h"
#include "SparseRepresentation.h"

namespace ZGeom {

struct MCAoptions
{
	int nIter;
};

void singleChannelMCA(const VecNd& vSignal, const std::vector<Dictionary>& dicts, MCAoptions* opts);

void multiChannelMCA(const std::vector<Dictionary>& dicts);

}

#endif