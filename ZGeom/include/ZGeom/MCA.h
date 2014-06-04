#ifndef ZGEOM_MCA_H
#define ZGEOM_MCA_H
#include <vector>
#include "VecN.h"
#include "SparseRepresentation.h"

namespace ZGeom {
	
struct MCAoptions
{
	enum ThresholdingMode { SOFT_THRESH, HARD_THRESH };
	enum ThresholdingStrategy {MEAN_OF_MAX, SECOND_TO_MAX, MIN_OF_MAX, MAX_OF_MAX};
	
	ThresholdingMode threshMode;
	ThresholdingStrategy threshStrategy;
	int nIter;

	MCAoptions() : nIter(500), threshMode(HARD_THRESH), threshStrategy(MEAN_OF_MAX) {}
};

void singleChannelMCA(const VecNd& vSignal, const std::vector<const Dictionary*>& vDicts, std::vector<SparseCoding>& vCodings, MCAoptions* opts);

void multiChannelMCA(const std::vector<Dictionary>& dicts);

}

#endif