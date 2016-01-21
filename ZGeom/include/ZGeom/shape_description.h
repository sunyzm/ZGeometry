#ifndef SHAPE_DESCRIPTION_H
#define SHAPE_DESCRIPTION_H
#include "mesh_primitives.h"

namespace ZGeom {

class RegionDescriptor
{
public:
    int num_contours;
    int num_scales;
    std::vector<double> coeff_stat;

    RegionDescriptor() {}
    RegionDescriptor(const RegionDescriptor& rd) = default;
    RegionDescriptor(RegionDescriptor&& rd) { *this = std::move(rd); }
};

RegionDescriptor computeCoeffStatRegionDescriptor(const BandedMeshRegions& banded_region, const std::vector<std::vector<double>>& signal_coeff);

} // end of namespace

#endif