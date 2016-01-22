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
    std::vector<double> perimeter_stat;
    std::vector<double> distance_stat;


    RegionDescriptor() {}
    RegionDescriptor(const RegionDescriptor& rd) = default;
    RegionDescriptor(RegionDescriptor&& rd) { *this = std::move(rd); }
};

void computeRegionDescriptorCoeffStat(RegionDescriptor& descriptor, const BandedMeshRegion& banded_region, const std::vector<std::vector<double>>& signal_coeff);
void computeRegionDescriptorContourStat(RegionDescriptor& descriptor, const BandedMeshRegion& banded_region, const CMesh& mesh);

} // end of namespace

#endif