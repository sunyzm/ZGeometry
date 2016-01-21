#include "shape_description.h"
#include <cmath>

namespace ZGeom {

RegionDescriptor computeCoeffStatRegionDescriptor(const BandedMeshRegions& banded_region, const std::vector<std::vector<double>>& signal_coeff)
{
    RegionDescriptor descriptor;
    const int num_of_contours = banded_region.numOfRegions();
    const int num_of_scales = (int)signal_coeff.size();;
    descriptor.num_contours = num_of_contours;
    descriptor.num_scales = num_of_scales;
    descriptor.coeff_stat.resize(num_of_contours * num_of_scales);

    for (int s = 0; s < num_of_scales; ++s) {
        for (int c = 0; c < num_of_contours; ++c) {
            double sum(0);
            for (int vi : banded_region.band_verts[c]) {
                sum += std::fabs(signal_coeff[s][vi]);
            }
            descriptor.coeff_stat[num_of_contours * s + c] = sum;
        }
    }
    return descriptor;
}

} // end of namespace