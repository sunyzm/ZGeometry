#ifndef ZGEOM_GEODESIC_H
#define ZGEOM_GEODESIC_H
#include "Mesh.h"

namespace ZGeom {
    double calGeodesic(const CMesh& mesh, int s, int t);
    double calGeodesicToBoundary(CMesh& mesh, int s);
}
#endif