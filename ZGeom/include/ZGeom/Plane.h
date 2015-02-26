#ifndef ZGEOM_PLANE_H
#define ZGEOM_PLANE_H
#include "Vec3.h"

namespace ZGeom {

class Plane3
{
public:
    // Construction and destruction.  The default constructor sets the normal
    // to (0,...,0,1) and the constant to zero (plane z = 0).
    Plane3() : normal(0, 0, 1), constant(0) {}

    // Specify U and c directly.
    Plane3(const Vec3d& inNormal, double inConstant) : normal(inNormal), constant(inConstant) {}

    // U is specified, c = Dot(U,p) where p is a point on the hyperplane.
    Plane3(const Vec3d& inNormal, const Vec3d& p) : normal(inNormal), constant(p.dot(inNormal)) {}

    // Public member access.
    Vec3d normal;
    double constant;
};

}

#endif