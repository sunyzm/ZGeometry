#ifndef ZGEOM_TRIGONOMETRY
#define ZGEOM_TRIGONOMETRY

#include <cmath>
#include "Vec3.h"

namespace ZGeom
{

template<typename T> inline
T cotVec3(const Vec3<T>& v1, const Vec3<T>& v2)
{
	return ZGeom::dot(v1, v2) / ZGeom::cross(v1, v2).length(); 
}

template<typename T> inline
T cosVec3(const Vec3<T>& v1, const Vec3<T>& v2)
{
	return ZGeom::dot(v1, v2) / (v1.length() * v2.length());
}

/// compute cosine of the angle opposite to e3
///
inline double cosTriSides(double e1, double e2, double e3)
{
	assert(e1 > 0 && e2 > 0 && e3 > 0 && 
			e1 + e2 > e3 && e1 + e3 > e2 && e2 + e3 > e1);

	return (e1*e1 + e2*e2 - e3*e3) / (2.0*e1*e2);
}

inline double triArea(double e1, double e2, double e3)
{
	assert(e1 > 0 && e2 > 0 && e3 > 0 && 
		e1 + e2 > e3 && e1 + e3 > e2 && e2 + e3 > e1);

	double p = (e1 + e2 + e3) / 2.0;
	return std::sqrt(p*(p-e1)*(p-e2)*(p-e3));
}

}// end of namespace


#endif