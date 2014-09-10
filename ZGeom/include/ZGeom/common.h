#ifndef ZGEOM_COMMON_H
#define ZGEOM_COMMON_H
#include <functional>

namespace ZGeom {

const double PI = 3.14159265358979323846;	

typedef unsigned int uint;

template<typename T>
inline T sqr(T x) { return x*x; }

template<typename T>
inline T cubic(T x) { return x*x*x; }

} 


#endif