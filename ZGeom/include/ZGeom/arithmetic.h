#ifndef ZGEOM_ARITHMETIC_H
#define ZGEOM_ARITHMETIC_H
#include <vector>
#include <cmath>
#include <cassert>
#include "common.h"
#include "Vec3.h"

namespace ZGeom {

template<typename T>
inline T sqr(T x) { return x*x; }

template<typename T>
inline T cubic(T x) { return x*x*x; }

std::vector<double> linspace(double lo, double hi, int N);

template<typename T>
inline T mean(const std::vector<T>& vec)
{
	double sum(0);
	for (T a: vec) sum += a;
	return T(sum / (double)vec.size());
}

inline double sinc(double x)
{
	if (fabs(x) < 1e-10) return 1.0;
	else return std::sin(PI*x) / (PI*x);
}

inline double VectorScalarProduct( const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& s )
{
	assert(v1.size() == v2.size() && v1.size() == s.size());

	double sum(0);
	size_t size = v1.size();
	for (size_t i = 0; i < size; ++i)
		sum += v1[i] * s[i] * v2[i];

	return sum;
}

inline void VectorPointwiseProduct( const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v )
{
	assert(v1.size() == v2.size());
	v.resize(v1.size());
	for (size_t i = 0; i < v1.size(); ++i)
	{
		v[i] = v1[i] * v2[i];
	}
}

inline void VectorPointwiseDivide( const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v )
{
	assert(v1.size() == v2.size());
	v.resize(v1.size());
	for (size_t i = 0; i < v1.size(); ++i)
	{
		v[i] = v1[i] / v2[i];
	}
}

//// judge whether the given three lengths can form a triangle
inline bool validAsTriangle(double e1, double e2, double e3) {
	return e1 > 0 && e2 > 0 && e3 > 0 && 
		e1 + e2 > e3 && e1 + e3 > e2 && e2 + e3 > e1;
}

void triangleCot( double a, double b, double c, double &cotan_a, double &cotan_c );

//// compute cosine of the angle opposite to e3
double cosTriSides(double e1, double e2, double e3);

//// compute the area of the triangle of given side lengths
double triArea(double e1, double e2, double e3);

inline double triArea(Vec3d v1, Vec3d v2, Vec3d v3)
{
    return cross<double>(v2 - v1, v3 - v2).length() / 2.0;
}

inline Vec3d triNormal(Vec3d v1, Vec3d v2, Vec3d v3) { 
    return ((v3 - v1) ^ (v3 - v2)).normalize(); 
}

template<typename T> 
inline T cotVec3(const Vec3<T>& v1, const Vec3<T>& v2)
{
	return v1.dot(v2) / ZGeom::cross(v1, v2).length(); 
}

template<typename T> 
inline T cosVec3(const Vec3<T>& v1, const Vec3<T>& v2)
{
	return v1.dot(v2) / (v1.length() * v2.length());
}

inline double vecAngle(Vec3d v1, Vec3d v2) { return acos(cosVec3(v1, v2)); }

double calMixedTriArea(double a, double b, double c);
double calMixedTriArea(double a, double b, double c, double& cotan_a, double& cotan_c);
double calHalfMixedTriArea(double a, double b, double c, double& cotan_a);
void quadricForm(int dim1, int dim2, double* mat1, double* diag, double *matResult);
std::pair<Vec3d, double> circumcenter(Vec3d p1, Vec3d p2, Vec3d p3);

//////////////////////////////////////////////////////////////////////////
//// fast numerical computation with AMP
void quadricFormAMP(int dim1, int dim2, double* mat1, double* diag, double *matResult);
void matVecMulAMP(int dim1, int dim2, double *mat, double *vec, double *vResult);

} // end of namespace

#endif
