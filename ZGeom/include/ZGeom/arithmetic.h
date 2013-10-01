#ifndef ZGEOM_ARITHMETIC_H
#define ZGEOM_ARITHMETIC_H

#include <vector>
#include <cmath>
#include <cassert>
#include "common.h"
#include "Vec3.h"

namespace ZGeom
{
	inline double sinc(double x)
	{
		if (fabs(x) < 1e-10) return 1.0;
		else return std::sin(PI*x) / (PI*x);
	}

	inline double VectorDotProduct( const std::vector<double>& v1, const std::vector<double>& v2 )
	{
		if(v1.size() != v2.size())
			throw std::logic_error("incompatible eigenfunctions!");

		double sum(0);
		size_t size = v1.size();
		for (size_t i = 0; i < size; ++i)
			sum += v1[i] * v2[i];

		return sum;
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

	/// judge whether the given three lengths can form a triangle
	inline bool validAsTriangle(double e1, double e2, double e3) {
		return e1 > 0 && e2 > 0 && e3 > 0 && 
			e1 + e2 > e3 && e1 + e3 > e2 && e2 + e3 > e1;
	}

	void triangleCotan( double a, double b, double c, double &cotan_a, double &cotan_c );

	/// compute cosine of the angle opposite to e3
	double cosTriSides(double e1, double e2, double e3);

	/// compute the area of the triangle of given side lengths
	double triArea(double e1, double e2, double e3);

	template<typename T> 
	inline T cotVec3(const Vec3<T>& v1, const Vec3<T>& v2)
	{
		return ZGeom::dot(v1, v2) / ZGeom::cross(v1, v2).length(); 
	}

	template<typename T> 
	inline T cosVec3(const Vec3<T>& v1, const Vec3<T>& v2)
	{
		return ZGeom::dot(v1, v2) / (v1.length() * v2.length());
	}
}



#endif