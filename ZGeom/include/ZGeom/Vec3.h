#ifndef ZGEOM_VEC3_H
#define ZGEOM_VEC3_H
#include <cmath>
#include <cassert>
#include <stdexcept>
#include "common.h"

namespace ZGeom {

template<typename T>
class Vec3
{
public:
	Vec3() : x(0.), y(0.), z(0.) {}
	Vec3(T x1, T y1, T z1) : x(x1), y(y1), z(z1) {}
	Vec3(const Vec3<T>& v2) : x(v2.x), y(v2.y), z(v2.z) {}
	const Vec3<T>& operator = (const Vec3<T>& v2) { x = v2.x; y = v2.y; z = v2.z; return *this; }

	Vec3<T> operator + (const Vec3<T>& v2) const;
	Vec3<T> operator - (const Vec3<T>& v2) const;
	Vec3<T> operator * (T lambda) const;
	Vec3<T> operator / (T lambda) const;
    Vec3<T> operator ^ (const Vec3<T>& v2) const;
	const Vec3<T>& operator += (const Vec3<T>& v2);
	const Vec3<T>& operator -= (const Vec3<T>& v2);
	const Vec3<T>& operator *= (T coeff);
	const Vec3<T>& operator /= (T coeff);
	T operator [] (uint i) const;
	T& operator [] (uint idx);
	T length() const;
    T dot(const Vec3<T>& v2) { return x*v2.x + y*v2.y + z*v2.z;  }

	friend Vec3<T> operator - (const Vec3<T>& v);
	friend T dot(const Vec3<T>& v1, const Vec3<T>& v2);
	friend Vec3<T> cross(const Vec3<T>& v1, const Vec3<T>& v2);

	template<typename U>
    operator Vec3<U>() { return Vec3<U>(this->x, this->y, this->z); }

public:
	T x, y, z;
};

template<typename T>
inline T Vec3<T>::operator[]( uint i ) const
{
	switch (i) {
	case 0: return this->x;
	case 1: return this->y;
	case 2: return this->z;
	default: throw std::logic_error("Invalid Vec3 subscript");
	}
}

template<typename T>
inline T& Vec3<T>::operator[]( uint idx )
{
	switch (idx) {
	case 0: return this->x;
	case 1: return this->y;
	case 2: return this->z;
	default: throw std::logic_error("Invalid Vec3 subscript");
	}
}

template<typename T>
inline T Vec3<T>::length() const
{
	return std::sqrt(x*x + y*y + z*z);
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator /= ( T coeff )
{
	this->x /= coeff;
	this->y /= coeff;
	this->z /= coeff;
	return *this;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator *= ( T coeff )
{
	this->x *= coeff;
	this->y *= coeff;
	this->z *= coeff;
	return *this;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator -= ( const Vec3<T>& v2 )
{
	this->x -= v2.x;
	this->y -= v2.y;
	this->z -= v2.z;
	return *this;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator += ( const Vec3<T>& v2 )
{
	this->x += v2.x;
	this->y += v2.y;
	this->z += v2.z;
	return *this;
}

template<typename T>
inline Vec3<T> operator -(const Vec3<T>& v)
{
	return Vec3<T>(-v.x, -v.y, -v.z);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator +(const Vec3<T>& v2) const
{
	return Vec3<T>(x + v2.x, y + v2.y, z + v2.z);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator -(const Vec3<T>& v2) const
{
	return Vec3<T>(x - v2.x, y - v2.y, z - v2.z);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator *(T coeff) const
{
	return Vec3<T>(x * coeff, y * coeff, z * coeff);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator /(T coeff) const
{
	return Vec3<T>(x / coeff, y / coeff, z / coeff);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator ^(const Vec3<T>& v2) const
{
    return Vec3<T>(this->y*v2.z - this->z*v2.y, this->z*v2.x - this->x*v2.z, this->x*v2.y - this->y*v2.x);
}

template<typename T> 
inline T dot(const Vec3<T>& v1, const Vec3<T>& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<typename T>
inline Vec3<T> cross( const Vec3<T>& v1, const Vec3<T>& v2 )
{
	return Vec3<T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

typedef Vec3<float>  Vec3s;
typedef Vec3<double> Vec3d;

} //end of namespace ZGeom

#endif

