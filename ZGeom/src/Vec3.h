#ifndef ZGEOM_VEC3_H
#define ZGEOM_VEC3_H

#include <cmath>
#include <cassert>

namespace ZGeom {
template<typename T>
class Vec3
{
public:
	Vec3() : x(0.), y(0.), z(0.) {}
	Vec3(T x1, T y1, T z1) : x(x1), y(y1), z(z1) {}

	friend Vec3<T> operator + (const Vec3<T>& v1, const Vec3<T>& v2);
	friend Vec3<T> operator - (const Vec3<T>& v1, const Vec3<T>& v2);
	friend Vec3<T> operator * (const Vec3<T>& v, T lambda);
	friend Vec3<T> operator / (const Vec3<T>& v, T lambda);
	friend Vec3<T> operator - (const Vec3<T>& v);
	friend T dot(const Vec3<T>& v1, const Vec3<T>& v2);
	friend Vec3<T> cross(const Vec3<T>& v1, const Vec3<T>& v2);

	T x() const { return this->x; }
	T y() const { return this->y; }
	T z() const { return this->z; }

	const Vec3<T>& operator += (const Vec3<T>& v2);
	const Vec3<T>& operator -= (const Vec3<T>& v2);
	const Vec3<T>& operator *= (T lambda);
	const Vec3<T>& operator /= (T lambda);
	T& operator [] (unsigned i) const;
	T length() const;
		
	template<typename U>
	operator Vec3<U>()
	{
		return Vec3<U>(this->x, this->y, this->z);
	}

private:
	T x, y, z;
};

template<typename T>
T& Vec3<T>::operator[]( unsigned i ) const
{
	assert(i == 0 || i == 1 || i == 2);
	switch (i) {
	case 0: return this->x;
	case 1: return this->y;
	case 2: return this->z;
	}
}

template<typename T>
T Vec3<T>::length() const
{
	return std::sqrt(x*x + y*y + z*z);
}

template<typename T>
const Vec3<T>& Vec3<T>::operator/=( T lambda )
{
	this->x /= lambda;
	this->y /= lambda;
	this->z /= lambda;
	return *this;
}

template<typename T>
const Vec3<T>& Vec3<T>::operator*=( T lambda )
{
	this->x *= lambda;
	this->y *= lambda;
	this->z *= lambda;
	return *this;
}

template<typename T>
const Vec3<T>& Vec3<T>::operator-=( const Vec3<T>& v2 )
{
	this->x -= v2.x;
	this->y -= v2.y;
	this->z -= v2.z;
	return *this;
}

template<typename T>
const Vec3<T>& Vec3<T>::operator+=( const Vec3<T>& v2 )
{
	this->x += v2.x;
	this->y += v2.y;
	this->z += v2.z;
	return *this;
}

template<typename T>
Vec3<T> operator - (const Vec3<T>& v)
{
	return Vec3<T>(-v.x, -v.y, -v.z);
}

template<typename T>
Vec3<T> operator + (const Vec3<T>& v1, const Vec3<T>& v2)
{
	return Vec3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template<typename T>
Vec3<T> operator - (const Vec3<T>& v1, const Vec3<T>& v2)
{
	return Vec3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

template<typename T>
Vec3<T> operator * (const Vec3<T>& v, T lambda)
{
	return Vec3<T>(v1.x * lambda, v1.y * lambda, v1.z * lambda);
}

template<typename T>
Vec3<T> operator / (const Vec3<T>& v, T lambda)
{
	return Vec3<T>(v1.x / lambda, v1.y / lambda, v1.z / lambda);
}

template<typename T> inline
T dot(const Vec3<T>& v1, const Vec3<T>& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<typename T> inline
Vec3<T> cross( const Vec3<T>& v1, const Vec3<T>& v2 )
{
	return Vec3<T>(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}


typedef Vec3<float> Vec2s;
typedef Vec3<double> Vec2d;

} //end of namespace ZGeom

#endif

