#ifndef ZGEOM_VEC3_H
#define ZGEOM_VEC3_H
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include "common.h"

namespace ZGeom {

template<typename T>
class Vec3
{
public:
    static Vec3<T> ZERO;

public:
	Vec3() : x(0.), y(0.), z(0.) {}
	Vec3(T x1, T y1, T z1) : x(x1), y(y1), z(z1) {}
	Vec3(const Vec3<T>& v2) : x(v2.x), y(v2.y), z(v2.z) {}	
    const Vec3<T>& operator = (const Vec3<T>& v2);

	Vec3<T> operator + (const Vec3<T>& v2) const;
	Vec3<T> operator - (const Vec3<T>& v2) const;
	Vec3<T> operator * (T coeff) const;
	Vec3<T> operator / (T coeff) const;
    Vec3<T> operator ^ (const Vec3<T>& v2) const;
	const Vec3<T>& operator += (const Vec3<T>& v2);
	const Vec3<T>& operator -= (const Vec3<T>& v2);
	const Vec3<T>& operator *= (T coeff);
	const Vec3<T>& operator /= (T coeff);
	T operator [] (uint i) const;
	T& operator [] (uint idx);

    bool equals(const Vec3<T>& v2, T eps = (T)1e-6) const;
    T length() const;
    T length2() const;
    T distanceTo(const Vec3<T>& v2) const { return (*this - v2).length(); }
    T dot(const Vec3<T>& v2) const { return x*v2.x + y*v2.y + z*v2.z;  }
    Vec3<T> cross(const Vec3<T> &v2) const;
    const Vec3<T>& normalize();

	template<typename U>
    operator Vec3<U>() { return Vec3<U>(this->x, this->y, this->z); }

public:
	T x, y, z;
};

typedef Vec3<float>  Vec3s;
typedef Vec3<double> Vec3d;

template<typename T>
Vec3<T> Vec3<T>::ZERO = Vec3<T>(0, 0, 0);

template<typename T>
const Vec3<T>& Vec3<T>::operator = (const Vec3<T>& v2) 
{ 
    x = v2.x; 
    y = v2.y; 
    z = v2.z; 
    return *this; 
}

template<typename T>
inline T Vec3<T>::operator [] ( uint i ) const
{
	switch (i) {
	case 0: return this->x;
	case 1: return this->y;
	case 2: return this->z;
	default: throw std::logic_error("Invalid Vec3 subscript");
	}
}

template<typename T>
inline T& Vec3<T>::operator [] ( uint idx )
{
	switch (idx) {
	case 0: return this->x;
	case 1: return this->y;
	case 2: return this->z;
	default: throw std::logic_error("Invalid Vec3 subscript");
	}
}

template<typename T>
bool Vec3<T>::equals(const Vec3<T>& v2, T eps = (T)1e-7) const 
{ 
    return fabs(x - v2.x) <= eps && fabs(y - v2.y) <= eps && fabs(z - v2.z) <= eps; 
}

template<typename T>
inline T Vec3<T>::length() const
{
	return std::sqrt(x*x + y*y + z*z);
}

template<typename T>
inline T Vec3<T>::length2() const
{
    return x*x + y*y + z*z;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::normalize() {
    if (x == 0 && y == 0 && z == 0) return *this;
    T len = length();
    x /= len; y /= len; z /= len;
    return *this;
}

template<typename T>
inline Vec3<T> operator -(const Vec3<T>& v) { return Vec3<T>(-v.x, -v.y, -v.z); }

template<typename T>
inline const Vec3<T>& Vec3<T>::operator += ( const Vec3<T>& v2 )
{
	this->x += v2.x;
	this->y += v2.y;
	this->z += v2.z;
	return *this;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator -= (const Vec3<T>& v2)
{
    *this += -v2;
    return *this;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator *= (T coeff)
{
    this->x *= coeff;
    this->y *= coeff;
    this->z *= coeff;
    return *this;
}

template<typename T>
inline const Vec3<T>& Vec3<T>::operator /= (T coeff)
{
    *this *= T(1.0) / coeff;
    return *this;
}

template<typename T>
inline Vec3<T> Vec3<T>::operator + (const Vec3<T>& v2) const
{
	return Vec3<T>(x + v2.x, y + v2.y, z + v2.z);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator - (const Vec3<T>& v2) const
{
	return Vec3<T>(x - v2.x, y - v2.y, z - v2.z);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator * (T coeff) const
{
	return Vec3<T>(x * coeff, y * coeff, z * coeff);
}

template<typename T>
inline Vec3<T> operator * (T coeff, const Vec3<T>& v1) { return v1 * coeff; }

template<typename T>
inline Vec3<T> Vec3<T>::operator / (T coeff) const
{
	return Vec3<T>(x / coeff, y / coeff, z / coeff);
}

template<typename T>
inline Vec3<T> Vec3<T>::operator ^ (const Vec3<T>& v2) const { return this->cross(v2); }

template<typename T>
inline Vec3<T> Vec3<T>::cross(const Vec3<T> &v2) const
{
    return Vec3<T>(this->y * v2.z - this->z * v2.y, 
                   this->z * v2.x - this->x * v2.z, 
                   this->x * v2.y - this->y * v2.x);
}

template<typename T> 
inline T dot(const Vec3<T>& v1, const Vec3<T>& v2) { return v1.dot(v2); }

template<typename T>
inline Vec3<T> cross( const Vec3<T>& v1, const Vec3<T>& v2 ) { return v1.cross(v2); }


template<typename T>
inline std::ostream& operator << (std::ostream& os, const Vec3<T>& vec)
{
    os << "(" << vec[0] << "," << vec[1] << "," << vec[2] << ")";
    return os;
}



/////////////////////////////////////////////////////////////
// Matrix3 : 3*3 matrix
/////////////////////////////////////////////////////////////

class Matrix3
{
private:
    double data[3][3];
public:
    Matrix3();
    Matrix3(const Matrix3& mat);
    Matrix3& operator =(const Matrix3& mat);

    friend Matrix3 vector3DMultiply(const Vec3d& v1, const Vec3d& v2);

    //operators overloading
    double& operator()(int i, int j);
    Matrix3& operator +=(const Matrix3& mat1);
    Matrix3& operator *=(double c);

    friend Matrix3 operator+(const Matrix3& mat1, const Matrix3& mat2);
    friend Matrix3 operator-(const Matrix3& mat1, const Matrix3& mat2);
    friend Vec3d operator*(const Vec3d& vec, const Matrix3& mat1);
    friend Vec3d operator*(const Matrix3& mat1, const Vec3d& vec);
};


} //end of namespace ZGeom

#endif

