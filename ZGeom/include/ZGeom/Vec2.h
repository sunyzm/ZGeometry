#ifndef ZMATH_VEC2_H
#define ZMATH_VEC2_H
#include <cmath>
#include <cassert>

namespace ZGeom {

template<typename T>
class Vec2
{
public:
	Vec2() : m_x(0.), m_y(0.) {}
	Vec2(T x1, T y1) : m_x(x1), m_y(y1) {}

	friend Vec2<T> operator + (const Vec2<T>& v1, const Vec2<T>& v2);
	friend Vec2<T> operator - (const Vec2<T>& v1, const Vec2<T>& v2);
	friend Vec2<T> operator * (const Vec2<T>& v, T lambda);
	friend Vec2<T> operator / (const Vec2<T>& v, T lambda);
	friend Vec2<T> operator - (const Vec2<T>& v);
	friend T dot(const Vec2<T>& v1, const Vec2<T>& v2);

	T x() const { return m_x; }
	T& x() { return m_x; }
	T y() const { return m_y; }
	T& y() { return m_y; }
	const Vec2<T>& operator += (const Vec2<T>& v2);
	const Vec2<T>& operator -= (const Vec2<T>& v2);
	const Vec2<T>& operator *= (T lambda);
	const Vec2<T>& operator /= (T lambda);
	T& operator [] (unsigned i) const;
	T length() const;

private:
	T m_x, m_y;
};

template<typename T>
T& Vec2<T>::operator[]( unsigned i ) const
{
	assert(i == 0 || i == 1);
	if (i == 0) return this->m_x;
	else if (i == 1) return this->m_y;
}

template<typename T>
T Vec2<T>::length() const
{
	return std::sqrt(m_x*m_x + m_y*m_y);
}

template<typename T>
const Vec2<T>& Vec2<T>::operator/=( T lambda )
{
	this->m_x /= lambda;
	this->m_y /= lambda;
	return *this;
}

template<typename T>
const Vec2<T>& Vec2<T>::operator*=( T lambda )
{
	this->m_x *= lambda;
	this->m_y *= lambda;
	return *this;
}

template<typename T>
const Vec2<T>& Vec2<T>::operator-=( const Vec2<T>& v2 )
{
	this->m_x -= v2.x;
	this->m_y -= v2.y;
	return *this;
}

template<typename T>
const Vec2<T>& Vec2<T>::operator+=( const Vec2<T>& v2 )
{
	this->m_x += v2.x;
	this->m_y += v2.y;
	return *this;
}

template<typename T>
Vec2<T> operator - (const Vec2<T>& v)
{
	return Vec2<T>(-v.x, -v.y);
}

template<typename T>
Vec2<T> operator + (const Vec2<T>& v1, const Vec2<T>& v2)
{
	return Vec2<T>(v1.x + v2.x, v1.y + v2.y);
}

template<typename T>
Vec2<T> operator - (const Vec2<T>& v1, const Vec2<T>& v2)
{
	return Vec2<T>(v1.x - v2.x, v1.y - v2.y);
}

template<typename T>
Vec2<T> operator * (const Vec2<T>& v, T lambda)
{
	return Vec2<T>(v1.x * lambda, v1.y * lambda);
}

template<typename T>
Vec2<T> operator / (const Vec2<T>& v, T lambda)
{
	return Vec2<T>(v1.x / lambda, v1.y / lambda);
}

template<typename T> inline
T dot(const Vec2<T>& v1, const Vec2<T>& v2)
{
	return v1.x * v2.x + v1.y * v2.y;
}

typedef Vec2<float> Vec2s;
typedef Vec2<double> Vec2d;

} //end of namespace ZGeom


#endif