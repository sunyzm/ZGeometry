#ifndef ZMATH_VEC2_H
#define ZMATH_VEC2_H

#include <cmath>

namespace ZGeom {
	template<typename T>
	class Vec2
	{
	public:
		Vec2() : x(0.), y(0.) {}
		Vec2(T x1, T y1) : x(x1), y(y1) {}

		friend Vec2<T> operator + (const Vec2<T>& v1, const Vec2<T>& v2);
		friend Vec2<T> operator - (const Vec2<T>& v1, const Vec2<T>& v2);
		friend Vec2<T> operator * (const Vec2<T>& v, T lambda);
		friend Vec2<T> operator / (const Vec2<T>& v, T lambda);
		friend Vec2<T> operator - (const Vec2<T>& v);
		friend T dot(const Vec2<T>& v1, const Vec2<T>& v2);

		const Vec2<T>& operator += (const Vec2<T>& v2);
		const Vec2<T>& operator -= (const Vec2<T>& v2);
		const Vec2<T>& operator *= (T lambda);
		const Vec2<T>& operator /= (T lambda);
		T& operator [] (unsigned i) const;
		T length() const;

		T x, y;
	};

	template<typename T>
	T& Vec2<T>::operator[]( unsigned i ) const
	{
		assert(i == 0 || i == 1);
		if (i == 0) return this->x;
		else if (i == 1) return this->y;
	}

	template<typename T>
	T Vec2<T>::length() const
	{
		return std::sqrt(x*x + y*y);
	}

	template<typename T>
	const Vec2<T>& Vec2<T>::operator/=( T lambda )
	{
		this->x /= lambda;
		this->y /= lambda;
		return *this;
	}

	template<typename T>
	const Vec2<T>& Vec2<T>::operator*=( T lambda )
	{
		this->x *= lambda;
		this->y *= lambda;
		return *this;
	}

	template<typename T>
	const Vec2<T>& Vec2<T>::operator-=( const Vec2<T>& v2 )
	{
		this->x -= v2.x;
		this->y -= v2.y;
		return *this;
	}

	template<typename T>
	const Vec2<T>& Vec2<T>::operator+=( const Vec2<T>& v2 )
	{
		this->x += v2.x;
		this->y += v2.y;
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
}


#endif