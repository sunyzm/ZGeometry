// Vector3D.cpp: implementation of the Vector3D class.
//
//////////////////////////////////////////////////////////////////////
#include "Vector3D.h"
#include <cassert>
#include <cstdio>
#include <memory>
#include <queue>
#include <functional>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Vector3D::operator std::string() const
{
	char buffer[50];
	sprintf_s(buffer, 50, "(%f,%f,%f)", x, y, z);
	return std::string(buffer);
}

/////////////////////////////////////////////////////////////
// Vector3D : 3D vector
/////////////////////////////////////////////////////////////
Vector3D& Vector3D::operator=(const Vector3D& v)
{	
	x = v.x;	y = v.y;	z = v.z;
	return (*this);
}
Vector3D& Vector3D::operator+=(const Vector3D& v)
{	
	x += v.x;	y += v.y;	z += v.z;	
	return (*this);	
}
Vector3D& Vector3D::operator-=(const Vector3D& v)
{	
	x -= v.x;	y -= v.y;	z -= v.z;	
	return (*this);	
}
Vector3D& Vector3D::operator*=(double u)
{	
	x *= u;		y *= u;		z *= u;		
	return (*this);	
}
Vector3D& Vector3D::operator/=(double u)
{	
	if (!EQUALZERO(u))
	{x /= u;		y /= u;		z /= u;}
	return(*this);
}

Vector3D& Vector3D::operator^=(const Vector3D& v)
{	
	double xx = y*v.z - z*v.y;	
	double yy = z*v.x - x*v.z;	
	double zz = x*v.y - y*v.x;	
	x = xx; y = yy; z = zz; 
	return (*this);	
}

double& Vector3D::operator[]( int i )
{
	assert(0 <= i && i <= 2);
	if (i == 0) return x;
	else if (i == 1) return y;
	else return z;
}

double Vector3D::operator[]( int i ) const
{
	assert(0 <= i && i <= 2);
	if (i == 0) return x;
	else if (i == 1) return y;
	else return z;
}

Vector3D Vector3D::normalize()
{
	double len = length();
	if(!EQUALZERO(len)) 
	{
		x /= len; y /= len; z /= len;
	} 
	return *this;
}

Vector3D operator-(const Vector3D& v)
{
	Vector3D ret;
	ret.x = -v.x;
	ret.y = -v.y;
	ret.z = -v.z;
	return ret;
}

Vector3D operator+(const Vector3D& lv, const Vector3D& rv)
{
	Vector3D rel = lv;
	rel += rv;
	return rel;
}

Vector3D operator-(const Vector3D& lv, const Vector3D& rv)
{
	Vector3D rel = lv;
	rel -= rv;
	return rel;
}

Vector3D operator*(const double u, const Vector3D& rv)
{
	Vector3D rel = rv;
	rel *= u;
	return rel;
}

Vector3D operator*(const Vector3D& lv, const double u)
{
	Vector3D rel = lv;
	rel *= u;
	return rel;
}

Vector3D operator/(const Vector3D& lv, const double u)
{
	Vector3D rel = lv;
	rel /= u;
	return rel;
}

double operator*(const Vector3D& lv, const Vector3D& rv)
{
	return lv.x*rv.x + lv.y*rv.y + lv.z*rv.z;
}

Vector3D operator^(const Vector3D& lv, const Vector3D& rv)
{
	Vector3D rel = lv;
	rel ^= rv;
	return rel;
}

double dotProduct3D( const Vector3D& v1, const Vector3D& v2 )
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

bool crossProduct3D( const Vector3D&v1, const Vector3D&v2, Vector3D& v3 )
{
	v3.x = v1.y*v2.z - v1.z*v2.y;
	v3.y = v1.z*v2.x - v1.x*v2.z;
	v3.z = v1.x*v2.y - v1.y*v2.x;
	return true;
}

Vector3D cross3D(const Vector3D& v1, const Vector3D& v2) 
{
	Vector3D v3;
	v3.x = v1.y*v2.z - v1.z*v2.y;
	v3.y = v1.z*v2.x - v1.x*v2.z;
	v3.z = v1.x*v2.y - v1.y*v2.x;
	return v3;
}


/////////////////////////////////////////////////////////////
// Matrix3 : 3*3 matrix
/////////////////////////////////////////////////////////////
Matrix3::Matrix3()
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			data[i][j] = 0.0;
}

Matrix3::Matrix3( const Matrix3& mat )
{
	std::memcpy((void*)data, (void*)mat.data, sizeof(data));
}

Matrix3& Matrix3::operator=( const Matrix3& mat )
{
	std::memcpy((void*)data, (void*)mat.data, sizeof(data));
	return *this;
}

double& Matrix3::operator()( int i, int j )
{
	assert(0 <= i && i <= 2 && 0 <= j && j <= 2);
	return data[i][j];
}

Matrix3& Matrix3::operator*=( double c )
{
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			data[i][j] *= c;
	return *this;
}

Matrix3& Matrix3::operator+=( const Matrix3& mat1 )
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			data[i][j] += mat1.data[i][j];
		}
	}
	return *this;
}

Matrix3 operator+(const Matrix3& mat1, const Matrix3& mat2)
{
	Matrix3 mat;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			mat.data[i][j] = mat1.data[i][j] + mat2.data[i][j];
	return mat;
}

Matrix3 operator-(const Matrix3& mat1, const Matrix3& mat2)
{
	Matrix3 mat;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			mat.data[i][j] = mat1.data[i][j] - mat2.data[i][j];
	return mat;
}

Vector3D operator*( const Vector3D& vec, const Matrix3& mat1 )
{
	Vector3D ret;
	Matrix3 mat(mat1);

	ret.x = vec.x * mat(0,0) + vec.y * mat(1,0) + vec.z * mat(2,0);
	ret.y = vec.x * mat(0,1) + vec.y * mat(1,1) + vec.z * mat(2,1);
	ret.z = vec.x * mat(0,2) + vec.y * mat(1,2) + vec.z * mat(2,2); 

	return ret;
}

Vector3D operator*( const Matrix3& mat1, const Vector3D& vec )
{
	Vector3D ret;
	Matrix3 mat(mat1);

	ret.x = mat(0,0) * vec.x + mat(0,1) * vec.y + mat(0,2) * vec.z;
	ret.x = mat(1,0) * vec.x + mat(1,1) * vec.y + mat(1,2) * vec.z;
	ret.x = mat(2,0) * vec.x + mat(2,1) * vec.y + mat(2,2) * vec.z;
	return ret;
}

Matrix3 vector3DMultiply( const Vector3D& v1, const Vector3D& v2 )
{
	Matrix3 mat;
	
	mat.data[0][0] = v1.x * v2.x;
	mat.data[0][1] = v1.x * v2.y;
	mat.data[0][2] = v1.x * v2.z;
	mat.data[1][0] = v1.y * v2.x;
	mat.data[1][1] = v1.y * v2.y;
	mat.data[1][2] = v1.y * v2.z;
	mat.data[2][0] = v1.z * v2.x;
	mat.data[2][1] = v1.z * v2.y;
	mat.data[2][2] = v1.z * v2.z;
	
	return mat;
}

Vector3D TriAreaNormal( const Vector3D& v1, const Vector3D& v2, const Vector3D& v3 )
{
	Vector3D v12 = v2 - v1, v23 = v3 - v2;
	return (v12 ^ v23);
}

double TriArea( const Vector3D& v1, const Vector3D& v2, const Vector3D& v3 )
{
	return TriAreaNormal(v1, v2, v3).length() / 2.;
}
