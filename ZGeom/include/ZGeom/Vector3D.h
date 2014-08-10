// Vector3D.h: classes of 3D vector and matrix3
//
//////////////////////////////////////////////////////////////////////
#ifndef ZMESH_GEOMETRY_H
#define ZMESH_GEOMETRY_H
#include <cmath>
#include <cstdio>
#include <string>

const double DOUBLE_EPS = 1e-16;
static inline bool EQUALZERO(double x) { return std::fabs((x)) < DOUBLE_EPS; }

/////////////////////////////////////////////////////////////
// Vector3D : 3D vector
/////////////////////////////////////////////////////////////

class Matrix3;

class Vector3D 
{
public :
	double x, y, z;

	// constructions
	Vector3D()									{ x = 0;	y = 0;		z = 0;	}
	Vector3D(double xx, double yy, double zz)	{ x = xx;	y = yy;		z = zz; }
	Vector3D(const Vector3D& v)					{ x = v.x;	y = v.y;	z = v.z;}

	// operator
	double	  length() const		{	return std::sqrt(x*x + y*y + z*z);	}
	double	  length2()	const	{	return x*x + y*y + z*z;	}
	Vector3D  normalize();
	double	  distantFrom( const Vector3D& v ) const { return std::sqrt( (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y) + (z-v.z)*(z-v.z) ); }
	bool	  equals(const Vector3D& v2) const { return abs(x-v2.x) < 1e-6 && abs(y-v2.y) < 1e-6 && abs(z-v2.z) < 1e-6; }
	Vector3D& operator =(const Vector3D& v);
	Vector3D& operator +=(const Vector3D& v);
	Vector3D& operator -=(const Vector3D& v);
	Vector3D& operator *=(double u);
	Vector3D& operator /=(double u);
	Vector3D& operator ^=(const Vector3D& v);
	double& operator[](int i);
	double  operator[](int i) const;
	operator std::string() const;

	friend Vector3D operator-(const Vector3D& v);
	friend Vector3D operator+(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator-(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator*(const double u, const Vector3D& rv);
	friend Vector3D operator*(const Vector3D& lv, const double u);
	friend Vector3D operator/(const Vector3D& lv, const double u);
	friend double   operator*(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator^(const Vector3D& lv, const Vector3D& rv);
	friend Vector3D operator*(const Vector3D& vec, const Matrix3& mat1);
	friend Vector3D operator*(const Matrix3& mat1, const Vector3D& vec);
};

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
	
	friend Matrix3 vector3DMultiply(const Vector3D& v1, const Vector3D& v2);
	
	//operators overloading
	double& operator()(int i, int j);
	Matrix3& operator +=(const Matrix3& mat1);
	Matrix3& operator *=(double c);

	friend Matrix3 operator+(const Matrix3& mat1, const Matrix3& mat2);
	friend Matrix3 operator-(const Matrix3& mat1, const Matrix3& mat2);
	friend Vector3D operator*(const Vector3D& vec, const Matrix3& mat1);
	friend Vector3D operator*(const Matrix3& mat1, const Vector3D& vec);
};

/////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////

double dotProduct3D(const Vector3D& v1, const Vector3D& v2);
bool crossProduct3D(const Vector3D&v1, const Vector3D&v2, Vector3D& v3);
Vector3D cross3D(const Vector3D& v1, const Vector3D& v2); 
Matrix3 vector3DMultiply(const Vector3D& v1, const Vector3D& v2);
Vector3D TriAreaNormal(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);
double TriArea(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);

#endif