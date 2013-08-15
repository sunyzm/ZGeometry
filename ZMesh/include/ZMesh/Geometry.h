// Geometry.h: classes of 2D vector and 3D vector and matrix3
//
//////////////////////////////////////////////////////////////////////
#ifndef ZMESH_GEOMETRY_H
#define ZMESH_GEOMETRY_H

#include <cmath>
#include <cstdio>
#include <string>

const double DOUBLE_EPS = 1e-16;
static inline bool	EQUALZERO(double x) { return std::fabs((x)) < DOUBLE_EPS; }

/////////////////////////////////////////////////////////////
// Vector2D : 2D vector
/////////////////////////////////////////////////////////////
class Vector2D  
{
public :
	double x, y;

public:
	Vector2D(){	x = 0;	y = 0;}
	// constructions
	Vector2D(double xx, double yy)	{ x = xx; y = yy; }
	Vector2D(const Vector2D& v)	{ x = v.x; y = v.y; }

  
	// operator
	double	  length()		{ return std::sqrt(x*x + y*y); }
	double	  length2()		{ return x*x + y*y;	}
	double	  normalize()	{ double len = length(); if (!EQUALZERO(len)) {x/=len; y/=len;}	return len;	}
	Vector2D& operator=(const Vector2D& v);
	Vector2D& operator+=(const Vector2D& v);
	Vector2D& operator-=(const Vector2D& v);
	Vector2D& operator*=(double u);
	Vector2D& operator/=(double u);
//	Vector2D& operator^=(const Vector2D& v);

	bool	Intersect(Vector2D v1,Vector2D v2,Vector2D v3,Vector2D v4);
	bool	Intersect(Vector2D v1,Vector2D v2);

	friend Vector2D operator+(const Vector2D& lv, const Vector2D& rv);
	friend Vector2D operator-(const Vector2D& lv, const Vector2D& rv);
	friend Vector2D operator*(const double u, const Vector2D& rv);
	friend Vector2D operator*(const Vector2D& lv, const double u);
	friend Vector2D operator/(const Vector2D& lv, const double u);
	friend double   operator*(const Vector2D& lv, const Vector2D& rv);
//	friend Vector2D operator^(const Vector2D& lv, const Vector2D& rv);

	short	AtWhere(Vector2D v0,Vector2D v1);
	bool	AtRight(Vector2D v0,Vector2D v1);
	bool	AtLeft(Vector2D v0,Vector2D v1);
	bool	OnLine(Vector2D v0,Vector2D v1);
	double	GetArea(Vector2D v);
	

};

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
	double	  length()		{	return std::sqrt(x*x + y*y + z*z);	}
	double	  length2()		{	return x*x + y*y + z*z;	}
	Vector3D  normalize();
	double	  distantFrom( const Vector3D& v )  { return std::sqrt( (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y) + (z-v.z)*(z-v.z) ); }
	bool	  equals(const Vector3D& v2) const { return abs(x-v2.x) < 1e-6 && abs(y-v2.y) < 1e-6 && abs(z-v2.z) < 1e-6; }
	Vector3D& operator =(const Vector3D& v);
	Vector3D& operator +=(const Vector3D& v);
	Vector3D& operator -=(const Vector3D& v);
	Vector3D& operator *=(double u);
	Vector3D& operator /=(double u);
	Vector3D& operator ^=(const Vector3D& v);
	double& operator[](int i);
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
// VectorND : n-dim vector
/////////////////////////////////////////////////////////////
class VectorND
{
public:
	int m_size;
	double* m_vec;

public:
	VectorND() { m_size = 0; m_vec = NULL; }
	VectorND(int n) { m_size = n; m_vec = new double[m_size]; }
	virtual ~VectorND() { if(m_vec) delete []m_vec; }
	VectorND& operator= (const VectorND& vnd); 
public:
	void resize(int n);	// equivalent to reserve
	void reserve(int n);
	bool empty() { return !m_vec; }
	virtual void clear() { m_size = 0; if(!empty()) delete m_vec; m_vec=NULL; }
	double length();
	double length2();
	double normalize();
	double normalizedDifference(const VectorND& v) const;

	void   calSubtract(VectorND& v1, VectorND& v2);
	double calDistance2(const VectorND& v) const;
	double calDistance(const VectorND& v) const;
	double calDistance(const VectorND& v, int method) const;
	double norm_2() const;
	double norm_inf() const;
};

VectorND operator-(const VectorND& v1, const VectorND& v2);

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
Matrix3 vector3DMultiply(const Vector3D& v1, const Vector3D& v2);

Vector3D TriAreaNormal(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);
double TriArea(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);

#endif