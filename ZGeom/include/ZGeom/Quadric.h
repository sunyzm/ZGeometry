#ifndef ZGEOM_QUADRIC_H
#define ZGEOM_QUADRIC_H
#include "Vector3D.h"

class Quadric
{
public:
	Quadric();
	Quadric(double a, double b, double c, double d, double r);
	Quadric(const Quadric& q);
	Quadric& operator =(const Quadric& q);

	double		eval(const Vector3D& v);
	Matrix3		getTensor() const;
	Vector3D	getVector() const;
	double		getOffset() const {return d2;}
	double		getArea() const {return area;}
	Quadric&	operator +=(const Quadric& q1);
	Quadric&	operator -=(const Quadric& q1);
	Quadric&	operator *=(double coeff);
	
	friend Quadric operator +(const Quadric& q1, const Quadric& q2);
	friend Quadric operator -(const Quadric& q1, const Quadric& q2);
	friend Quadric operator *(const Quadric& q1, double coeff);

private:
	void initInternal();

	// Attributes
	double a2, ab, ac, ad;
	double b2, bc, bd;
	double c2, cd;
	double d2;
	double area;

	double* pData[11];
};

#endif