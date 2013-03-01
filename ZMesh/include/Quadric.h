#pragma once
#include "Geometry.h"

class Quadric
{
	// Attributes
private:
	double a2, ab, ac, ad;
	double b2, bc, bd;
	double c2, cd;
	double d2;
	double area;

	double* pData[11];
	
	// Methods
public:
	//constructors
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
};



