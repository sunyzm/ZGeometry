#include "Quadric.h"

using namespace std;

double Quadric::eval( const Vector3D& v )
{
	double x = v.x, y = v.y, z = v.z;
	double ret = x*x*a2	+ 2*x*y*ab	+ 2*x*z*ac	+ 2*x*ad
						+ y*y*b2	+ 2*y*z*bc	+ 2*y*bd
									+ z*z*c2	+ 2*z*cd
												+ d2;
	
	return ret;
}

Quadric& Quadric::operator *= ( double coeff )
{
	for (int i = 0; i < 10; ++i)		//area doesn't scale
		*pData[i] *= coeff;
	return *this;
}

Quadric::Quadric()
{
	initInternal();
	for (int i = 0; i < 11; ++i)
		*pData[i] = 0;
}

Quadric::Quadric( const Quadric& q )
{
	initInternal();
	for (int i = 0; i < 11; ++i)
		*pData[i] = *q.pData[i];
}

Quadric::Quadric( double a, double b, double c, double d, double r )
{
	initInternal();
	a2 = a*a;	ab = a*b;	ac = a*c;	ad = a*d;
	b2 = b*b;	bc = b*c;	bd = b*d;
	c2 = c*c;	cd = c*d;
	d2 = d*d;
	area = r;
}

Quadric& Quadric::operator += ( const Quadric& q1 )
{
	for (int i = 0; i < 11; ++i)
	{
		*pData[i] += *q1.pData[i];
	}
	return *this;
}

Quadric& Quadric::operator -= ( const Quadric& q1 )
{
	for (int i = 0; i < 11; ++i)
	{
		*pData[i] -= *q1.pData[i];
	}
	return *this;
}

void Quadric::initInternal()
{
	pData[0] = &a2; pData[1] = &ab; pData[2] = &ac; pData[3] = &ad;
	pData[4] = &b2; pData[5] = &bc; pData[6] = &bd;
	pData[7] = &c2; pData[8] = &cd;
	pData[9] = &d2;
	pData[10] = &area;
}

Matrix3 Quadric::getTensor() const
{
	Matrix3 tensor;
	tensor(0,0) = a2;
	tensor(0,1) = ab;
	tensor(0,2) = ac;
	tensor(1,0) = ab;
	tensor(1,1) = b2;
	tensor(1,2) = bc;
	tensor(2,0) = ac;
	tensor(2,1) = bc;
	tensor(2,2) = c2;
	return tensor;
}

Vector3D Quadric::getVector() const
{
	Vector3D vec;
	vec[0] = ad;
	vec[1] = bd;
	vec[2] = cd;
	return vec;
}

Quadric& Quadric::operator = ( const Quadric& q )
{
	for (int i = 0; i < 11; ++i)
		*pData[i] = *q.pData[i];
	return *this;
}


Quadric operator + (const Quadric& q1, const Quadric& q2)
{
	Quadric q(q1);
	q += q2;
	return q;
}

Quadric operator - (const Quadric& q1, const Quadric& q2)
{
	Quadric q(q1);
	q -= q2;
	return q;
}

Quadric operator * (const Quadric& q1, double coeff)
{
	Quadric q;
	q *= coeff;
	return q;
}