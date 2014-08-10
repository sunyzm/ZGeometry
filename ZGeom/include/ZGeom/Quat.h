#ifndef ZGEOM_QUAT_H
#define ZGEOM_QUAT_H
#include "Vector3D.h"

typedef double mat44[4][4];

class CQrot
{
public:
	CQrot(){m_w = 1; m_x = m_y = m_z = 0;};

	CQrot( double w, double x, double y, double z)
	{
		m_w = w;
		m_x = x;
		m_y = y;
		m_z = z;
	};

	CQrot( const CQrot & q )
	{
		m_w = q.m_w; 
		m_x =  q.m_x; 
		m_y =  q.m_y; 
		m_z =  q.m_z;
	}

	CQrot & operator=(const CQrot & q)
	{
		m_w = q.m_w;
		m_x  = q.m_x;
		m_y  = q.m_y;
		m_z  = q.m_z;
		return *this;
	}

	~CQrot(){};
  
	CQrot& operator^(double p);
	void convert( double * ) const;
	Vector3D operator*( const Vector3D & v );

	//multiplication
	friend CQrot    operator*(const CQrot &p, const CQrot & q );
	friend CQrot    operator*(const CQrot &p, double        s );

	friend CQrot    operator+(const CQrot &p, const CQrot & q );
	friend CQrot    operator-(const CQrot &p, const CQrot & q );
	friend double   operator^(const CQrot &p, const CQrot & q ); //dot product
	//power function
	friend CQrot    operator^(const CQrot &q, double p );

public:
	double m_w, m_x, m_y, m_z;
 
private:
	void normalize();
};

#endif