#include "arithmetic.h"
#include <stdexcept>

double ZGeom::triArea( double e1, double e2, double e3 )
{
	if(!validAsTriangle(e1, e2, e3)) throw std::runtime_error("Invalid triangle side lengths");

	double p = (e1 + e2 + e3) / 2.0;
	return std::sqrt(p*(p-e1)*(p-e2)*(p-e3));
}

double ZGeom::cosTriSides( double e1, double e2, double e3 )
{
	if(!validAsTriangle(e1, e2, e3)) throw std::runtime_error("Invalid triangle side lengths");

	return (e1*e1 + e2*e2 - e3*e3) / (2.0*e1*e2);
}

void ZGeom::triangleCotan( double a, double b, double c, double &cotan_a, double &cotan_c )
{
	if(!validAsTriangle(a, b, c)) throw std::runtime_error("Invalid triangle side lengths");

	double cosa = (b*b+c*c-a*a)/(2.0*b*c);
	double cosc = (b*b+a*a-c*c)/(2.0*b*a);
	cotan_a = cosa / std::sqrt(1-cosa*cosa);
	cotan_c = cosc / std::sqrt(1-cosc*cosc);
}
