#include <util/arithmetic.h>
#include <cassert>
#include <stdexcept>

double VectorDotProduct( const std::vector<double>& v1, const std::vector<double>& v2 )
{
	if(v1.size() != v2.size())
		throw std::logic_error("incompatible eigenfunctions!");

	double sum(0);
	int size = v1.size();
	for (int i = 0; i < size; ++i)
		sum += v1[i] * v2[i];

	return sum;
}

double VectorScalarProduct( const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& s )
{
	assert(v1.size() == v2.size() && v1.size() == s.size());

	double sum(0);
	int size = v1.size();
	for (int i = 0; i < size; ++i)
		sum += v1[i] * s[i] * v2[i];

	return sum;
}

void VectorPointwiseProduct( const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v )
{
	assert(v1.size() == v2.size());
	v.resize(v1.size());
	for (size_t i = 0; i < v1.size(); ++i)
	{
		v[i] = v1[i] * v2[i];
	}
}

void VectorPointwiseDivide( const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v )
{
	assert(v1.size() == v2.size());
	v.resize(v1.size());
	for (size_t i = 0; i < v1.size(); ++i)
	{
		v[i] = v1[i] / v2[i];
	}
}

void triangleCotan( double a, double b, double c, double &cotan_a, double &cotan_c )
{
	double cosa = (b*b+c*c-a*a)/(2.0*b*c);
	double cosc = (b*b+a*a-c*c)/(2.0*b*a);
	cotan_a = cosa / std::sqrt(1-cosa*cosa);
	cotan_c = cosc / std::sqrt(1-cosc*cosc);
}
