#pragma once
#include <vector>

const double PI = 3.14159265358979323846;

inline double sinc(double x)
{
	if (fabs(x)<1e-10) return 1.0;
	else return std::sin(PI*x) / (PI*x);
}

double VectorDotProduct(const std::vector<double>& v1, const std::vector<double>& v2);

double VectorScalarProduct(const std::vector<double>& v1, const std::vector<double>& v2, const std::vector<double>& s);

void VectorPointwiseProduct(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v);