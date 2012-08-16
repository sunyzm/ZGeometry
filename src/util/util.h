#pragma once
#include <string>
#include <sstream>
#include <cmath>
#include "color.h"
#define PI 3.14159265358979323846

inline std::string Int2String(int i)
{
	std::ostringstream ostr;
	ostr << i << std::flush;
	return ostr.str();
}

inline std::string Double2String(double f)
{
	std::ostringstream ostr;
	ostr << f << std::flush;
	return ostr.str();
}

inline double sinc(double x)
{
	if (fabs(x)<1e-8) return 1.0;
	else return std::sin(PI*x) / (PI*x);
}