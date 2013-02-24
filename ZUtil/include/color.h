#pragma once

class RGBf 
{
public :
	float r, g, b;
	// constructions
	RGBf()	{ r = 0; g = 0; b = 0; } 
	RGBf(float x, float y, float z)	 {r = x; g = y; b = z; }
	RGBf(const RGBf& v)	{ r = v.r; g = v.g; b = v.b; }
	RGBf(float grayscale);

	// operator
	float rgb2gray() const	
	{
		return float(0.2989 * r + 0.5870 * g + 0.1140 * b); 
	}

	operator float () const { return rgb2gray(); }

	RGBf& operator=(const RGBf& v) {r = v.r; g = v.g; b = v.b; return (*this);}
};

class FalseColorMap
{
public:
	FalseColorMap();
	double RedMap[256], GreenMap[256], BlueMap[256];
	void BuildLUT();
};
