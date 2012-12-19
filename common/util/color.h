#ifndef _COLOR_H_
#define _COLOR_H_

class RGBf 
{
public :
	float r, g, b;
	// constructions
	RGBf()	{ r = 0; g = 0; b = 0; } 
	RGBf(float x, float y, float z)	{r = x; g = y; b = z; }
	RGBf(const RGBf& v)	{r = v.r; g = v.g; b = v.b; }

	// operator
	float rgb2gray()	
	{
		return float(0.2989 * r + 0.5870 * g + 0.1140 * b); 
	}

	RGBf& operator=(const RGBf& v) {r = v.r; g = v.g; b = v.b; return (*this);}
};

class FalseColorMap
{
public:
	FalseColorMap();
	static double RedMap[256], GreenMap[256], BlueMap[256];
	static void BuildLUT();
};


#endif