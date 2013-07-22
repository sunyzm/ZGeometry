#ifndef ZUTIL_COLOR_H
#define ZUTIL_COLOR_H

enum ColorID {COLOR_RED, COLOR_GREEN, COLOR_BLUE, COLOR_YELLOW, COLOR_PURPLE, COLOR_MAGENTA, COLOR_CYAN};

const float RGBColors[][3] = {	
	/*red*/{1.f, 0, 0},
	/*green*/{0, 1.f, 0},	
	/*blue*/{0, 0, 1.f},
	/*yellow*/{1.f, 1.f, 0},
	/*purple*/{159.f/255.f, 0, 197.f/255.f},
	/*magenta*/{1.f, 0, 1.f},
	/*cyan*/{0, 1.f, 1.f}
};

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
	float RedMap[256], GreenMap[256], BlueMap[256];
	void BuildLUT();
};



#endif