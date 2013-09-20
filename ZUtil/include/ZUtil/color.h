#ifndef ZUTIL_COLOR_H
#define ZUTIL_COLOR_H

#include <cassert>

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
	RGBf() : r(0), g(0), b(0) {} 
	RGBf(float x, float y, float z) : r(x), g(y), b(z) {}
	RGBf(const RGBf& v)	: r(v.r), g(v.g), b(v.b) {}
	RGBf(float grayscale) {
		assert(0 <= grayscale && grayscale <= 1);
		r = grayscale;
		if ( grayscale < 0.5) g = grayscale * 2.0f;
		else g = 2.0f * (1.0f - grayscale);
		b = 1.0f - grayscale;
	}

	float rgb2gray() const	
	{
		return float(0.2989 * r + 0.5870 * g + 0.1140 * b); 
	}
	operator float () const { return rgb2gray(); }
	const RGBf& operator = (const RGBf& v) { r = v.r; g = v.g; b = v.b; return (*this); }

public:
	float r, g, b;
};


class FalseColorMap
{
public:
	FalseColorMap();
	float RedMap[257], GreenMap[257], BlueMap[257];
};

inline FalseColorMap::FalseColorMap()
{
	//initialize false color map
	for (int i = 0; i < 256; i++)
	{
		RedMap[i] = i / 255.0f;
		if (i < 128) GreenMap[i] = i / 128.0f;
		else GreenMap[i] = 1.0f - (i-128.0f) / 127.0f;
		BlueMap[i] = 1.0f - i / 255.0f;
	}

	RedMap[256] = RedMap[255];
	GreenMap[256] = GreenMap[255];
	BlueMap[256] = BlueMap[255];
}

#endif