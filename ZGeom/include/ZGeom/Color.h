#ifndef ZGEOM_COLOR_H
#define ZGEOM_COLOR_H
#include "Vec4.h"
#include <assert.h>

namespace ZGeom
{

const float ColorBlack[]   = {0, 0, 0, 1.f};
const float ColorWhite[]   = {1.f, 1.f, 1.f, 1.f};
const float ColorMagenta[] = {1.f, 0, 1.f, 1.f};
const float ColorGreen[]   = {0, 1.f, 0, 1.f};
const float ColorRed[]     = {1.f, 0, 0, 1.f};
const float ColorBlue[]    = {0, 0, 1.f, 1.f};
const float ColorPurple[]  = {159.f/255.f, 0, 197.f/255.f, 1.f};
const float ColorCyan[]    = {0, 1.f, 1.f, 1.f};
const float ColorYellow[]  = {1.f, 1.f, 0, 1.f};
const float ColorOrange[]  = {1.f, 0.5f, 0, 1.f};
const float ColorAzure[]   = {0, 0.5f, 1.f, 1.f};

const float* RGBColors[];

class Colorf : public Vec4s
{
public:
	Colorf();
	Colorf(float r , float g , float b, float a = 1);
	Colorf(const float *c) { setAs(c); }

	float r() const { return mVal[0]; }
	float g() const { return mVal[1]; }
	float b() const { return mVal[2]; }
	float a() const { return mVal[3]; }
	float toGrayscale() const { return 0.2989f * r() + 0.5870f * g() + 0.1140f * b(); }
	void falseColor(float gray, float alpha = 1.f);
	void posNegColor(float val, const float* colorPos = ColorOrange, const float* colorNeg = ColorAzure);

	void setAs(const float *c);
	operator const float* () const { return mVal; }
};

class FalseColorMap
{
public:
	FalseColorMap();
	float RedMap[257], GreenMap[257], BlueMap[257];

	static float red(float gray) { return gray; }
	static float green(float gray) { return (gray < 0.5f) ? gray * 2.f : (2.f * (1.f - gray)); }
	static float blue(float gray) { return 1.f - gray; }
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

} //end of namespace

#endif
