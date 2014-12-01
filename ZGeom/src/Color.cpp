#include "Color.h"
#include <cassert>
#include <algorithm>

namespace ZGeom {
	
inline void interpolateColor(const float *color1, const float *color2, float coeff1, Colorf &color3)
{
	float coeff2 = 1 - coeff1;
	for (int i = 0; i < 3; ++i) {
		color3[i] = color1[i] * coeff2 + color2[i] * coeff1;
	}		
}

Colorf::Colorf() : mVal(4, 0)
{
    mVal[3] = 1.f;
}

Colorf::Colorf(float r, float g, float b, float a) : mVal(4, 0)
{
	assert( 0 <= r && r <= 1 && 0 <= g && g <= 1 &&
			0 <= b && b <= 1 &&	0 <= a && a <= 1 );
	mVal[0] = r;
	mVal[1] = g;
	mVal[2] = b;
	mVal[3] = a;
}

Colorf::Colorf(const float *c) : mVal(4, 0)
{
    setAs(c);
}


void Colorf::setAs( const float *c )
{
	assert( 0 <= c[0] && c[0] <= 1 && 0 <= c[1] && c[1] <= 1 &&
		0 <= c[2] && c[2] <= 1 &&	0 <= c[3] && c[3] <= 1 );
	std::copy_n(c, 4, mVal.begin());
}

void Colorf::falseColor( float gray, float alpha, ColorMapType cmt)
{
	//assert(0 <= gray && gray <= 1.f);
	if (gray < 0) gray = 0;
	else if (gray > 1) gray = 1.f;
    int idx = int(gray * 255.999999);

	if (cmt == CM_JET) {
        mVal[0] = (float)ColorMap::jet[idx * 3];
        mVal[1] = (float)ColorMap::jet[idx * 3 + 1];
        mVal[2] = (float)ColorMap::jet[idx * 3 + 2];
	} 
    else if (cmt == CM_PARULA) {
        mVal[0] = (float)ColorMap::parula[idx * 3];
        mVal[1] = (float)ColorMap::parula[idx * 3 + 1];
        mVal[2] = (float)ColorMap::parula[idx * 3 + 2];
    }
	else if (cmt == CM_COOL) {
		interpolateColor(ColorCyan, ColorMagenta, gray, *this);			
	} 
	else if (cmt == CM_HOT) {
		if (gray < 0.375)
			interpolateColor(ColorBlack, ColorRed, gray / 0.375f, *this);
		else if(gray < 0.75)
			interpolateColor(ColorRed, ColorYellow, (gray - 0.375f)/0.375f, *this);
		else 
			interpolateColor(ColorYellow, ColorWhite, (gray - 0.75f)/0.25f, *this);
	}

	mVal[3] = alpha;
}

void Colorf::posNegColor( float val, const float* colorPos /*= ColorOrange*/, const float* colorNeg /*= ColorBlue*/ )
{
	assert(-1.f <= val && val <= 1.f);
	if (val >= 0) {
		for (int i = 0; i < 3; ++i) mVal[i] = val * colorPos[i] + (1.f-val) * ColorWhite[i];
	}
	else {
		for (int i = 0; i < 3; ++i) mVal[i] = -val * colorNeg[i] + (1.f+val) * ColorWhite[i];
	}
	mVal[3] = 1.f;
}

}   // end of namespace