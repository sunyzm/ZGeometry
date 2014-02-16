#include "Color.h"
#include <cassert>
#include <algorithm>

namespace ZGeom
{
	const float* RGBColors[] = { ColorRed, ColorGreen, ColorBlue, ColorYellow, ColorPurple, ColorMagenta, ColorCyan };	
	
	void interpolateColor(const float *color1, const float *color2, float coeff1, float *color3)
	{
		float coeff2 = 1 - coeff1;
		for (int i = 0; i < 3; ++i) {
			color3[i] = color1[i] * coeff2 + color2[i] * coeff1;
		}		
	}

	Colorf::Colorf()
	{
		mVal[0] = mVal[1] = mVal[2] = 0;
		mVal[1] = 1.f;
	}

	Colorf::Colorf( float r, float g, float b, float a )
	{
		assert( 0 <= r && r <= 1 && 0 <= g && g <= 1 &&
				0 <= b && b <= 1 &&	0 <= a && a <= 1 );
		mVal[0] = r;
		mVal[1] = g;
		mVal[2] = b;
		mVal[3] = a;
	}

	void Colorf::setAs( const float *c )
	{
		assert( 0 <= c[0] && c[0] <= 1 && 0 <= c[1] && c[1] <= 1 &&
			0 <= c[2] && c[2] <= 1 &&	0 <= c[3] && c[3] <= 1 );
		std::copy_n(c, 4, mVal);
	}

	void Colorf::falseColor( float gray, float alpha )
	{
		assert(0 <= gray && gray <= 1.f);
		/*
		mVal[0] = gray;
		mVal[1] = (gray < 0.5f) ? gray * 2.f : ((1.f - gray) * 2.f);
		mVal[2] = 1.f - gray;
		*/
		if (gray < 0.25) {
			interpolateColor(ColorBlue, ColorCyan, gray / 0.25f, mVal);
		} else if (gray < 0.5) {
			interpolateColor(ColorCyan, ColorGreen, (gray-0.25f)/0.25f, mVal);
		} else if (gray < 0.75) {
			interpolateColor(ColorGreen, ColorYellow, (gray-0.5f)/0.25f, mVal);
		} else {
			interpolateColor(ColorYellow, ColorRed, (gray - 0.75f)/0.25f, mVal);
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

}