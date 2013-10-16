#include "Color.h"
#include <cassert>
#include <algorithm>

namespace ZGeom
{
	const float* RGBColors[] = { ColorRed, ColorGreen, ColorBlue, ColorYellow, ColorPurple, ColorMagenta, ColorCyan };	

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

	void Colorf::falseColor( float gray, float alpha )
	{
		assert(0 <= gray && gray <= 1.f);
		mVal[0] = gray;
		mVal[1] = (gray < 0.5f) ? gray * 2.f : ((1.f - gray) * 2.f);
		mVal[2] = 1.f - gray;
		mVal[3] = alpha;
	}

	void Colorf::setAs( const float *c )
	{
		assert( 0 <= c[0] && c[0] <= 1 && 0 <= c[1] && c[1] <= 1 &&
				0 <= c[2] && c[2] <= 1 &&	0 <= c[3] && c[3] <= 1 );
		std::copy_n(c, 4, mVal);
	}

}