#include "Color.h"
#include <cassert>
#include <algorithm>

namespace ZGeom {

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
		0 <= c[2] && c[2] <= 1 && 0 <= c[3] && c[3] <= 1 );
	std::copy_n(c, 4, mVal.begin());
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