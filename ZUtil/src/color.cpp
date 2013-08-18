#include "color.h"
#include <cassert>

void FalseColorMap::BuildLUT()
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

FalseColorMap::FalseColorMap()
{
	BuildLUT();
}

RGBf::RGBf( float grayscale )
{
	assert(0 <= grayscale && grayscale <= 1);

	r = grayscale;
	if ( grayscale < 0.5) g = grayscale * 2.0f;
	else g = 2.0f * (1.0f - grayscale);
	b = 1.0f - grayscale;
}
