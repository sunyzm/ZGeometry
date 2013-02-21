#include <util/color.h>
#include <cassert>

void FalseColorMap::BuildLUT()
{
	//initialize false color map
	for (int i = 0; i < 256; i++)
	{
		RedMap[i] = i / 255.0;
		if (i < 128) GreenMap[i] = i / 128.0;
		else GreenMap[i] = 1.0 - (i-128.0) / 127.0;
		BlueMap[i] = 1.0 - i / 255.0;
	}
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
