#include "color.h"

double FalseColorMap::RedMap[256], FalseColorMap::GreenMap[256], FalseColorMap::BlueMap[256];

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
