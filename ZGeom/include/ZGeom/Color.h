#ifndef ZGEOM_COLOR_H
#define ZGEOM_COLOR_H
#include <vector>

namespace ZGeom {

const float ColorBlack[4]		= {0, 0, 0, 1.f};
const float ColorWhite[4]		= {1.f, 1.f, 1.f, 1.f};
const float ColorMagenta[4]		= {1.f, 0, 1.f, 1.f};
const float ColorGreen[4]		= {0, 1.f, 0, 1.f};
const float ColorRed[4]			= {1.f, 0, 0, 1.f};
const float ColorDarkRed[4]		= {1.28f/2.55f, 0, 0, 1.f};
const float ColorBlue[4]		= {0, 0, 1.f, 1.f};
const float ColorDarkBlue[4]	= {0, 0, 1.43f/2.55f};
const float ColorPurple[4]		= {1.59f/2.55f, 0, 1.97f/2.55f, 1.f};
const float ColorCyan[4]		= {0, 1.f, 1.f, 1.f};
const float ColorYellow[4]		= {1.f, 1.f, 0, 1.f};
const float ColorOrange[4]		= {1.f, 0.5f, 0, 1.f};
const float ColorAzure[4]		= {0, 0.5f, 1.f, 1.f};
const float ColorPaleVioletRed2[4] = { 238.f / 255.f, 121.f / 255.f, 159.f / 255.f, 1.f };
const float ColorPink[4]        = { 1.f, 192.f / 255.f, 203.f / 255.f };
const float ColorSlateGray[4]   = { 112.f / 255.f, 128.f / 255.f, 144.f / 255.f };
const float ColorLightGray[4]   = { 119.f / 255.f, 136.f / 255.f, 153.f / 255.f };

const float* const RGBPresetColors[] = { ColorRed, ColorGreen, ColorBlue, ColorYellow, ColorPurple, ColorMagenta, ColorCyan };	

const float ColorMesh1[4] = { 0.53f, 0.70f, 0.93f, 1.0f };
const float ColorMesh2[4] = { 0.99f, 0.73f, 0.62f, 1.0f };
const float ColorMesh3[4] = { 0.63f, 0.78f, 0.63f, 1.0f };
const float* const MeshPresetColors[] = { ColorLightGray, ColorMesh1, ColorMesh2, ColorMesh3 };

class Colorf
{
public:
	Colorf();
	Colorf(float r , float g , float b, float a = 1);
	Colorf(const float *c);

	float r() const { return mVal[0]; }
	float g() const { return mVal[1]; }
	float b() const { return mVal[2]; }
	float a() const { return mVal[3]; }
	float toGrayscale() const { return 0.2989f * r() + 0.5870f * g() + 0.1140f * b(); }
	void posNegColor(float val, const float* colorPos = ColorOrange, const float* colorNeg = ColorAzure);
	void setAs(const float *c);
    float& operator[] (int i) { return mVal[i]; }    
    float operator[] (int i) const { return mVal[i]; }

    static void interpolateColor(const float *color1, const float *color2, float coeff1, Colorf &color3)
    {
        float coeff2 = 1 - coeff1;
        for (int i = 0; i < 3; ++i) {
            color3[i] = color1[i] * coeff2 + color2[i] * coeff1;
        }
    }

private:
    std::vector<float> mVal;
};

}   // end of namespace

#endif
