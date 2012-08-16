#ifndef _COLOR_H_
#define _COLOR_H_

class RGBf {
public :
	float r, g, b;
	// constructions
	RGBf()	{r=0; g=0; b=0;}
	RGBf(float x, float y, float z)	{r=x; g=y; b=z;}
	RGBf(const RGBf& v)	{r=v.r; g=v.g; b=v.b;}

	// operator
	float rgb2gray()	{
		return float(0.2989 * r + 0.5870 * g + 0.1140 * b); 
		//return (r+g+b)/3.0f;	
	}

	RGBf& operator=(const RGBf& v) {r=v.r; g=v.g; b=v.b; return (*this);}
};

#endif