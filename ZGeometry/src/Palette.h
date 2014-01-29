#pragma once
#include <vector>
#include <QColor>
#include <ZGeom/Color.h>

class Palette
{
public:
	void generatePalette(int n) 
	{
		mPaletteColors.resize(n);
		
		for (int i = 0; i < n; ++i) {
			QColor rgbColor = QColor::fromHsvF((i+1.0)/(n+1.0), 0.75, 0.75).convertTo(QColor::Rgb);
			mPaletteColors[i][0] = rgbColor.redF();
			mPaletteColors[i][1] = rgbColor.greenF();
			mPaletteColors[i][2] = rgbColor.blueF();
			mPaletteColors[i][3] = rgbColor.alphaF();
		}
	}
	int totalColors() const { return mPaletteColors.size(); }
	const ZGeom::Colorf& getColor(int i) { return mPaletteColors[i]; }

private:
	std::vector<ZGeom::Colorf> mPaletteColors;
};