#pragma once
#include <vector>
#include <QColor>
#include <ZGeom/Color.h>

const float PresetMeshColor1[] = {0.53, 0.70, 0.93, 1.0};

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

		if (n == 1) mPaletteColors[0].setAs(PresetMeshColor1);
	}
	int totalColors() const { return mPaletteColors.size(); }
	const ZGeom::Colorf& getColor(int i) const { return mPaletteColors[i]; }

private:
	std::vector<ZGeom::Colorf> mPaletteColors;
};