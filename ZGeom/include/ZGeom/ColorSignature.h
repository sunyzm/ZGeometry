#ifndef ZGEOM_COLOR_SIGNATURE_H
#define ZGEOM_COLOR_SIGNATURE_H
#include <vector>
#include "Color.h"
#include "ColorMap.h"

namespace ZGeom {

enum SignatureMode {
    SM_Normalized = 0, SM_LogNormalized, SM_MarkNegNormalized,
    SM_AbsNormalized, SM_BandCurved, SM_PosNegPlot,
    SM_CountSigModes
};

class ColorSignature
{
public:
    ColorSignature() : mCurrentColorMap(CM_JET) {}
    ColorSignature(const std::vector<Colorf>& vCol);
    ColorSignature(const std::vector<double>& vVals, ColorMapType cmt = CM_JET, bool requireNormalize = true);
    
    const std::vector<Colorf>& getColorSig() const { return mColors; }
    std::vector<Colorf>& getColors() { return mColors; }
    std::vector<double>& getVals() { return mValues; }
    int size() const { return (int)mColors.size(); }

    void changeColorMap(ColorMapType cmt);
    void curve(double sMin, double sMax);
    void changeSignatureMode(SignatureMode sigMode);
    Colorf& operator [] (int i) { return mColors[i]; }

private:
    std::vector<double> mOriginalValues;
    std::vector<double> mValues;
    std::vector<Colorf> mColors;
    ColorMapType mCurrentColorMap;

    std::vector<double> vec_normalize(const std::vector<double>& vec);
};

}

#endif