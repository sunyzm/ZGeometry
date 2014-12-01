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
    ColorSignature() : mCurrentColorMap(CM_PARULA) {}
    ColorSignature(const std::vector<Colorf>& vCol);
    ColorSignature(const std::vector<double>& vVals, ColorMapType cmt = CM_PARULA, bool requireNormalize = true);
    
    const std::vector<Colorf>& getColorSig() const { return mColors; }
    std::vector<Colorf>& getColors() { return mColors; }
    std::vector<double>& getVals() { return mValues; }
    int size() const { return (int)mColors.size(); }

    void changeColorMap(ColorMapType cmt, bool requireNormalize = true);
    void changeSignatureMode(SignatureMode sigMode);

private:
    std::vector<Colorf> mColors;
    std::vector<double> mValues;
    ColorMapType mCurrentColorMap;
};

}

#endif