#ifndef ZGEOM_COLOR_SIGNATURE_H
#define ZGEOM_COLOR_SIGNATURE_H
#include <vector>
#include "Color.h"
#include "ColorMap.h"

namespace ZGeom {

class ColorSignature
{
public:
    ColorSignature() : mCurrentColorMap(CM_JET) {}
    ColorSignature(const std::vector<Colorf>& vCol);
    ColorSignature(const std::vector<double>& vVals, ColorMapType cmt = CM_JET, bool requireNormalize = true);
    
    const std::vector<Colorf>& getColorSig() const { return mColors; }
    std::vector<Colorf>& getColors() { return mColors; }
    const std::vector<double>& getOriginalVals() const { return mOriginalValues; }
    std::vector<double>& getVals() { return mNormalizedValues; }
    int size() const { return (int)mColors.size(); }
    double getCurveMax() const { return mCurveMax; }
    double getCurveMin() const { return mCurveMin; }

    void changeColorMap(ColorMapType cmt);
    void curve(double sMin, double sMax);
    Colorf& operator [] (int i) { return mColors[i]; }
    bool hasOriginalValues() const { return !mOriginalValues.empty(); }
    void setColorBuffer(float* dst) const;

private:
    std::vector<double> mOriginalValues;
    double mCurveMin, mCurveMax;
    std::vector<double> mNormalizedValues;
    std::vector<Colorf> mColors;
    ColorMapType mCurrentColorMap;

    std::vector<double> vec_normalize(const std::vector<double>& vec);
};

}

#endif