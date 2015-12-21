#include "ColorSignature.h"
#include <algorithm>

namespace ZGeom {

ColorSignature::ColorSignature(const std::vector<Colorf>& vCol)
{
    mColors = vCol;
    mOriginalValues.resize(mColors.size());
    for (size_t i = 0; i < mColors.size(); ++i)
        mOriginalValues[i] = mColors[i].toGrayscale();
    
    mNormalizedValues = mOriginalValues;
    mCurrentColorMap = CM_JET;
}

ColorSignature::ColorSignature(const std::vector<double>& vVals, ColorMapType cmt /*= CM_PARULA*/, bool requireNormalize /*= true*/)
{
    mNormalizedValues = mOriginalValues = vVals;    
    if (requireNormalize) {
        auto iResult = std::minmax_element(mOriginalValues.begin(), mOriginalValues.end());
        mCurveMin = *(iResult.first);
        mCurveMax = *(iResult.second);        
    } else {
        mCurveMin = 0.;
        mCurveMax = 1.;
    }
    mCurrentColorMap = cmt;
    curve(mCurveMin, mCurveMax);
}

void ColorSignature::changeColorMap(ColorMapType cmt)
{
    mCurrentColorMap = cmt;
    if (mNormalizedValues.empty()) return;
    mColors.resize(mNormalizedValues.size());
    for (size_t i = 0; i < mNormalizedValues.size(); ++i) 
        mColors[i] = ZGeom::ColorMap::falseColor((float)mNormalizedValues[i], 1.f, cmt);
}

std::vector<double> ColorSignature::vec_normalize(const std::vector<double>& vec)
{
    std::vector<double> result(vec.size());
    auto iResult = std::minmax_element(vec.begin(), vec.end());
    const double sMin = *(iResult.first);
    const double sMax = *(iResult.second);
    double maxDiff = sMax - sMin;
    if (maxDiff < 1e-5) {
        result.resize(vec.size(), 0);
        return result;
    }    

    for (int i = 0; i < (int)result.size(); ++i) {
        result[i] = (vec[i] - sMin) / maxDiff;
    }
    return result;
}

void ColorSignature::curve(double sMin, double sMax)
{
    if (sMax <= sMin) {
        throw std::runtime_error("Color signature range error!");
    }

    const double maxDiff = sMax - sMin;
    mNormalizedValues.resize(mOriginalValues.size());
    for (int i = 0; i < (int)mOriginalValues.size(); ++i) {
        if (mOriginalValues[i] <= sMin) mNormalizedValues[i] = 0;
        else if (mOriginalValues[i] >= sMax) mNormalizedValues[i] = 1.;
        else mNormalizedValues[i] = (mOriginalValues[i] - sMin) / maxDiff;
    }
    mCurveMin = sMin;
    mCurveMax = sMax;
    changeColorMap(mCurrentColorMap);
}

void ColorSignature::setColorBuffer(float* dst) const
{
    for (size_t i = 0; i < mColors.size(); ++i) {
        for (int c = 0; c < 4; ++c) {
            dst[4 * i + c] = mColors[i][c];
        }
    }
}

}   // end of namespace