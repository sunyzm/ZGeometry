#include "ColorSignature.h"
#include <algorithm>

namespace ZGeom {

ColorSignature::ColorSignature(const std::vector<Colorf>& vCol)
{
    mColors = vCol;
    mOriginalValues.resize(mColors.size());
    for (size_t i = 0; i < mColors.size(); ++i)
        mOriginalValues[i] = mColors[i].toGrayscale();
    
    mValues = mOriginalValues;
    mCurrentColorMap = CM_JET;
}

ColorSignature::ColorSignature(const std::vector<double>& vVals, ColorMapType cmt /*= CM_PARULA*/, bool requireNormalize /*= true*/)
{
    mValues = mOriginalValues = vVals;
    if (requireNormalize)
        mValues = vec_normalize(mOriginalValues);
    else curve(0, 1.);
    changeColorMap(cmt);
}

void ColorSignature::changeColorMap(ColorMapType cmt)
{
    mCurrentColorMap = cmt;

    if (mValues.empty()) return;
    mColors.resize(mValues.size());
    for (size_t i = 0; i < mValues.size(); ++i) 
        mColors[i].falseColor((float)mValues[i], 1.f, cmt);
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
    const double maxDiff = sMax - sMin;
    mValues.resize(mOriginalValues.size());
    for (int i = 0; i < (int)mOriginalValues.size(); ++i) {
        if (mOriginalValues[i] <= sMin) mValues[i] = 0;
        else if (mOriginalValues[i] >= sMax) mValues[i] = 1.;
        else mValues[i] = (mOriginalValues[i] - sMin) / maxDiff;
    }
    changeColorMap(mCurrentColorMap);
}

void ColorSignature::changeSignatureMode(SignatureMode smode)
{
    if (mValues.empty()) return;
    size_t vSize = mValues.size();
    std::vector<double> tmpSig = mValues;

    if (smode == ZGeom::SM_AbsNormalized)
        for (double& v : mValues) v = std::fabs(v);

    if (smode == ZGeom::SM_LogNormalized) {
        double sMin = *(std::minmax_element(mValues.begin(), mValues.end()).first);
        if (sMin <= 0) return;
        for (double& v : mValues) v = std::log(v);
    }

    if (smode == ZGeom::SM_PosNegPlot)
    {
        auto mmp = std::minmax_element(mValues.begin(), mValues.end());
        double sMin = *mmp.first, smax = *mmp.second;
        double absMax = std::max(std::fabs(sMin), std::fabs(smax));
        for (int i = 0; i < vSize; ++i)
            mColors[i].posNegColor(float(mValues[i] / absMax), ZGeom::ColorOrange, ZGeom::ColorAzure);
    }
    else if (smode == ZGeom::SM_Normalized || smode == ZGeom::SM_AbsNormalized || smode == ZGeom::SM_BandCurved ||
            smode == ZGeom::SM_LogNormalized || smode == ZGeom::SM_MarkNegNormalized)
    {
        changeColorMap(mCurrentColorMap);
    }

    if (smode == ZGeom::SM_MarkNegNormalized) {
        for (int i = 0; i < vSize; ++i) {
            if (mValues[i] < 0) mColors[i].setAs(ZGeom::ColorBlack);
        }
    }

    mValues = tmpSig;
}

}   // end of namespace