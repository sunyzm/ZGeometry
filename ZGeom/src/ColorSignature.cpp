#include "ColorSignature.h"
#include <algorithm>

namespace ZGeom {

ColorSignature::ColorSignature(const std::vector<Colorf>& vCol)
{
    mColors = vCol;
    mValues.resize(mColors.size());
    for (size_t i = 0; i < mColors.size(); ++i)
        mValues[i] = mColors[i].toGrayscale();
    
    mCurrentColorMap = CM_PARULA;
}

ColorSignature::ColorSignature(const std::vector<double>& vVals, ColorMapType cmt /*= CM_PARULA*/, bool requireNormalize /*= true*/)
{
    mValues = vVals;
    changeColorMap(cmt, requireNormalize);
}

void ColorSignature::changeColorMap(ColorMapType cmt, bool requireNormalize)
{
    if (mValues.empty()) return;
    if (requireNormalize) {
        auto iResult = std::minmax_element(mValues.begin(), mValues.end());
        double sMin = *(iResult.first);
        double sMax = *(iResult.second);
        mColors.resize(mValues.size());

        for (size_t i = 0; i < mValues.size(); ++i) 
            mColors[i].falseColor(float((mValues[i] - sMin) / (sMax - sMin)), 1.f, cmt);        
    }
    else {
        for (size_t i = 0; i < mValues.size(); ++i) 
            mColors[i].falseColor((float)mValues[i], 1.f, cmt);
    }

    mCurrentColorMap = cmt;
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

}