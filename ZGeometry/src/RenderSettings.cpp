#include <algorithm>
#include <cassert>
#include <tuple>
#include "RenderSettings.h"
using namespace std;


void RenderSettings::normalizeSignatureFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	

	importSignature(vFrom);
	normalizeSignature();
}

void RenderSettings::logNormalizeSignatureFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	

	importSignature(vFrom);
	logNormalizeSignature();
}

void RenderSettings::bandCurveSignatureFrom( const std::vector<double>& vFrom, double lowend, double highend )
{
	assert(lowend < highend);
	if (vFrom.empty()) return;	

	importSignature(vFrom);
	bandCurveSignature(lowend, highend);
}

void RenderSettings::importSignature( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;

	vOriginalSignature = vFrom;
	auto iResult = minmax_element(vOriginalSignature.begin(), vOriginalSignature.end());
	sigMin = *(iResult.first); 
	sigMax = *(iResult.second);
}

void RenderSettings::normalizeSignature()
{
	assert(!vOriginalSignature.empty());

	const double sMin = sigMin, sMax = sigMax;

	vDisplaySignature.resize(vOriginalSignature.size());
	std::transform(vOriginalSignature.begin(), vOriginalSignature.end(), vDisplaySignature.begin(),
		[=](double v){ return (v-sMin)/(sMax-sMin); });
}

void RenderSettings::logNormalizeSignature()
{
	assert(!vOriginalSignature.empty());

	std::vector<double> vLog(vOriginalSignature.size());
	std::transform(vOriginalSignature.begin(), vOriginalSignature.end(), vLog.begin(), [](double v){ return std::log(v+1); });
	const double sMin = std::log(sigMin+1), sMax = std::log(sigMax+1);

	vDisplaySignature.resize(vLog.size());
	std::transform(vLog.begin(), vLog.end(), vDisplaySignature.begin(),
		[=](double v){ return (v-sMin)/(sMax-sMin); });
}

void RenderSettings::bandCurveSignature(double lowend, double highend)
{
	assert(!vOriginalSignature.empty());
	assert(lowend < highend);

	this->vDisplaySignature.resize(vOriginalSignature.size());
	transform(vOriginalSignature.begin(), vOriginalSignature.end(), vDisplaySignature.begin(), [=](double v)->double
	{
		if (v <= lowend)
			return 0.0;
		else if (v >= highend)
			return 1.0;
		else 
			return (v - lowend) / (highend - lowend);
	});
}