#include <algorithm>
#include <cassert>
#include "RenderSettings.h"
using namespace std;


void RenderSettings::normalizeSignatureFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	

	auto iResult = minmax_element(vFrom.begin(), vFrom.end());
	const double sMin = *iResult.first, sMax = *iResult.second;

	sigMin = sMin, sigMax = sMax;

	vDisplaySignature.resize(vFrom.size());
	std::transform(vFrom.begin(), vFrom.end(), vDisplaySignature.begin(),
		[=](double v){ return (v-sMin)/(sMax-sMin); });

// 	this->vDisplaySignature.clear();
// 	for (vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
// 	{
// 		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
// 	}
}

void RenderSettings::logNormalizeSignatureFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	

	std::vector<double> vLog(vFrom.size());
	std::transform(vFrom.begin(), vFrom.end(), vLog.begin(), [](double v){ return std::log(v+1); } );

// 	vLog.resize(vFrom.size());
// 	for (std::vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
// 	{
// 		vLog.push_back(std::log(*iter + 1));
// 	}

	auto iResult = minmax_element(vLog.begin(), vLog.end());
	const double sMin = *iResult.first, sMax = *iResult.second;
	sigMin = sMin, sigMax = sMax;

	vDisplaySignature.resize(vLog.size());
	std::transform(vLog.begin(), vLog.end(), vDisplaySignature.begin(),
		[=](double v){ return (v-sMin)/(sMax-sMin); });

// 	this->vDisplaySignature.clear();
// 	this->vDisplaySignature.reserve(vLog.size());
// 	for (std::vector<double>::const_iterator iter = vLog.begin(); iter != vLog.end(); ++iter)
// 	{
// 		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
// 	}
}

void RenderSettings::bandCurveSignatureFrom( const std::vector<double>& vFrom, double lowend, double highend )
{
	assert(lowend < highend);

	sigMin = lowend, sigMax = highend;

	this->vDisplaySignature.clear();
	this->vDisplaySignature.reserve(vFrom.size());
	for (std::vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		if (*iter <= lowend)
			vDisplaySignature.push_back(0.0);
		else if (*iter >=highend)
			vDisplaySignature.push_back(1.0);
		else 
			vDisplaySignature.push_back((*iter - lowend)/(highend - lowend));
	}
}
