#include <algorithm>
#include <cassert>
#include "RenderSettings.h"
using namespace std;


void RenderSettings::normalizeSignatureFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	
	assert(vFrom.size() == m_size);

	auto iResult = minmax_element(vFrom.begin(), vFrom.end());
	double sMin = *iResult.first, sMax = *iResult.second;

	sigMin = sMin, sigMax = sMax;

	this->vDisplaySignature.clear();
	for (vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
	}
}

void RenderSettings::logNormalizeSignatureFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	
	assert(vFrom.size() == m_size);

	std::vector<double> vLog;
	vLog.reserve(vFrom.size());
	for (std::vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		vLog.push_back(std::log(*iter + 1));
	}

	auto iResult = minmax_element(vLog.begin(), vLog.end());
	double sMin = *iResult.first, sMax = *iResult.second;
	sigMin = sMin, sMax = sMax;

	this->vDisplaySignature.clear();
	this->vDisplaySignature.reserve(vLog.size());
	for (std::vector<double>::const_iterator iter = vLog.begin(); iter != vLog.end(); ++iter)
	{
		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
	}
}

void RenderSettings::bandCurveSignatureFrom( const std::vector<double>& vFrom, double lowend, double highend )
{
	assert(lowend < highend);
	assert(vFrom.size() == m_size);

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
