#pragma once

#include "MeshProcessor.h"

class WaveletCoefficients
{
public:
	
};

class WaveletMeshProcessor : public MeshProcessor
{
public:
	void computeMexicanHatWavelet(std::vector<double>& vMHW, int scale, int wtype = 1);
	void computeExperimentalWavelet(std::vector<double>& vExp, int scale);
	void computeDWTCoefficient(std::vector<double>& vCoeff, const std::vector<double>& vScales, const std::vector<double>& vfunc);
	void calGeometryDWT();
	void reconstructExperimental1(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
};