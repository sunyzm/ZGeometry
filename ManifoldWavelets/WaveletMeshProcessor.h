#pragma once

#include "MeshProcessor.h"

class WaveletMeshProcessor : public MeshProcessor
{
private:
	std::vector<double> m_vTimescales;
	std::vector<std::vector<double> > m_vSGW;
public:
	void computeMexicanHatWavelet(std::vector<double>& vMHW, double scale, int wtype = 1);
	void computeExperimentalWavelet(std::vector<double>& vExp, double scale);
	void computeDWTCoefficient(std::vector<double>& vCoeff, const std::vector<double>& vScales, const std::vector<double>& vfunc);
	void calGeometryDWT();
	void reconstructExperimental1(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructBySGW(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false);
	void reconstructByDifferential(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint = false) const;
	void reconstructByMHB(int approxN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
	void computeSGW(const std::vector<double>& timescales);
	void getSGWSignature(double timescale, std::vector<double>& values) const;
};