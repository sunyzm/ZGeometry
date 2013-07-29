#pragma once
#include <engine.h>
#include <vector>

class MatlabWrapper
{
public:
	void setEngine(Engine* ep) { m_ep = ep; }
	void SparseBiConjugateGradient(int m_size, int n_size, const std::vector<int>& vI, const std::vector<int>& vJ, const std::vector<double>& vS, const std::vector<double>& vY, std::vector<double>& vX);
	void SparseMatVecMul(int m_size, int n_size, const std::vector<int>& vI, const std::vector<int>& vJ, const std::vector<double>& vS, const std::vector<double>& vX, std::vector<double>& vY);	// AX=Y
	void DenseConjugateGradient(const std::vector<std::vector<double> >& A, const std::vector<double>& b, std::vector<double>& sol);
	
private:
	Engine* m_ep;
};

void matlab_cgls(Engine* ep, const std::vector<std::vector<double> >& SGW, const std::vector<double>& b, std::vector<double>& sol);

void matlab_scgls(Engine* ep, const std::vector<std::vector<double> >& SGW, const std::vector<double>& b, std::vector<double>& sol);