#include "EigenSystem.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace ZGeom {
void EigenSystem::setSize(int order, int nev)
{
	mOrder = order;
	mEvCount = nev;
	mEigVals.resize(mEvCount);
	mEigVecs.resize(mEvCount);
	for (uint k = 0; k < mEvCount; ++k)
		mEigVecs[k].resize(mOrder);
}

void EigenSystem::setValues(int index, double eigVal, double eigVec[])
{
	mEigVals[index] = eigVal;
	std::copy(eigVec, eigVec + mOrder, mEigVecs[index].c_ptr());
}

void EigenSystem::setValues(int index, double eigVal, const std::vector<double>& eigVec)
{
	mEigVals[index] = eigVal;
	std::copy(eigVec.begin(), eigVec.end(), mEigVecs[index].c_ptr());
}

void EigenSystem::print(const std::string& file1, const std::string& file2) const
{
	std::ofstream ofs1(file1.c_str()), ofs2(file2.c_str());
	for (uint k = 0; k < mEvCount; ++k) {
		ofs1 << mEigVals[k] << '\n';
		for (uint l = 0; l < mOrder; ++l) {
			ofs2 << mEigVecs[k][l] << ' ';
		}
		ofs2 << '\n';
	}
	ofs1.close();
	ofs2.close();
}

void EigenSystem::inverseEigVals()
{
	for (uint k = 0; k < mEigVals.size(); ++k) mEigVals[k] = 1.0 / mEigVals[k];
}

void EigenSystem::evalError(MatVecFunctor* A, std::vector<double>& vErrors) const
{
	double *in = new double[mOrder];
	double *left = new double[mOrder];
	double *diff = new double[mOrder];
	vErrors.resize(mEvCount);

	std::cout << "Standard eigen-decomposition error: ";

	for (uint k = 0; k < mEvCount; ++k) {
		double *in = mEigVecs[k].c_ptr();
		(*A)(in, left);
		for (uint i = 0; i < mOrder; ++i)
			diff[i] = std::abs(left[i] - mEigVals[k] * in[i]);
		double maxDiff = *std::max_element(diff, diff + mOrder);
		vErrors[k] = maxDiff;
		std::cout << maxDiff << ' ';
	}

	std::cout << std::endl;

	delete []in;
	delete []left;
	delete []diff;
}

void EigenSystem::evalGeneralError(MatVecFunctor* A, MatVecFunctor* M, std::vector<double>& vErrors) const
{
	double *in = new double[mOrder];
	double *left = new double[mOrder];
	double *right = new double[mOrder];
	double *diff = new double[mOrder];
	vErrors.resize(mEvCount);

	std::cout << "General eigen-decomposition error: ";

	for (uint k = 0; k < mEvCount; ++k) {
		double* in = mEigVecs[k].c_ptr();
		(*A)(in, left);
		(*M)(in, right);
		for (uint i = 0; i < mOrder; ++i)
			diff[i] = std::abs(left[i] - mEigVals[k] * right[i]);
		double maxDiff = *std::max_element(diff, diff + mOrder);
		vErrors[k] = maxDiff;
		std::cout << maxDiff << ' ';
	}

	std::cout << std::endl;

	delete []in;
	delete []left;
	delete []right;
	delete []diff;
}

void EigenSystem::printEigVals( const std::string& file1 ) const
{
	std::ofstream ofs1(file1.c_str());
	for (uint k = 0; k < mEvCount; ++k) {
		ofs1 << mEigVals[k] << '\n';
	}
	ofs1.close();
}

void EigenSystem::save(const std::string& file) const
{
	std::ofstream ofs(file.c_str(), std::ios::trunc|std::ios::binary);
	ofs.write((char*)&mEvCount, sizeof(int));
	ofs.write((char*)&mOrder, sizeof(int));

	ofs.write((char*)&mEigVals[0], sizeof(double)*mEvCount);
	for (uint i = 0; i < mEvCount; ++i)
		ofs.write((char*)mEigVecs[i].c_ptr(), sizeof(double)*mOrder);

	ofs.close();
}

void EigenSystem::load(const std::string& file)
{
	std::ifstream ifs(file.c_str(), std::ios::binary);
	ifs.read((char*)&mEvCount, sizeof(int));
	ifs.read((char*)&mOrder, sizeof(int));

	setSize(mOrder, mEvCount);
	ifs.read((char*)&mEigVals[0], sizeof(double)*mEvCount);
	for (uint i = 0; i < mEvCount; ++i)
		ifs.read((char*)mEigVecs[i].c_ptr(), sizeof(double)*mOrder);

	ifs.close();
}

double EigenSystem::heatKernel( uint x, uint y, double t ) const
{
	double sum(0);
	for (uint k = 0; k < mEvCount; ++k) {
		double lambda = mEigVals[k];
		double* phi = mEigVecs[k].c_ptr();
		double term = std::exp(-lambda * t) * phi[x] * phi[y];
		sum += term;
	}
	return sum;
}

double EigenSystem::kernel( uint x, uint y, double t, std::function<double(double, double)> transferFunc ) const
{
	double sum(0);
	for (uint k = 0; k < mEvCount; ++k) {
		double lambda = mEigVals[k];
		double *phi = mEigVecs[k].c_ptr();
		double term = transferFunc(lambda, t) * phi[x] * phi[y];
		sum += term;
	}
	return sum;
}

double EigenSystem::kernel( uint x, uint y, std::function<double(double)> transferFunc ) const
{
	double sum(0);
	for (uint k = 0; k < mEvCount; ++k) {
		double lambda = mEigVals[k];
		double *phi = mEigVecs[k].c_ptr();
		double term = transferFunc(lambda) * phi[x] * phi[y];
		sum += term;
	}
	return sum;
}

void EigenSystem::ToDenseMatrix(ZGeom::DenseMatrixd& matEigVec) const
{
	matEigVec.resize(mEvCount, mOrder);
	double *pMatEigVec = matEigVec.raw_ptr();
	for (int i = 0; i < (int)mEvCount; ++i)
		matEigVec.setRow(i, mEigVecs[i]);
}


}	// end of namespace 
