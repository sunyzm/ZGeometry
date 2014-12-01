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
	for (int k = 0; k < mEvCount; ++k) mEigVecs[k].resize(mOrder);
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
	for (int k = 0; k < eigVecCount(); ++k) {
		ofs1 << mEigVals[k] << '\n';
		for (int l = 0; l < mOrder; ++l) {
			ofs2 << mEigVecs[k][l] << ' ';
		}
		ofs2 << '\n';
	}
	ofs1.close();
	ofs2.close();
}

void EigenSystem::printEigVals( const std::string& file1 ) const
{
	std::ofstream ofs1(file1.c_str());
	for (int k = 0; k < eigVecCount(); ++k) {
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
	for (int i = 0; i < mEvCount; ++i)
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
	for (int i = 0; i < mEvCount; ++i)
		ifs.read((char*)mEigVecs[i].c_ptr(), sizeof(double)*mOrder);

	ifs.close();
}

ZGeom::DenseMatrixd EigenSystem::toDenseMatrix() const
{
    DenseMatrixd matEigVec;
	matEigVec.resize(eigVecCount(), mOrder);
	for (int i = 0; i < (int)eigVecCount(); ++i)
		matEigVec.setRow(i, mEigVecs[i]);

    return matEigVec;
}


}	// end of namespace 
