#include "EigenSystem.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace ZGeom {

void EigenSystem::setSize(int order, int nev)
{
	mEigVals.resize(nev);
	mEigVecs.resize(nev);
    for (auto& vk : mEigVecs) vk.resize(order);
}

void EigenSystem::print(const std::string& file1, const std::string& file2) const
{
	std::ofstream ofs1(file1.c_str()), ofs2(file2.c_str());
	for (int k = 0; k < eigVecCount(); ++k) {
		ofs1 << mEigVals[k] << '\n';
		for (int l = 0; l < eigVecSize(); ++l) {
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
    int nOrder = eigVecSize();
    int nev = eigVecCount();

	std::ofstream ofs(file.c_str(), std::ios::trunc|std::ios::binary);
	ofs.write((char*)&nev, sizeof(int));
	ofs.write((char*)&nOrder, sizeof(int));
	ofs.write((char*)mEigVals.data(), sizeof(double)*nev);
	for (int i = 0; i < nev; ++i)
		ofs.write((char*)mEigVecs[i].c_ptr(), sizeof(double)*nOrder);
	ofs.close();
}

void EigenSystem::load(const std::string& file)
{
    int nev, nOrder;
	std::ifstream ifs(file.c_str(), std::ios::binary);
	ifs.read((char*)&nev, sizeof(int));
	ifs.read((char*)&nOrder, sizeof(int));

	setSize(nOrder, nev);
	ifs.read((char*)&mEigVals[0], sizeof(double)*nev);	
    for (int i = 0; i < nev; ++i)
		ifs.read((char*)mEigVecs[i].c_ptr(), sizeof(double)*nOrder);
	ifs.close();
}

ZGeom::DenseMatrixd EigenSystem::toDenseMatrix() const
{
    DenseMatrixd matEigVec;
	matEigVec.resize(eigVecCount(), eigVecSize());
	for (int i = 0; i < (int)eigVecCount(); ++i)
		matEigVec.setRow(i, mEigVecs[i]);

    return matEigVec;
}


}	// end of namespace 
