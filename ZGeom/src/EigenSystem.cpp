#include "EigenSystem.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include "dataio.h"
#include "MatlabEngineWrapper.h"
#include "MatVecArithmetic.h"

using namespace std::tr2::sys;

namespace ZGeom {

void EigenSystem::setSize(int order, int nev)
{
	mEigVals.resize(nev);
	mEigVecs.resize(nev);
    for (auto& vk : mEigVecs) {
        vk.resize(order);
    }
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
    DenseMatrixd mat_eig_vec;
    // eigenvectors arranged row by row
	mat_eig_vec.resize(eigVecCount(), eigVecSize());
    for (int i = 0; i < (int)eigVecCount(); ++i) {
        mat_eig_vec.setRow(i, mEigVecs[i]);
    }
    return mat_eig_vec;
}

void EigenSystem::resize(int m)
{
    if (m >= mEigVals.size() || m < 0) return;
    mEigVals.erase(mEigVals.begin() + m, mEigVals.end());
    mEigVecs.erase(mEigVecs.begin() + m, mEigVecs.end());
}

const ZGeom::DenseMatrixd& EigenSystem::getEigenMat()
{
    if (mEigenMat.empty()) mEigenMat = toDenseMatrix();
    return mEigenMat;
}

void EigenSystem::clear()
{
    mEigVals.clear(); mEigVecs.clear(); mEigenMat.clear();
}

void EigenSystem::computeEigenMat()
{
    mEigenMat.clear(); 
    mEigenMat = toDenseMatrix();
}

void EigenSystem::validate()
{
    std::vector<double> v0 = getEigVec(1).toStdVector(), v1 = getEigVec(10).toStdVector();
    
    if (hasInducingMat()) {
        std::vector<double> diag = getInducingMat().getDiagonal();
        std::cout << "<v0,v1> = " << inducedInnerProduct(v0, v1, diag) << std::endl;
        std::cout << "<v0,v0> = " << inducedInnerProduct(v0, v0, diag) << std::endl;
        std::cout << "<v1,v1> = " << inducedInnerProduct(v1, v1, diag) << std::endl;
    }
    else {
        std::cout << "<v0,v1> = " << innerProduct(v0, v1) << std::endl;
        std::cout << "<v0,v0> = " << innerProduct(v0, v0) << std::endl;
        std::cout << "<v1,v1> = " << innerProduct(v1, v1) << std::endl;
    }
}

void EigenSystem::saveToMat(MatlabEngineWrapper* eng, const std::string& eigen_file_path) const
{
    eng->addDoubleScalar(eigVecCount(), "eigen_num");
    eng->addDoubleScalar(eigVecSize(), "mesh_size");
    eng->addColVec(VecNd(getAllEigVals()), "eigenvals");
    eng->addDenseMat(hasEigenMat() ? getEigenMat_const() : toDenseMatrix(), "eigenvecs", true);
    eng->addSparseMat(getInducingMat(), "inducing_mat");
    eng->eval("save " + eigen_file_path + " eigen_num mesh_size eigenvals eigenvecs inducing_mat");
}

void EigenSystem::loadFromMat(MatlabEngineWrapper* eng, const std::string& file)
{
    if (!fileExist(file)) throw std::runtime_error("Eigen mat file not exist!");
    clear();    

    auto pwd = initial_path<path>();
    eng->eval("cd " + pwd.string());

    eng->eval("load " + file);
    int eigen_num = int(*eng->getDblVariablePtr("eigen_num"));
    int mesh_size = int(*eng->getDblVariablePtr("mesh_size"));
    setSize(mesh_size, eigen_num);

    double* eig_vals_ptr = eng->getDblVariablePtr("eigenvals");
    for (int k = 0; k < eigen_num; ++k) this->mEigVals[k] = eig_vals_ptr[k];

    this->mEigenMat = eng->getDenseMat("eigenvecs", true);
    for (int k = 0; k < eigen_num; ++k) {
        mEigVecs[k] = mEigenMat.getRowVec(k);
    }
    
    eng->getSparseMat("inducing_mat", this->mInducingMat);
}

}	// end of namespace 
