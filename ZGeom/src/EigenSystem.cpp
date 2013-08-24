#include "EigenSystem.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace ZGeom
{
    void EigenSystem::setSize(int order, int nev)
    {
        mOrder = order;
        mEvCount = nev;
        mEigVals.resize(mEvCount);
        mEigVecs.resize(mEvCount);
        for (int k = 0; k < mEvCount; ++k)
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
        for (int k = 0; k < mEvCount; ++k) {
            ofs1 << mEigVals[k] << '\n';
            for (int l = 0; l < mOrder; ++l) {
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

        for (int k = 0; k < mEvCount; ++k) {
            double *in = mEigVecs[k].c_ptr();
            (*A)(in, left);
            for (int i = 0; i < mOrder; ++i)
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

        for (int k = 0; k < mEvCount; ++k) {
            double* in = mEigVecs[k].c_ptr();
            (*A)(in, left);
            (*M)(in, right);
            for (int i = 0; i < mOrder; ++i)
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

} //end of namespace 
