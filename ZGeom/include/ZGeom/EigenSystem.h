#ifndef ZGEOM_EIGEN_SYSTEM_H
#define ZGEOM_EIGEN_SYSTEM_H

#include <vector>
#include "VecN.h"
#include "MatVecFunctor.h"

namespace ZGeom
{
    class EigenSystem
    {
    public:
        friend class EigenCompute;

        void setSize(int order, int nev);
        void setValues(int index, double eigVal, double eigVec[]);
        void setValues(int index, double eigVal, const std::vector<double>& eigVec);
        double getEigVal(int index) const { return mEigVals[index]; }
        const VecNd& getEigVec(int index) const { return mEigVecs[index]; }
        VecNd& getEigVec(int index) { return mEigVecs[index]; }
        void print(const std::string& file1, const std::string& file2) const;

        void evalError(MatVecFunctor* A, std::vector<double>& vErrors) const;
        void evalGeneralError(MatVecFunctor* A, MatVecFunctor* M, std::vector<double>& vErrors) const;
        void inverseEigVals();

    private:
        int mOrder;
        int mEvCount;
        std::vector<double> mEigVals;
        std::vector<VecNd> mEigVecs;
    };
}



#endif