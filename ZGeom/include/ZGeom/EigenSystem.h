#ifndef ZGEOM_EIGEN_SYSTEM_H
#define ZGEOM_EIGEN_SYSTEM_H

#include <vector>
#include "VecN.h"

namespace ZGeom
{
    class EigenSystem
    {
    public:


    private:
        int mOrder;
        int mEigCount;
        std::vector<double> mEigVals;
        std::vector<VecNd>  mEigVecs;
    };
}



#endif