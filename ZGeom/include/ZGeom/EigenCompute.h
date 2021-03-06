#ifndef ZGEOM_EIGEN_COMPUTE_H
#define ZGEOM_EIGEN_COMPUTE_H
#include "MatlabEngineWrapper.h"
#include "EigenSystem.h"
#include "SparseMatrix.h"

namespace ZGeom {

class EigenCompute
{
public:
    EigenCompute(MatlabEngineWrapper* engineWrapper) : m_ep(engineWrapper) {}
    void solveStdSym(const SparseMatrix<double>& matA, int nEig, EigenSystem& eigSys);
    void solveGenSym(const SparseMatrix<double>& matA, const SparseMatrix<double>& matB, int nEig, EigenSystem& eigSys);

private:
    MatlabEngineWrapper* m_ep;
};

}	// end of namespace

#endif