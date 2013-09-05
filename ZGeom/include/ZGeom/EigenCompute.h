#include "MatlabEngineWrapper.h"
#include "EigenSystem.h"
#include "SparseMatrix.h"

namespace ZGeom
{
    class EigenCompute
    {
    public:
        EigenCompute(const MatlabEngineWrapper* engineWrapper) : m_ep(engineWrapper) {}
        void solveStdSym(const SparseMatrix<double>& matA, int nEig, EigenSystem& eigSys);
        void solveGenSym(const SparseMatrix<double>& matA, const SparseMatrix<double>& matB, int nEig, EigenSystem& eigSys);

    private:
        const MatlabEngineWrapper* m_ep;
    };
}