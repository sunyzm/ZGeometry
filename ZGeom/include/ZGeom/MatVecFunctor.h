#ifndef ZGEOM_MAT_VEC_MUL_FUNCTOR_H
#define ZGEOM_MAT_VEC_MUL_FUNCTOR_H

namespace ZGeom
{

class MatVecFunctor
{
public:
    virtual void operator() (double *in, double* out) = 0;
    int getOrder() const { return mOrder; }

protected:
    int mOrder;
};

}

#endif