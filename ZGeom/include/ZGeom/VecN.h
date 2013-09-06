#ifndef ZGEOM_VECND_H
#define ZGEOM_VECND_H

#include <cassert>
#include <algorithm>
#include <vector>
#include "common.h"

#define USE_PPL

#ifdef USE_PPL
#include <ppl.h>
#endif

namespace ZGeom
{

template<typename T> class VecN;
template<typename T> class SparseMatrix;
template<typename T> VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

template<typename T>
class VecN
{
public:
    friend class SparseMatrix<T>;

	VecN() : mVec(NULL), mDim(0) {}
	VecN(uint n) : mVec(NULL) { resize(n); }
    VecN(uint n, T val) : mVec(NULL) { resize(n, val); }
    VecN(const VecN<T>& v2);
    VecN(T* vec, uint n);
    VecN(const std::vector<T>& vec);
	~VecN() { delete []mVec; }

    T operator [] (uint i) const { return mVec[i]; }
    T operator () (uint i) const { return mVec[i-1]; }
	T& operator [] (uint i) { return mVec[i]; }
    T& operator () (uint i) { return mVec[i-1]; }
    const VecN<T>& operator=(const VecN<T>& v2);

    T* c_ptr() const { return mVec; }
    T* c_ptr_end() const { return mVec + mDim; }
    std::vector<T> toStdVector() const { return std::vector<T>(mVec, mVec + mDim); }
    std::vector<T> operator() () const { return std::vector<T>(mVec, mVec + mDim); }
    uint size() const { return mDim; }
	void resize(int n);
    void resize(int n, T val);

    const VecN<T>& operator += (const VecN<T>& v2);
    const VecN<T>& operator -= (const VecN<T>& v2) { return (*this) += -v2; }
    const VecN<T>& operator += (T trans);
    const VecN<T>& operator -= (T trans) { return (*this) += -trans; }
    const VecN<T>& operator *= (T scale);
    const VecN<T>& operator /= (T scale) { return (*this) *= 1.0 / scale;}

    VecN<T> operator - () const;
    VecN<T> operator + (const VecN<T>& v2) const { VecN<T> v3(*this); v3 += v2; return v3; }
    VecN<T> operator - (const VecN<T>& v2) const { return *this + (-v2); }
    VecN<T> operator * (T coeff) const { VecN<T> v3(*this); v3 *= coeff; return v3; }
    VecN<T> operator / (T coeff) const { return (*this) * coeff; }
	friend VecN<T> operator * (T t, const VecN<T>& v1);    

    T dot(const VecN<T>& v2) const;
    T normEuclidean() const;
    friend VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

private:
	T *mVec;
	int mDim;
};

/* defining constructors */
template<typename T>
inline VecN<T>::VecN(const VecN<T>& v2) : mDim(v2.mDim)
{
    mVec = new T[mDim];
    for (int i = 0; i < mDim; ++i) mVec[i] = v2.mVec[i];
}

template<typename T>
inline VecN<T>::VecN(T* vec, uint n) 
    : mDim(n)
{
    mVec = new T[mDim];
    std::copy_n(vec, n, mVec);
}

template<typename T>
inline VecN<T>::VecN(const std::vector<T>& vec)
{
    mDim = vec.size();
    mVec = new T[mDim];
    std::copy_n(vec.begin(), mDim, mVec);
}
/* end of defining constructors */

template<typename T>
inline const VecN<T>& VecN<T>::operator = (const VecN<T>& v2)
{
    delete []mVec;
    mDim = v2.size();
    mVec = new T[mDim];
    std::copy_n(v2.c_ptr(), mDim, mVec);
    return *this;
}

template<typename T>
inline void VecN<T>::resize(int n)
{
	delete []mVec;
	this->mDim = n;
	mVec = new T[mDim];
}

template<typename T>
inline void VecN<T>::resize(int n, T val)
{
    resize(n);
    for (int i = 0; i < mDim; ++i) mVec[i] = val;
}

template<typename T>
inline VecN<T> VecN<T>::operator - () const
{
    VecN<T> v2(mDim);
    for (int i = 0; i < mDim; ++i) v2.mVec[i] = -mVec[i];
    return v2;
}

template<typename T>
inline VecN<T> operator * (T t, const VecN<T>& v1)
{
	return v1 * t;
}


/* definitions for operator+= */
#ifdef USE_PPL
template<typename T>
inline const VecN<T>& VecN<T>::operator += (const VecN<T>& v2)
{
    assert(mDim == v2.mDim);
    Concurrency::parallel_for(0, mDim, [&](int i) {
        mVec[i] += v2.mVec[i];
    });
    return *this;
}
template<typename T>
inline const VecN<T>& VecN<T>::operator += (T trans)
{
    Concurrency::parallel_for(0, mDim, [&](int i){
        mVec[i] += trans;
    });
    return *this;
}
#else
template<typename T>
inline const VecN<T>& VecN<T>::operator += (const VecN<T>& v2)
{
    assert(mDim == v2.mDim);
    for(int i = 0; i < mDim, ++i) {
        mVec[i] += v2.mVec[i];
    }
    return *this;
}
template<typename T>
inline const VecN<T>& VecN<T>::operator += (T trans)
{
    for(int i = 0; i < mDim; ++i){
        mVec[i] += trans;
    }
    return *this;
}
#endif


/* definitions for operator*= */
#ifdef USE_PPL
template<typename T>
inline const VecN<T>& VecN<T>::operator *= (T scale)
{
    Concurrency::parallel_for(0, mDim, [&](int i){
        mVec[i] *= scale;
    });
    return *this;
}
#else
template<typename T>
inline const VecN<T>& VecN<T>::operator *= (T scale)
{
    for(int i = 0; i < mDim; ++i){
        mVec[i] *= scale;
    }
    return *this;
}
#endif


/* definitions for dot product */
template<typename T>
inline T VecN<T>::dot(const VecN<T>& v2) const
{
    assert(mDim == v2.mDim);
    T sum(0.0);
    for (int i = 0; i < mDim; ++i) sum += mVec[i] * v2.mVec[i];
    return sum;
}

typedef VecN<float>  VecNs;
typedef VecN<double> VecNd;

}// end of namespace ZGeom

#endif