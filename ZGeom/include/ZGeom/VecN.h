#ifndef ZGEOM_VECND_H
#define ZGEOM_VECND_H

#include <cassert>
#include <algorithm>
#include <vector>
#include "types.h"

#ifdef USE_MKL
#include <mkl.h>
#endif
#ifdef USE_PPL
#include <ppl.h>
#endif

namespace ZGeom
{

template<typename T>
class VecN
{
public:
	VecN() : mVec(NULL), mDim(0), mInternalStorage(true) {}
	VecN(uint n) : mVec(NULL), mInternalStorage(true) { resize(n); }
    VecN(uint n, T val) : mVec(NULL), mInternalStorage(true) { resize(n, val); }
    template<typename F> VecN(const VecN<F>& v2);
    VecN(T* vec, uint n, bool copyValues = true);
    VecN(const std::vector<T>& vec);

	~VecN() { if (mInternalStorage) delete []mVec; }
	T& operator [] (uint i) { return mVec[i]; }
    T& operator () (uint i) { return mVec[i-1]; }
    T& c_ptr() const { return mVec; } 
	void resize(int n);
    void resize(int n, T val);
    uint size() const { return mDim; }

    T normEuclidean() const;

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

    std::vector<T> toStdVector() const { return std::vector<T>(mVec, mVec + mDim); }
    std::vector<T> operator() () const { return std::vector<T>(mVec, mVec + mDim); }

private:
	T *mVec;
	int mDim;
    bool mInternalStorage;
};

/* defining constructors */
template<typename T>
template<typename F>
inline VecN<T>::VecN(const VecN<F>& v2) : mVec(NULL), mInternalStorage(true)
{
    resize(v2.mDim);
    for (int i = 0; i < mDim; ++i) mVec[i] = static_cast<T>(v2.mVec[i]);
}

template<typename T>
inline VecN<T>::VecN(T* vec, uint n, bool copyValues /* = true */) 
    : mVec(NULL), mInternalStorage(copyValues)
{
    if (mInternalStorage) {
        resize(n);
        std::copy(vec, vec + n, mVec);
    }
    else {
        mDim = n;
        mVec = vec;
    }
}

template<typename T>
inline VecN<T>::VecN(const std::vector<T>& vec) : mVec(NULL), mInternalStorage(true)
{
    resize(vec.size());
    std::copy(vec.begin(), vec.end(), mVec);
}
/* end of defining constructors */


template<typename T>
inline void VecN<T>::resize(int n)
{
	if (mInternalStorage) delete []mVec;
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
#ifdef USE_MKL
template<>
inline const VecN<double>& VecN<double>::operator *= (double scale)
{
    cblas_dscal(mDim, scale, mVec, 1);
    return *this;
}
template<>
inline const VecN<float>& VecN<float>::operator *= (float scale)
{
    cblas_sscal(mDim, scale, mVec, 1);
    return *this;
}
#else
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
#endif

/* definitions for dot product */
#ifdef USE_MKL
template<>
inline float VecN<float>::dot(const VecN<float>& v2) const
{
    assert(mDim == v2.mDim);
    return cblas_sdot(mDim, mVec, 1, v2.mVec, 1);
}
template<>
inline double VecN<double>::dot(const VecN<double>& v2) const
{
    assert(mDim == v2.mDim);
    return cblas_ddot(mDim, mVec, 1, v2.mVec, 1);
}
#else
template<typename T>
inline T VecN<T>::dot(const VecN<T>& v2) const
{
    assert(mDim == v2.mDim);
    T sum(0.0);
    for (int i = 0; i < mDim; ++i) sum += mVec[i] * v2.mVec[i];
    return sum;
}
#endif

typedef VecN<double> VecNd;
typedef VecN<float>  VecNf;

}// end of namespace ZGeom



#endif