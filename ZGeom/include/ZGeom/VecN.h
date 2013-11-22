#ifndef ZGEOM_VECN_H
#define ZGEOM_VECN_H
#include <cassert>
#include <algorithm>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <iterator>
#include <functional>
#include "common.h"

//#define USE_PPL
#ifdef USE_PPL
#include <ppl.h>
#endif

namespace ZGeom
{

template<typename T> class VecN;
template<typename T> class SparseMatrix;
template<typename T> VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);
typedef std::function<double(const VecN<double>&, const VecN<double>&)> InnerProdcutFunc;

template<typename T>
class VecN
{
public:
	friend class SparseMatrix<T>;
	class iterator;

	VecN() : mVec(NULL), mDim(0) {}
	VecN(const VecN<T>& v2);
	VecN(VecN<T>&& v2);
	VecN(uint n) : mVec(NULL), mDim(0) { resize(n, 0); }
	VecN(uint n, T val) : mVec(NULL), mDim(0) { resize(n, val); }
	VecN(T* vec, uint n);
	VecN(const std::vector<T>& vec);
	~VecN() { delete []mVec; }
	const VecN<T>& operator = (const VecN<T>& v2);
	VecN<T>& operator = (VecN<T>&& v2);

	T operator [] (uint i) const { return mVec[i]; }
	T operator () (uint i) const { return mVec[i-1]; }
	T& operator [] (uint i) { return mVec[i]; }
	T& operator () (uint i) { return mVec[i-1]; }
	T at(int i) const { 
		if (i < 0 || i >= mDim) throw std::runtime_error("Invalid VecN subscript!"); 
		return mVec[i]; 
	}
	T& at(int i) {
		if (i < 0 || i >= mDim) throw std::runtime_error("Invalid VecN subscript!");
		return mVec[i]; 
	}

	T* c_ptr() const { return mVec; }
	T* c_ptr_end() const { return mVec + mDim; }
	std::vector<T> toStdVector() const { return std::vector<T>(mVec, mVec + mDim); }
	std::vector<T> operator() () const { return std::vector<T>(mVec, mVec + mDim); }
	uint size() const { return mDim; }
	void resize(int n);
	void resize(int n, T val);
	void expandTo(int n);
	void copyElements(const VecN<T>& v2, int startingPos = 0);

	const VecN<T>& add(T* v2);

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
	friend VecN<T> operator * (T t, const VecN<T>& v1) {
		return v1 * t;
	}

	T norm1() const;
	T norm2() const;
	T normEuclidean() const { return norm2(); }
	T pNorm(double p) const;
	T inftyNorm() const;

	T dot(const VecN<T>& v2) const;
	T sum() const;
	T partial_sum(int a, int b) const;
	friend VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

	void normalize(double p);
	void normalize(const InnerProdcutFunc& innerProdFunc);

	/* iterator operations */
	iterator begin() { return iterator(this, 0); };
	iterator end() { return iterator(this, mDim); }

private:
	T *mVec;
	int mDim;
};

template<typename T>
class VecN<T>::iterator : public std::iterator<std::bidirectional_iterator_tag, T>
{
public:
	iterator() : mVecN(nullptr), mPos(0) {}
	iterator(const iterator& iter) : mVecN(iter.mVecN), mPos(iter.mPos) {}

	iterator(VecN<T>* vec, int pos) : mVecN(vec) { 
		if (pos < 0 || pos > mVecN->mDim) throw std::runtime_error("Invalid VecN subscript"); 
		mPos = pos;
	}

	iterator& operator = (const iterator iter) { 
		mVecN = iter.mVecN; 
		mPos = iter.mPos;
		return *this; 
	}

	bool operator == (const iterator& iter) const {
		return mVecN == iter.mVecN && mPos == iter.mPos;
	}

	bool operator != (const iterator& iter) const {
		return mVecN != iter.mVecN || mPos != iter.mPos;
	}

	T& operator * () {  return mVecN->at(mPos); }
	bool operator () () { return mPos < mVecN->mDim; }

	iterator operator ++ (int) { //postfix ++
		VecN<T>::iterator iter(*this); 
		if (mPos < mVecN->mDim) mPos++; 
		return iter; 
	}

	iterator& operator ++ () { //prefix --	
		if (mPos < mVecN->mDim) ++mPos; 
		return *this; 
	}

	iterator operator -- (int) { //postfix --
		iterator iter(*this); 
		if (mPos > 0) mPos--;
		return iter; 
	}

	iterator& operator -- () { //prefix --
		if (mPos > 0) --mPos; 
		return *this; 
	}

private:
	VecN<T>* mVecN;
	int mPos;
};

/* defining constructors and assignment operator */
template<typename T>
inline VecN<T>::VecN(const VecN<T>& v2) : mDim(v2.mDim)
{
	mVec = new T[mDim];
	std::copy_n(v2.mVec, v2.mDim, mVec);
}

template<typename T>
inline VecN<T>::VecN(VecN<T>&& v2) : mDim(0), mVec(nullptr)
{
	mDim = v2.mDim;
	mVec = v2.mVec;

	v2.mDim = 0;
	v2.mVec = nullptr;
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
	mDim = (int)vec.size();
	mVec = new T[mDim];
	std::copy_n(vec.begin(), mDim, mVec);
}

template<typename T>
inline const VecN<T>& VecN<T>::operator = (const VecN<T>& v2)
{
	if (this != &v2) {
		delete []mVec;
		mDim = v2.mDim;
		mVec = new T[mDim];
		std::copy_n(v2.c_ptr(), mDim, mVec);
	}
	return *this;
}

template<typename T>
inline VecN<T>& VecN<T>::operator = (VecN<T>&& v2)
{
	if (this != &v2) {
		delete []mVec;
		mDim = v2.mDim;
		mVec = v2.mVec;
		v2.mDim = 0;
		v2.mVec = nullptr;
	}
	return *this;
}
template<typename T>
inline void VecN<T>::resize(int n)
{
	if (n == mDim && mVec != NULL) return;

	delete []mVec;
	mDim = n;
	mVec = new T[mDim];
}

template<typename T>
inline void VecN<T>::resize(int n, T val)
{
	resize(n);
	for (int i = 0; i < mDim; ++i) mVec[i] = val;
}

/* end of defining constructors and assignment operator */



template<typename T>
inline void VecN<T>::copyElements( const VecN<T>& v2, int startingPos )
{
	if (mDim - startingPos < v2.size())
		throw std::runtime_error("Error VecN::copyElements: Insufficient size in destination VecN");
	std::copy_n(v2.mVec, v2.size(), mVec + startingPos);
}

template<typename T>
inline void VecN<T>::expandTo( int n )
{
	if (n < mDim) throw runtime_error("VecN cannot be expanded to a smaller size than before");
	else if (n == mDim) return;
	else {
		T* data = new T[n];
		if (mDim > 0) std::copy_n(mVec, mDim, data);
		delete []mVec;
		mVec = data;
		mDim = n;
	}
}

template<typename T>
inline VecN<T> VecN<T>::operator - () const
{
	VecN<T> v2(mDim);
	for (int i = 0; i < mDim; ++i) v2.mVec[i] = -mVec[i];
	return v2;
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
	for(int i = 0; i < mDim; ++i) {
		mVec[i] += v2.mVec[i];
	}
	return *this;
}
template<typename T>
inline const VecN<T>& VecN<T>::operator += (T trans)
{
	for(int i = 0; i < mDim; ++i) {
		mVec[i] += trans;
	}
	return *this;
}
#endif

template<typename T>
inline const VecN<T>& VecN<T>::add(T* v2)
{
	for (int i = 0; i < mDim; ++i)
		mVec[i] += v2[i];
	return *this;
}

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
	for(int i = 0; i < mDim; ++i) {
		mVec[i] *= scale;
	}
	return *this;
}
#endif

template<typename T>
inline T VecN<T>::norm1() const 
{
	T sum(0);
	for (int i = 0; i < mDim; ++i) sum += std::fabs(mVec[i]);
	return sum;
}

template<typename T>
inline T VecN<T>::norm2() const 
{
	T sum(0);
	for (int i = 0; i < mDim; ++i) sum += std::pow(mVec[i], 2);
	return std::sqrt(sum);
}

template<typename T>
inline T VecN<T>::pNorm( double p ) const
{
	if( p < 0) return inftyNorm();

	double sum(0);
	for (int i = 0; i < mDim; ++i) sum += std::pow(std::fabs(mVec[i]), p);
	return std::pow(sum, 1.0/p);
}


template<typename T>
inline T ZGeom::VecN<T>::inftyNorm() const
{
	double maxNorm(0);
	for (int i = 0; i < mDim; ++i) 
		if (std::fabs(mVec[i]) > maxNorm) maxNorm = std::fabs(mVec[i]);
	return maxNorm;
}


template<typename T>
inline void VecN<T>::normalize(double p)
{
	assert(p > 0);
	double norm = this->pNorm(p);
	for (int i = 0; i < mDim; ++i) mVec[i] /= norm;
}

template<>
inline void VecN<double>::normalize(const InnerProdcutFunc& innerProdFunc)
{
	double norm = std::sqrt(innerProdFunc(*this, *this));
	for (int i = 0; i < mDim; ++i) mVec[i] /= norm;
}

/* definitions for dot product */
template<typename T>
inline T VecN<T>::dot(const VecN<T>& v2) const
{
	assert(mDim == v2.mDim);	
	return std::inner_product(mVec, mVec + mDim, v2.mVec, 0.0);
}

template<typename T>
inline T VecN<T>::sum() const
{
	double retval(0);
	for (int i = 0; i < mDim; ++i) retval += mVec[i];
	return retval;
}

template<typename T>
inline T VecN<T>::partial_sum(int a, int b) const 
{
	assert(a >= 0 && a <= b && b <= mDim);
	return std::accumulate(mVec + a, mVec + b, 0);
}

typedef VecN<float>  VecNs;
typedef VecN<double> VecNd;

}// end of namespace ZGeom

#endif