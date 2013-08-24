#ifndef ZGEOM_SPARSE_MATRIX_CSR_H
#define ZGEOM_SPARSE_MATRIX_CSR_H

#include <cassert>
#include <algorithm>
#include <vector>
#include "types.h"

namespace ZGeom
{

template<typename T, typename U> class SparseMatrixCSR;

template<typename T, typename U> 
void addCSR(const SparseMatrixCSR<T,U>& m1, const SparseMatrixCSR<T,U>& m2, T beta, SparseMatrixCSR<T,U>& m3); 

template<typename T, typename U>
class SparseMatrixCSR
{
public:
	SparseMatrixCSR() : mVal(NULL), mColIdx(NULL), mRowPtr(NULL), mRowCount(0), mNonzeroCount(0) {}
	SparseMatrixCSR(const std::vector<T>& vVal, const std::vector<U>& vCol, const std::vector<U>& vRow);
	~SparseMatrixCSR();

	void initialize(const std::vector<T>& vVal, const std::vector<U>& vCol, const std::vector<U>& vRow);
	void clear();
	T* nzVal() const { return mVal; }
	U* colIdx() const { return mColIdx; }
	U* rowPtr() const { return mRowPtr; }
	uint rowCount() const { return mRowCount; }
	uint nonzeroCount() const { return mNonzeroCount; }

	// compute CSR addition m3 = m1 + beta*m2
	friend void addCSR<T,U>(const SparseMatrixCSR<T,U>& m1, const SparseMatrixCSR<T,U>& m2, T beta, SparseMatrixCSR<T,U>& m3); 

private:
	T* mVal;    
	U* mColIdx;
	U* mRowPtr;    
	uint mRowCount;
	uint mNonzeroCount;
};


template<typename T, typename U>
inline SparseMatrixCSR<T,U>::SparseMatrixCSR(const std::vector<T>& vVal, const std::vector<U>& vCol, const std::vector<U>& vRow) 
{
	mVal = NULL;
	mColIdx = NULL;
	mRowPtr = NULL;    
	initialize(vVal, vCol, vRow);
}

template<typename T, typename U>
SparseMatrixCSR<T,U>::~SparseMatrixCSR() 
{
	clear();
}

template<typename T, typename U>
inline void SparseMatrixCSR<T,U>::initialize(const std::vector<T>& vVal, const std::vector<U>& vCol, const std::vector<U>& vRow)
{
	clear();

	mRowCount = vRow.size() - 1;
	mNonzeroCount = vRow[mRowCount] - 1;
	mVal = new T[mNonzeroCount];
	mColIdx = new U[mNonzeroCount];
	mRowPtr = new U[mRowCount+1];

	std::copy(vVal.begin(), vVal.end(), mVal);
	std::copy(vCol.begin(), vCol.end(), mColIdx);
	std::copy(vRow.begin(), vRow.end(), mRowPtr);    
}

template<typename T, typename U>
inline void	SparseMatrixCSR<T,U>::clear()
{
	delete []mVal;
	delete []mColIdx;
	delete []mRowPtr;
}

} //end of namespace ZGeom


#ifdef __INTEL_MKL__

typedef ZGeom::SparseMatrixCSR<double, MKL_INT> SparseCSRd;
typedef ZGeom::SparseMatrixCSR<float, MKL_INT>  SparseCSRf;

// compute CSR addition m3 = m1 + beta*m2
//
template<>
inline void ZGeom::addCSR<double,MKL_INT>(const SparseCSRd& m1, const SparseCSRd& m2, double beta, SparseCSRd& m3)
{
	assert(m1.rowCount() == m2.rowCount());    

	MKL_INT rowCount = m1.rowCount();
	MKL_INT nzmax = m1.nonzeroCount() + m2.nonzeroCount();       //maximum count of non-zero elements

	m3.mVal = new double[nzmax];
	m3.mColIdx = new MKL_INT[nzmax];
	m3.mRowPtr = new MKL_INT[rowCount+1]; 
	m3.mRowCount = rowCount;

	char trans = 'N';
	MKL_INT request = 0;
	MKL_INT sort = 3;
	MKL_INT m = rowCount, n = rowCount;
	MKL_INT info;

	mkl_dcsradd(&trans, &request, &sort, &m, &n, m1.nzVal(), m1.colIdx(), m1.rowPtr(), &beta, 
		m2.nzVal(), m2.colIdx(), m2.rowPtr(), m3.mVal, m3.mColIdx, m3.mRowPtr, &nzmax, &info); 

	m3.mNonzeroCount = m3.mRowPtr[rowCount] - 1;
}

template<>
inline void ZGeom::addCSR<float,MKL_INT>(const SparseCSRf& m1, const SparseCSRf& m2, float beta, SparseCSRf& m3)
{
	assert(m1.rowCount() == m2.rowCount());    

	MKL_INT rowCount = m1.rowCount();
	MKL_INT nzmax = m1.nonzeroCount() + m2.nonzeroCount();       //maximum count of non-zero elements

	m3.mVal = new float[nzmax];
	m3.mColIdx = new MKL_INT[nzmax];
	m3.mRowPtr = new MKL_INT[rowCount+1];  
	m3.mRowCount = rowCount; 

	char trans = 'N';
	MKL_INT request = 0;
	MKL_INT sort = 3;
	MKL_INT m = rowCount, n = rowCount;
	MKL_INT info;

	mkl_scsradd(&trans, &request, &sort, &m, &n, m1.nzVal(), m1.colIdx(), m1.rowPtr(), &beta, 
		m2.nzVal(), m2.colIdx(), m2.rowPtr(), m3.mVal, m3.mColIdx, m3.mRowPtr, &nzmax, &info);    

	m3.mNonzeroCount = m3.mRowPtr[rowCount] - 1;
}

#endif //__INTEL_MKL__


#endif