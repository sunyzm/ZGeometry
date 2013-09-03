#ifndef ZGEOM_SPARSE_MATRIX_CSR_H
#define ZGEOM_SPARSE_MATRIX_CSR_H

#include <cassert>
#include <algorithm>
#include <vector>
#include "common.h"

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

// compute CSR addition m3 = m1 + beta*m2
//
template<typename T, typename U>
void addCSR(const SparseMatrixCSR<T, U>& m1, const SparseMatrixCSR<T, U>& m2, T beta, SparseMatrixCSR<T, U>& m3);

} //end of namespace ZGeom




#endif