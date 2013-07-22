#ifndef ZGEOM_SPARSE_MATRIX_H
#define ZGEOM_SPARSE_MATRIX_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "types.h"
#include "SparseMatrixCSR.h"


namespace ZGeom 
{

template<typename T> class MatElem;
template<typename T> class SparseMatrix;
template<typename T> class Laplacian;
template<typename T> bool operator< (const MatElem<T>& t1, const MatElem<T>& t2);

template<typename T>
class MatElem
{
public:
	friend class SparseMatrix<T>;
	friend bool operator< <T>(const MatElem<T>& t1, const MatElem<T>& t2);

	MatElem() : mRow(0), mCol(0), mVal(0.) {}
	MatElem(uint ii, uint jj, T vv) : mRow(ii), mCol(jj), mVal(vv) {}
	uint row() const { return mRow; }
	uint col() const { return mCol; }
	T val() const { return mVal; }

private:   
	uint mRow, mCol;
	T mVal;
};

template<typename T> 
inline bool operator < (const MatElem<T>& t1, const MatElem<T>& t2) {
	return t1.mRow < t2.mRow || (t1.mRow == t2.mRow && t1.mCol < t2.mCol);
}


template<typename T>
class SparseMatrix
{
public:
	friend class Laplacian<T>;    
	enum MatrixForm {MAT_UPPER, MAT_LOWER, MAT_FULL};

	uint rowCount() const;
	uint colCount() const;
	uint nonzeroCount() const;
	const MatElem<T>& getElem(uint index) const;
	T getElemVal(uint index) const;

	template<typename U> 
	void convertFromCSR(uint rowCount, uint colCount, T nzVal[], U colIdx[], U rowPtr[]);

	template<typename U> 
	void convertToCSR(std::vector<T>& nzVal, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form = MAT_UPPER) const;

	template<typename U> 
	void convertToCSR(SparseMatrixCSR<T,U>& MatCSR, MatrixForm form = MAT_UPPER) const;

	template<typename F> 
	void convertFromFull(F* fullMat, double sparse_eps = 1e-10);

	template<typename F> 
	void convertToFull(F* fullMat, MatrixForm form = MAT_UPPER) const;

	void truncate(MatrixForm form = MAT_UPPER); // truncate matrix 
	void symmetrize();              // turn upper or lower matrix into symmetric one
	void fillEmptyDiagonal();       // fill empty digonal elements with 0
	void toIdentity(uint order);     // turn into an identity matrix
	void print(std::ostream& out) const;

private:
	bool testNoEmptyRow() const;
	bool testSymmetric() const;

	std::vector< MatElem<T> > mElements;
	uint mRowCount;
	uint mColCount;
	uint mNonzeroCount;
};

typedef SparseMatrix<double> SparseMatrixD;
typedef SparseMatrix<double> SparseMatrixF;


template<typename T>
inline uint	SparseMatrix<T>::rowCount() const
{
	return mRowCount;    
}

template<typename T>
inline uint	SparseMatrix<T>::colCount() const
{
	return mColCount;    
}

template<typename T>
inline uint	SparseMatrix<T>::nonzeroCount() const
{
	return mNonzeroCount;
}

template<typename T>
inline const MatElem<T>& SparseMatrix<T>::getElem(uint index) const
{   
	return mElements[index];
}

template<typename T>
inline T SparseMatrix<T>::getElemVal(uint index) const
{
	return mElements[index].val();
}

template<typename T>
inline void SparseMatrix<T>::toIdentity(uint order)
{
	mRowCount = mColCount = mNonzeroCount = order;
	mElements.clear();
	for (uint k = 0; k < order; ++k) {
		mElements.push_back(MatElem<T>(k+1, k+1, 1.0));   
	}    
}

template<typename T> 
inline void	SparseMatrix<T>::print(std::ostream& out) const
{
	for (typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin(); iter != mElements.end(); ++iter) {
		out << iter->i << ' ' << iter->j << ' ' << iter->v << std::endl;
	}    
}

template<typename T> 
inline bool SparseMatrix<T>::testNoEmptyRow() const
{
	std::vector<bool> rowIsEmpty(mRowCount, true);
	for (typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin(); iter != mElements.end(); ++iter) {
		rowIsEmpty[iter->mRow - 1] = false;
	}
	for (uint k = 0; k < mRowCount; ++k) {
		if (rowIsEmpty[k]) return false;
	}
	return true;
}

template<typename T>
inline bool SparseMatrix<T>::testSymmetric() const
{
	assert(mRowCount == mColCount);

	for (uint k = 0; k < mNonzeroCount; ++k) {
		if (mElements[k].row() >= mElements[k].col()) continue;
		bool symElemFound = false;
		for (uint l = 0; l < mNonzeroCount; ++l) {
			if (mElements[k].row() == mElements[l].col() || mElements[k].col() == mElements[l].row()) { // element in symmetric position found
				if (fabs(mElements[l].val() - mElements[k].val()) < 1e-7) symElemFound = true;
				break;
			}
		}
		if (!symElemFound) return false;
	}
	return true;
}

template<typename T> 
inline void SparseMatrix<T>::symmetrize()   
{
	assert(mRowCount == mColCount);

	for (uint k = 0; k < mNonzeroCount; ++k) {
		if (mElements[k].row() != mElements[k].col()) mElements.push_back(MatElem<T>(mElements[k].col(), mElements[k].row(), mElements[k].val()));
	}
	std::sort(mElements.begin(), mElements.end(), std::less<MatElem<T> >());
	mNonzeroCount = mElements.size();
}

template<typename T> 
inline void SparseMatrix<T>::fillEmptyDiagonal()
{
	assert(mRowCount == mColCount);

	std::vector<bool> emptyDiag(mRowCount, true);
	for (typename std::vector< MatElem<T> >::iterator iter = mElements.begin(); iter != mElements.end(); ++iter) {
		if (iter->row() == iter->col()) emptyDiag[iter->row() - 1] = false;
	}

	for (uint k = 0; k < mRowCount; ++k) {
		if (emptyDiag[k]) mElements.push_back(MatElem<T>(k+1,k+1,0.0)); 
	}

	mNonzeroCount = mElements.size();
}

template<typename T> 
inline void SparseMatrix<T>::truncate(MatrixForm form /*=MAT_UPPER*/)
{
	assert(mRowCount == mColCount);

	if (form == MAT_UPPER) {
		for (typename std::vector< MatElem<T> >::iterator iter = mElements.begin(); iter != mElements.end(); ) {
			if (iter->row() > iter->col()) iter = mElements.erase(iter);
			else iter++;
		}
	}
	else if (form == MAT_LOWER) {
		for (typename std::vector< MatElem<T> >::iterator iter = mElements.begin(); iter != mElements.end();) {
			if (iter->row() < iter->col()) iter = mElements.erase(iter);
			else iter++;
		}
	}

	mNonzeroCount = mElements.size();
}

template<typename T>
template<typename U>
inline void SparseMatrix<T>::convertToCSR(std::vector<T>& nzVal, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form /*=MAT_UPPER*/) const
{
	assert (mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR        
	assert (testNoEmptyRow());          // CSR format requires at least one element in each row

	nzVal.clear();
	colIdx.clear();
	rowPtr.clear();

	std::vector< MatElem<T> > sortedElements = mElements;
	std::sort(sortedElements.begin(), sortedElements.end(), std::less< MatElem<T> >());
	int prevRow = 0;
	typename std::vector< MatElem<T> >::const_iterator iter = sortedElements.begin();
	for (; iter != sortedElements.end(); ++iter) {
		if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only upper triangle elements
		if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only lower triangle elements

		uint curRow = iter->row();
		nzVal.push_back(iter->val());
		colIdx.push_back(iter->col());
		if (1 == curRow - prevRow) rowPtr.push_back(nzVal.size()); // new row
		prevRow = curRow;               
	}
	rowPtr.push_back(nzVal.size() + 1);     // trailing element of rowPtr stores NNZ + 1
}

template<typename T>
template<typename U>
inline void SparseMatrix<T>::convertToCSR(SparseMatrixCSR<T,U>& matCSR, MatrixForm form) const
{
	std::vector<T> nzVal;
	std::vector<U> colIdx;
	std::vector<U> rowPtr;

	convertToCSR(nzVal, colIdx, rowPtr, form);
	matCSR.initialize(nzVal, colIdx, rowPtr);
}

template<typename T>
template<typename F>
inline void SparseMatrix<T>::convertToFull(F* fullMat, MatrixForm form /*= MAT_UPPER*/) const
{
	assert (mRowCount == mColCount || form == MAT_FULL);  // only square matrix has upper or lower CSR   
	for (uint k = 0; k < mRowCount * mColCount; ++k) fullMat[k] = 0.0;

	typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin();
	for (; iter != mElements.end(); ++iter) {
		if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only upper triangle elements
		if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only lower triangle elements

		fullMat[iter->row()*mColCount + iter->col()] = iter->val();    
	}
}

template<typename T>
template<typename U> 
inline void SparseMatrix<T>::convertFromCSR(uint rowCount, uint colCount, T nzVal[], U colIdx[], U rowPtr[])
{
	mRowCount = rowCount;
	mColCount = colCount;
	mNonzeroCount = rowPtr[rowCount] - 1;

	mElements.resize(mNonzeroCount);
	for (uint r = 1; r <= mRowCount; ++r) {
		for (uint k = (uint)rowPtr[r-1]; k < (uint)rowPtr[r]; ++k) {
			mElements[k-1].mRow = r;
			mElements[k-1].mCol = colIdx[k-1];
			mElements[k-1].mVal = nzVal[k-1];
		}
	}
}

} //end of namespace ZGeom

#endif
