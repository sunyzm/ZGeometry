#ifndef ZMATH_SPARSE_MATRIX_H
#define ZMATH_SPARSE_MATRIX_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

namespace ZMath
{

template<typename T> class Laplacian;

template<typename T>
struct TupleCOO
{
    TupleCOO() : i(0), j(0), v(0.) {}
    TupleCOO(int ii, int jj, T vv) : i(ii), j(jj), v(vv){}
    
    unsigned int i, j;
    T v;
    template<typename U> friend bool operator < (const TupleCOO<U>& t1, const TupleCOO<U>& t2);
};

template<typename U>
bool operator < (const TupleCOO<U>& t1, const TupleCOO<U>& t2) {
    return t1.i < t2.i || (t1.i == t2.i && t1.j < t2.j);
}


template<typename T>
class SparseMatrix
{
public:
    friend class Laplacian<T>;    
    enum MatrixForm {MAT_UPPER, MAT_LOWER, MAT_FULL};
	SparseMatrix();

    unsigned rowCount() const { return mRowCount; }
    unsigned colCount() const { return mColCount; }
    unsigned nonzeroCount() const { return mNonzeroCount; }

    const std::vector<TupleCOO<T> >& getCOO_const() const { return mCOO; } 
	//std::vector<TupleCOO<T> >& getCOO { return mCOO; }
	template<typename U> void loadCOO(unsigned rowCount, unsigned colCount, const std::vector<U>& vRowIdx, const std::vector<U>& vColIdx, const std::vector<T>& vVal);
    template<typename U> void loadCSR(unsigned rowCount, unsigned colCount, T val[], U colIdx[], U rowPtr[]);
    template<typename U> void getCSR(std::vector<T>& val, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form = MAT_UPPER) const;

    void print(std::ostream& out) const;
    bool testNoEmptyRow() const;
	bool testSymmetric(T eps = 1e-7) const;
    void truncate(MatrixForm form = MAT_UPPER);
    void symmetrize();
    void fillEmptyDiagonal();       // fill empty diagonal element with 0
    
private:
    std::vector<TupleCOO<T>> mCOO;
    unsigned int mRowCount;
    unsigned int mColCount;
    unsigned int mNonzeroCount;	 
	bool		 isSingle;
};

template<typename T>
ZMath::SparseMatrix<T>::SparseMatrix()
{
	if (sizeof(T) == sizeof(float)) this->isSingle = true;
	else this->isSingle = false;
}



template<typename T>
template<typename U>
void SparseMatrix<T>::loadCOO( unsigned rowCount, unsigned colCount, const std::vector<U>& vRowIdx, const std::vector<U>& vColIdx, const std::vector<T>& vVal )
{
	assert(vRowIdx.size() == vColIdx.size() && vColIdx.size() == vVal.size());
	int nnz = vVal.size();

	this->mRowCount = rowCount; 
	this->mColCount = colCount;

	mCOO.resize(nnz);
	for (int k = 0; k < nnz; ++k) {
		mCOO[k].i = vRowIdx[k];
		mCOO[k].j = vColIdx[k];
		mCOO[k].v = vVal[k];
	}
}


template<typename T> void
SparseMatrix<T>::print(std::ostream& out) const
{
    for (typename std::vector<TupleCOO<T>>::const_iterator iter = mCOO.begin(); iter != mCOO.end(); ++iter) {
        out << iter->i << ' ' << iter->j << ' ' << iter->v << std::endl;
    }    
}

template<typename T> 
bool SparseMatrix<T>::testNoEmptyRow() const
{
    std::vector<bool> rowIsEmpty(mRowCount, true);
    for (typename std::vector<TupleCOO<T>>::const_iterator iter = mCOO.begin(); iter != mCOO.end(); ++iter) {
        rowIsEmpty[iter->i - 1] = false;
    }
    for (uint k = 0; k < mRowCount; ++k) {
        if (rowIsEmpty[k]) return false;
    }
    return true;
}

template<typename T>
bool SparseMatrix<T>::testSymmetric(T eps /*= 1e-7*/) const
{
    for (unsigned k = 0; k < mNonzeroCount; ++k) {
        if (mCOO[k].i >= mCOO[k].j) continue;
        bool symElemFound = false;
        for (unsigned l = 0; l < mNonzeroCount; ++l) {
            if (mCOO[l].i != mCOO[k].j || mCOO[l].j != mCOO[k].i) continue;
            if (abs(mCOO[l].v - mCOO[k].v) < eps) symElemFound = true;
            break;
        }
        if (!symElemFound) return false;
    }
    return true;
}

template<typename T> 
void SparseMatrix<T>::symmetrize() 
{
    for (unsigned k = 0; k < mNonzeroCount; ++k) {
        if (mCOO[k].i != mCOO[k].j) mCOO.push_back(TupleCOO<T>(mCOO[k].j, mCOO[k].i, mCOO[k].v));
    }
    std::sort(mCOO.begin(), mCOO.end(), std::less<TupleCOO<T>>());
    mNonzeroCount = mCOO.size();
}

template<typename T> 
void SparseMatrix<T>::fillEmptyDiagonal()
{
    assert(mRowCount == mColCount);

    std::vector<bool> emptyDiag(mRowCount, true);
    for (std::vector<TupleCOO<T>>::iterator iter = mCOO.begin(); iter != mCOO.end(); ++iter) {
        if (iter->i == iter->j) emptyDiag[iter->i - 1] = false;
    }
    
    for (uint k = 0; k < mRowCount; ++k) {
        if (emptyDiag[k]) mCOO.push_back(TupleCOO<T>(k+1,k+1,0.0)); 
    }
    
    mNonzeroCount = mCOO.size();
}

template<typename T> 
void SparseMatrix<T>::truncate(MatrixForm form /*=MAT_UPPER*/)
{
    assert(mRowCount == mColCount);
    
    if (form == MAT_UPPER) {
        for (typename std::vector<TupleCOO<T>>::iterator iter = mCOO.begin(); iter != mCOO.end(); ) {
            if (iter->i > iter->j) iter = mCOO.erase(iter);
            else iter++;
        }
    }
    else if (form == MAT_LOWER) {
        for (typename std::vector<TupleCOO<T>>::iterator iter = mCOO.begin(); iter != mCOO.end();) {
            if (iter->i < iter->j) iter = mCOO.erase(iter);
            else iter++;
        }
    }
    
    mNonzeroCount = mCOO.size();
}

template<typename T>
template<typename U>
void SparseMatrix<T>::getCSR(std::vector<T>& val, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form /*=MAT_UPPER*/) const 
{
    if (mRowCount != mColCount) {
        assert(form == MAT_FULL);       // only square matrix has upper or lower CSR
    }
    assert (testNoEmptyRow());          // CSR format requires at least one element in each row

    val.clear();
    colIdx.clear();
    rowPtr.clear();

    std::sort(mCOO.begin(), mCOO.end(), std::less<TupleCOO<T>>());
    int prevRow = 0;
    typename std::vector<TupleCOO<T>>::const_iterator iter = mCOO.begin();
    for (; iter != mCOO.end(); ++iter) {
        if (form == MAT_UPPER && iter->i > iter->j) continue;   // only upper triangle elements
        if (form == MAT_LOWER && iter->i < iter->j) continue;   // only lower triangle elements

        uint curRow = iter->i;
        val.push_back(iter->v);
        colIdx.push_back(iter->j);
        if (1 == curRow - prevRow) rowPtr.push_back(val.size()); // new row
        prevRow = curRow;               
    }
    rowPtr.push_back(val.size()+1);     // trailing element of rowPtr stores NNZ + 1
}

template<typename T>
template<typename U> 
void SparseMatrix<T>::loadCSR(unsigned rowCount, unsigned colCount, T val[], U colIdx[], U rowPtr[])
{
    mRowCount = rowCount;
    mColCount = colCount;
    mNonzeroCount = rowPtr[rowCount] - 1;

    mCOO.resize(mNonzeroCount);
    for (unsigned r = 1; r <= mRowCount; ++r) {
        for (unsigned k = (unsigned)rowPtr[r-1]; k < (unsigned)rowPtr[r]; ++k) {
            mCOO[k-1].i = r;
            mCOO[k-1].j = colIdx[k-1];
            mCOO[k-1].v = val[k-1];
        }
    }
}

typedef SparseMatrix<double> SparseMatrixD;
typedef SparseMatrix<float>	 SparseMatrixS;

} //end of namespace ZMath

#endif
