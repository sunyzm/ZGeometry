#ifndef ZGEOM_SPARSE_MATRIX_H
#define ZGEOM_SPARSE_MATRIX_H
#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <algorithm>
#include "common.h"
#include "SparseMatrixCSR.h"

namespace ZGeom {

template<typename T> class MatElem;
template<typename T> class VecN;
template<typename T> class SparseMatrix;    
template<typename T> VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

enum MatrixForm {MAT_UPPER, MAT_LOWER, MAT_FULL};

template<typename T>
class MatElem
{
public:
	friend class SparseMatrix<T>;

	MatElem() : mRow(0), mCol(0), mVal(0.) {}
	MatElem(uint ii, uint jj, T vv) : mRow(ii), mCol(jj), mVal(vv) {}

	uint row() const { return mRow; }
	uint col() const { return mCol; }
	T    val() const { return mVal; }
	T&   val() { return mVal; }
    bool operator< (const MatElem<T>& t2) const
    {
        return this->mRow < t2.mRow || (this->mRow == t2.mRow && this->mCol < t2.mCol);
    }

public:   
	uint mRow, mCol;
	T mVal;
};


template<typename T>
class SparseMatrix
{
public:
	SparseMatrix() : mRowCount(0), mColCount(0), mNonzeroCount(0) {}
	SparseMatrix(int row, int col) : mRowCount(row), mColCount(col), mNonzeroCount(0) {}
    SparseMatrix<T>& operator=(const SparseMatrix<T>& mat2) = default;
    SparseMatrix(const SparseMatrix<T>& mat2) = default;
    SparseMatrix<T>& operator=(SparseMatrix<T>&& mat2)
    {
        this->mColCount = mat2.mColCount;
        this->mRowCount = mat2.mRowCount;
        this->mNonzeroCount = mat2.mNonzeroCount;
        this->mElements = std::move(mat2.mElements);
        return *this;
    }
    SparseMatrix(SparseMatrix<T>&& mat2) { *this = std::move(mat2); }

    /* basic attributes and elements access */
	uint rowCount() const { return mRowCount; }
	uint colCount() const { return mColCount; }
	uint nonzeroCount() const { return mElements.size(); }
	bool empty() const { return mRowCount == 0 || mColCount == 0; }
    T getElemValByIndex(uint index) const;
    T& getElemValByIndex(uint index);
    T getElemVal(uint row, uint col) const;
    T& getElemVal(uint row, uint col);
    MatElem<T>& getElemByIndex(uint index) { return mElements[index]; }
    const MatElem<T>& getElemByIndex(uint index) const { return mElements[index]; }
    const std::vector<MatElem<T> >& allElements() const { return mElements; }
    std::vector<MatElem<T> >& allElements() { return mElements; }
    T operator() (uint row, uint col) const;
    T& operator() (uint row, uint col);

    /* modify elements */
    void resize(int row, int col);
    void clear() { resize(0, 0); }
	void insertElem(uint row, uint col, T val);
	void removeElem(uint row, uint col);   
	void copyElements(const SparseMatrix<T>& mat2);
	void computeSubMatrix(const std::vector<int>& vSelected, SparseMatrix<T>& subMat) const;
    const SparseMatrix<T>& scale(T scalar);
    const SparseMatrix<T>& operator *= (T coeff);
    void truncate(MatrixForm form); // truncate matrix to upper or lower triangle matrix
    void symmetrize();              // turn upper or lower matrix into symmetric one
    void fillEmptyDiagonal();       // fill empty diagonal elements with 0
    void setToIdentity(uint order); // make the matrix an identity matrix
    void sortElements();

    T frobeniusNorm() const;
    bool testNoEmptyRow() const;
    bool testSymmetric(double eps = 1e-7) const;
    void print(std::ostream& out) const;
    void print(const std::string& path) const;
    void read(std::istream& in);

    std::vector<T> getDiagonal() const;

    template<typename F>
	void convertFromDiagonal(const std::vector<F>& diag);    

	template<typename U, typename F>
	void convertFromCOO(uint rowCount, uint colCount, const std::vector<U>& rowInd, std::vector<U>& colInd, const std::vector<F>& val);

	template<typename U, typename F>
	void convertFromCOO(uint rowCount, uint colCount, const std::vector<std::tuple<U,U,F> >& vElem);

	template<typename U, typename F>
	void convertToCOO(std::vector<U>& rowInd, std::vector<U>& colInd, std::vector<F>& val, MatrixForm form) const;

	template<typename U, typename F>
	void convertToCOO(U* rowInd, U* colInd, F* val, int& nonzeroCount, MatrixForm form) const;    

	template<typename U, typename F> 
	void convertFromCSR(uint rowCount, uint colCount, F nzVal[], U colIdx[], U rowPtr[]);

	template<typename U, typename F> 
	void convertToCSR(std::vector<F>& nzVal, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form) const;

	template<typename U, typename F> 
	void convertToCSR(SparseMatrixCSR<F,U>& MatCSR, MatrixForm form) const;

	template<typename F> 
	void convertFromFull(F* fullMat, double sparse_eps = 1e-10);

	template<typename F> 
	void convertToFull(F* fullMat, MatrixForm form) const;

    friend VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

private:
	std::vector<MatElem<T> > mElements;
	uint mRowCount;
	uint mColCount;
	uint mNonzeroCount;
};

typedef SparseMatrix<double> SparseMatrixd;

template<typename T>
inline void SparseMatrix<T>::read(std::istream& in)
{
    in >> mRowCount >> mColCount;
    in >> mNonzeroCount;
    mElements.resize(mNonzeroCount);
    for (uint k = 0; k < mNonzeroCount; ++k) {
        in >> mElements[k].mRow >> mElements[k].mCol >> mElements[k].mVal;
    }
}

template<typename T>
inline void SparseMatrix<T>::print(std::ostream& out) const
{
    out << mRowCount << ' ' << mColCount << std::endl;
    out << mNonzeroCount << std::endl;
    for (auto iter = mElements.begin(); iter != mElements.end(); ++iter) {
        out << iter->mRow << ' ' << iter->mCol << ' ' << iter->mVal << std::endl;
    }
}

template<typename T>
inline void SparseMatrix<T>::print(const std::string& path) const
{
    std::ofstream ofs(path);
    print(ofs);
    ofs.close();
}

template<typename T>
void SparseMatrix<T>::resize(int row, int col)
{
    mRowCount = row;
    mColCount = col;
    mNonzeroCount = 0;
    mElements.clear();
}

template<typename T>
inline T SparseMatrix<T>::getElemValByIndex(uint index) const
{
    assert(index < mNonzeroCount);
    return mElements[index].mVal;
}

template<typename T>
inline T& SparseMatrix<T>::getElemValByIndex(uint index)
{
    assert(index < mNonzeroCount);
    return mElements[index].mVal;
}

template<typename T>
inline T SparseMatrix<T>::getElemVal(uint row, uint col) const
{
    assert(row <= mRowCount && col <= mColCount);

    for (uint k = 0; k < mElements.size(); ++k) {
        if (mElements[k].mRow == row && mElements[k].mCol == col)
            return mElements[k].mVal;
    }
    return 0;
}

template<typename T>
inline T& SparseMatrix<T>::getElemVal(uint row, uint col)
{
    assert(row <= mRowCount && col <= mColCount);

    for (uint k = 0; k < mElements.size(); ++k) {
        if (mElements[k].mRow == row && mElements[k].mCol == col)
            return mElements[k].mVal;
    }

    mElements.push_back(MatElem<T>(row, col, 0.0));
    return mElements.back().mVal;
}

template<typename T>
inline T& SparseMatrix<T>::operator() (uint row, uint col)
{
    return getElemVal(row, col);
}

template<typename T>
inline T SparseMatrix<T>::operator() (uint row, uint col) const
{
    return getElemVal(row, col);
}

template<typename T>
inline void SparseMatrix<T>::insertElem(uint row, uint col, T val)
{
    mElements.push_back(MatElem<T>(row, col, val));
    mNonzeroCount = (int)mElements.size();
}

template<typename T>
inline void SparseMatrix<T>::removeElem(uint row, uint col)
{
    for (auto iter = mElements.begin(); iter != mElements.end();) {
        if (iter->row() == row && iter->col() == col) {
            iter = mElements.erase(iter);
            continue;
        }
        ++iter;
    }
    mNonzeroCount = mElements.size();
}

template<typename T>
inline void SparseMatrix<T>::copyElements(const SparseMatrix<T>& mat2)
{
    for (auto elem : mat2.mElements) {
        if (elem.row() > mRowCount || elem.col() > mColCount)
            throw std::runtime_error("SparseMatrix::copyElements encounter invalid element!");
    }
    mElements = mat2.mElements;
    mNonzeroCount = mElements.size();
}

template<typename T>
T SparseMatrix<T>::frobeniusNorm() const
{
    T sum(0);
    for (auto e : mElements) sum += e.val() * e.val();
    return std::sqrt(sum);
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
inline bool SparseMatrix<T>::testSymmetric(double eps/*=1e-7*/) const
{
    for (uint k = 0; k < mNonzeroCount; ++k) {
        if (mElements[k].row() >= mElements[k].col()) continue;
        bool symElemFound = false;
        for (uint l = 0; l < mNonzeroCount; ++l) {
            if (mElements[k].row() == mElements[l].col() && mElements[k].col() == mElements[l].row()) { // element in symmetric position found
                if (fabs(mElements[l].val() - mElements[k].val()) < eps)
                    symElemFound = true;
                break;
            }
        }
        if (!symElemFound) return false;
    }
    return true;
}

template<typename T>
inline void SparseMatrix<T>::truncate(MatrixForm form /*=MAT_UPPER*/)
{
    assert(mRowCount == mColCount);

    if (form == MAT_UPPER) {
        for (typename std::vector< MatElem<T> >::iterator iter = mElements.begin(); iter != mElements.end();) {
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
const SparseMatrix<T>& SparseMatrix<T>::scale(T scalar)
{
    for (MatElem<T>& elem : mElements) elem.mVal *= scalar;
    return *this;
}

template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator *= (T coeff)
{
    for (MatElem<T>& e : mElements) e.mVal *= coeff;
    return *this;
}

template<typename T>
inline void SparseMatrix<T>::setToIdentity(uint order)
{
    mRowCount = mColCount = mNonzeroCount = order;
    mElements.resize(order);
    for (uint k = 0; k < order; ++k)
        mElements[k] = MatElem<T>(k + 1, k + 1, 1.0);
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
        if (emptyDiag[k]) mElements.push_back(MatElem<T>(k + 1, k + 1, 0));
    }

    mNonzeroCount = mElements.size();
}

template<typename T>
inline void SparseMatrix<T>::symmetrize()
{
    for (uint k = 0; k < mNonzeroCount; ++k) {
        if (mElements[k].row() != mElements[k].col()) mElements.push_back(MatElem<T>(mElements[k].col(), mElements[k].row(), mElements[k].val()));
    }
    std::sort(mElements.begin(), mElements.end(), std::less<MatElem<T> >());
    mNonzeroCount = mElements.size();
}

template<typename T>
void SparseMatrix<T>::sortElements()
{
    std::sort(mElements.begin(), mElements.end(), std::less<MatElem<T>>());
}

template<typename T>
inline void SparseMatrix<T>::computeSubMatrix(const std::vector<int>& vSelected, SparseMatrix<T>& subMat) const
{
    subMat.mRowCount = subMat.mColCount = (int)vSelected.size();
    subMat.mElements.clear();

    std::unordered_map<int, int> subIndexMap;
    for (size_t i = 0; i < vSelected.size(); ++i) {
        subIndexMap.insert(std::make_pair(vSelected[i], (int)i));
    }

    for (auto elem : this->mElements) {
        int rowIdx = elem.mRow - 1, colIdx = elem.mCol - 1;
        if (subIndexMap.find(rowIdx) != subIndexMap.end() && subIndexMap.find(colIdx) != subIndexMap.end()) {
            subMat.insertElem(subIndexMap[rowIdx] + 1, subIndexMap[colIdx] + 1, elem.mVal);
        }
    }

    subMat.mNonzeroCount = (int)subMat.mElements.size();
}

template<typename T>
std::vector<T> SparseMatrix<T>::getDiagonal() const
{
    assert(mRowCount == mColCount);
    std::vector<T> diag(mRowCount, (T)0);
    for (auto iter = mElements.begin(); iter != mElements.end(); ++iter) {
        if (iter->row() == iter->col())
            diag[iter->row() - 1] = iter->val();
    }
    return diag;
}

template<typename T>
template<typename F>
void SparseMatrix<T>::convertFromDiagonal(const std::vector<F>& diag)
{
    mRowCount = mColCount = mNonzeroCount = (uint)diag.size();
    mElements.resize(mRowCount);
    for (uint k = 0; k < mRowCount; ++k) {
        mElements[k].mRow = mElements[k].mCol = k + 1;
        mElements[k].mVal = diag[k];
    }
}

template<typename T>
template<typename U, typename F>
void SparseMatrix<T>::convertFromCOO(uint rowCount, uint colCount, const std::vector<U>& rowInd, std::vector<U>& colInd, const std::vector<F>& val)
{
    assert(rowInd.size() == colInd.size() && rowInd.size() == val.size());

    mRowCount = rowCount;
    mColCount = colCount;
    mNonzeroCount = (uint)val.size();

    mElements.clear();
    mElements.reserve(mNonzeroCount);
    for (uint k = 0; k < mNonzeroCount; ++k) {
        assert(rowInd[k] > 0 && colInd[k] > 0);
        assert(rowInd[k] <= (U)rowCount && colInd[k] <= (U)colCount);
        mElements.push_back(MatElem<T>(rowInd[k], colInd[k], val[k]));
    }
}

template<typename T>
template<typename U, typename F>
void SparseMatrix<T>::convertFromCOO(uint rowCount, uint colCount, const std::vector< std::tuple<U, U, F> >& vElem)
{
    int nnz = (int)vElem.size();
    std::vector<U> rowInd, colInd;
    std::vector<F> val;
    for (auto& elem : vElem) {
        rowInd.push_back(std::get<0>(elem));
        colInd.push_back(std::get<1>(elem));
        val.push_back(std::get<2>(elem));
    }

    convertFromCOO(rowCount, colCount, rowInd, colInd, val);
}

template<typename T>
template<typename U, typename F>
inline void SparseMatrix<T>::convertToCOO(std::vector<U>& rowInd, std::vector<U>& colInd, std::vector<F>& val, MatrixForm form) const
{
    assert(mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR

    rowInd.clear();
    colInd.clear();
    val.clear();

    typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin();
    for (; iter != mElements.end(); ++iter) {
        if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only upper triangle elements
        if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only lower triangle elements

        rowInd.push_back(static_cast<U>(iter->row()));
        colInd.push_back(static_cast<U>(iter->col()));
        val.push_back(iter->val());
    }
}

template<typename T>
template<typename U, typename F>
inline void SparseMatrix<T>::convertToCOO(U* rowInd, U* colInd, F* val, int& nonzeroCount, MatrixForm form) const
{
    assert(mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR

    nonzeroCount = 0;
    typename std::vector< MatElem<T> >::const_iterator iter = sortedElements.begin();
    for (; iter != sortedElements.end(); ++iter) {
        if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only preserve upper triangle elements
        if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only preserve lower triangle elements

        rowInd[nonzeroCount] = static_cast<U>(iter->row());
        colInd[nonzeroCount] = static_cast<U>(iter->col());
        val[nonzeroCount] = iter->val();

        nonzeroCount++;
    }
}

template<typename T>
template<typename U, typename F>
inline void SparseMatrix<T>::convertToCSR(std::vector<F>& nzVal, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form /*=MAT_UPPER*/) const
{
    assert(mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR        
    assert(testNoEmptyRow());          // CSR format requires at least one element in each row

    nzVal.clear();
    colIdx.clear();
    rowPtr.clear();

    std::vector< MatElem<T> > sortedElements = mElements;
    std::sort(sortedElements.begin(), sortedElements.end(), std::less< MatElem<T> >());
    uint prevRow = 0;
    typename std::vector< MatElem<T> >::const_iterator iter = sortedElements.begin();
    for (; iter != sortedElements.end(); ++iter) {
        if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only upper triangle elements
        if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only lower triangle elements

        uint curRow = iter->row();
        nzVal.push_back(iter->val());
        colIdx.push_back(iter->col());
        if (1 == curRow - prevRow) rowPtr.push_back((int)nzVal.size()); // new row
        prevRow = curRow;
    }
    rowPtr.push_back((int)nzVal.size() + 1);     // trailing element of rowPtr stores NNZ + 1
}

template<typename T>
template<typename U, typename F>
inline void SparseMatrix<T>::convertToCSR(SparseMatrixCSR<F, U>& matCSR, MatrixForm form) const
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
    assert(mRowCount == mColCount || form == MAT_FULL);  // only square matrix has upper or lower CSR   
    for (uint k = 0; k < mRowCount * mColCount; ++k) fullMat[k] = 0.0;

    typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin();
    for (; iter != mElements.end(); ++iter) {
        if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only upper triangle elements
        if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only lower triangle elements

        fullMat[(iter->row() - 1) * mColCount + (iter->col() - 1)] = iter->val();
    }
}

template<typename T>
template<typename U, typename F>
inline void SparseMatrix<T>::convertFromCSR(uint rowCount, uint colCount, F nzVal[], U colIdx[], U rowPtr[])
{
    mRowCount = rowCount;
    mColCount = colCount;
    mNonzeroCount = rowPtr[rowCount] - 1;

    mElements.resize(mNonzeroCount);
    for (uint r = 1; r <= mRowCount; ++r) {
        for (uint k = (uint)rowPtr[r - 1]; k < (uint)rowPtr[r]; ++k) {
            mElements[k - 1].mRow = r;
            mElements[k - 1].mCol = colIdx[k - 1];
            mElements[k - 1].mVal = nzVal[k - 1];
        }
    }
}

}   // end of namespace ZGeom

#endif
