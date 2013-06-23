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

    uint rowCount() const { return mRowCount; }
    uint colCount() const { return mColCount; }
    uint nonzeroCount() const { return mNonzeroCount; }

    const std::vector<TupleCOO<T>>& getCOO() const { return mCOO; } 
    template<typename U> void loadCSR(unsigned rowCount, unsigned colCount, T val[], U colIdx[], U rowPtr[]);
    template<typename U> void getCSR(std::vector<T>& val, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form = MAT_UPPER) const;

    void print(std::ostream& out) const;
    bool testNoEmptyRow() const;
	bool testSymmetric(double eps = 1e-7) const;
    void truncate(MatrixForm form = MAT_UPPER);
    void symmetrize();
    void fillEmptyDiagonal();       // fill empty diagonal element with 0
    
private:
    std::vector<TupleCOO<T>> mCOO;
    uint mRowCount;
    uint mColCount;
    uint mNonzeroCount;
};

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
bool SparseMatrix<T>::testSymmetric(double eps /*= 1e-7*/) const
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







template<typename T>
class SparseCSR;

template<typename T>
class SparseCOO
{
public:
    uint mRowNum, mColNum;
    uint mNNZ;  // number of non-zero elements
    
    std::vector<T> mVal;
    std::vector<uint> mRowInd, mColInd;
    
    int constructFromDiag(const std::vector<T>& vWW); 
    bool judgeSymmetric() const;
    int toPackedSymmetricU(T ap[]) const;
    int toCSR(SparseCSR<T>& csr) const;

    void symmetrize();

    template<typename U> friend std::ostream& operator<< (std::ostream& out, const SparseCOO<U> &coo);
};

template<typename U> 
std::ostream& operator<< (std::ostream& out, const SparseCOO<U> &coo)
{
    out << "values\t=(";
    for (std::vector<U>::const_iterator iter = coo.mVal.begin(); iter != coo.mVal.end(); ++iter) 
        out << *iter << ' ';
    out << ")\nRowIdx\t=(";
    for (std::vector<uint>::const_iterator iter = coo.mRowInd.begin(); iter != coo.mRowInd.end(); ++iter)
        out << *iter << ' ';
    out << ")\nColIdx\t=(";
    for (std::vector<uint>::const_iterator iter = coo.mColInd.begin(); iter != coo.mColInd.end(); ++iter)
        out << *iter << ' ';
    out << ")" << std::endl;
        
    return out;
}

template<typename T> bool
SparseCOO<T>::judgeSymmetric() const
{
    for (uint k = 0; k < mNNZ; ++k) {
        if (mRowInd[k] > mColInd[k]) continue;
        bool symElemFound = false;
        for (uint l = 0; l < mNNZ; ++l) {
            if (mRowInd[l] != mColInd[k] || mColInd[l] != mRowInd[k]) continue;
            if (abs(mVal[l] - mVal[k]) < 1e-7) symElemFound = true;
            break;
        }
        if (!symElemFound) return false;
   }
   return true;
}

template<typename T> void
SparseCOO<T>::symmetrize()
{
    for (uint k = 0; k < mNNZ; ++k) {
        if (mRowInd[k] != mColInd[k]) {
            mVal.push_back(mVal[k]);
            mRowInd.push_back(mColInd[k]);
            mColInd.push_back(mRowInd[k]);
        }
    }
    mNNZ = mVal.size();
}

template<typename T>
class SparseCSR
{
public:
    int constructCSRFromCOO(int nRow, int nCol, const std::vector<uint>& vII, const std::vector<uint>& vJJ, const std::vector<T>& vVV);
    int constructCSRFromDiag(const std::vector<T>& vDiag);
    
    template<typename U>
    int constructCSRFrom(uint nRow, uint nCol, T a[], U ja[], U ia[]) 
    {
        mRowNum = nRow;
        mColNum = nCol;
        mNNZ = ia[nRow] - 1;
    
        mVal.resize(mNNZ);
        std::copy(a, a+mNNZ, mVal.begin());

        mColInd.resize(mNNZ);
        for (uint k = 0; k < mNNZ; ++k) mColInd[k] = ja[k];
    
        mRowPtr.resize(nRow+1);
        for (uint k = 0; k <= nRow; ++k) mRowPtr[k] = ia[k];

        return 0;
    }
    int toCOO(SparseCOO<T> &coo) const;
    
    uint mRowNum, mColNum;
    uint mNNZ;

    std::vector<T> mVal;
    std::vector<uint> mColInd;
    std::vector<uint> mRowPtr;

    SparseCSR();

    template<typename U> friend std::ostream& operator<< (std::ostream &out, const SparseCSR<U> &csr);
};

template<typename U>
std::ostream& operator<< (std::ostream &out, const SparseCSR<U> &csr) {
    out << "NNZ = " << csr.mNNZ << std::endl;
    out << "Values=(";
    for (uint k = 0; k < csr.mNNZ; ++k) 
        out << csr.mVal[k] << (k < csr.mNNZ-1 ? " " : "");
    out << ")\nColIdx=(";
    for (uint k = 0; k < csr.mNNZ; ++k)
        out << csr.mColInd[k] << (k < csr.mNNZ-1 ? " " : "");
    out << ")\nRowPtr=(";
    for (uint k = 0; k < csr.mRowNum; ++k)
        out << csr.mRowPtr[k] << (k < csr.mRowNum-1 ? " " : "");
    out << ")" << std::endl; 
        
    return out;
}

template<typename T> int 
SparseCSR<T>::toCOO(SparseCOO<T> &coo) const
{
    coo.mRowNum = this->mRowNum;
    coo.mColNum = this->mColNum;
    coo.mNNZ = this->mNNZ;
    coo.mVal.resize(mNNZ);
    std::copy(this->mVal.begin(), this->mVal.end(), coo.mVal.begin());
    coo.mColInd.resize(mNNZ);
    coo.mRowInd.resize(mNNZ);
    for (uint k = 0; k < mNNZ; ++k) coo.mColInd[k] = this->mColInd[k];
    
    for (uint r = 1; r <= mRowNum; ++r) {
        for (uint k = mRowPtr[r-1]; k < mRowPtr[r]; ++k) {
            coo.mRowInd[k-1] = r;
        }         
    }       

    return 0;
}

/*
template<typename T> 
template<typename U> int 
SparseCSR<T>::constructCSRFrom(uint nRow, uint nCol, T a[], U ja[], U ia[])
{
    mRowNum = nRow;
    mColNum = nCol;
    mNNZ = ia[nRow] - 1;
    
    mVal.resize(mNNZ);
    std::copy(a, a+mNNZ, mVal.begin());

    mColInd.resize(mNNZ);
    for (uint k = 0; k < mNNZ; ++k) mColInd[k] = ja[k];
    
    mRowPtr.resize(nRow+1);
    for (uint k = 0; k <= nRow; ++k) mRowPtr[k] = ia[k];

    return 0;
}
*/

template<typename T> int
SparseCOO<T>::toPackedSymmetricU(T ap[]) const
{
    // ap should has the size of N(N+1)/2
    uint N = mRowNum;
    for (uint k = 0; k < N*(N+1)/2; ++k) ap[k] = 0.0;
    
    for (uint k = 0; k < mNNZ; ++k) {
        uint i = mRowInd[k];
        uint j = mColInd[k];
        if (i > j) continue;
        T val = mVal[k];
        
        ap[i + j*(j-1)/2 - 1] = val;
    }

    return 0;
}

template<typename T> int
SparseCOO<T>::toCSR(SparseCSR<T>& csr) const
{
    csr.constructCSRFromCOO(mRowNum, mColNum, mRowInd, mColInd, mVal);

    return 0;
}

template<typename T> int
SparseCOO<T>::constructFromDiag(const std::vector<T>& vWW) 
{
    mNNZ = vWW.size();
    mRowNum = mColNum = vWW.size();
    mRowInd.resize(mNNZ);
    mColInd.resize(mNNZ);
    mVal.resize(mNNZ);    

    for (unsigned k = 0; k < mNNZ; ++k) {
        mRowInd[k] = k + 1;
        mColInd[k] = k + 1;
        mVal[k] = vWW[k];
    }
    
    return 0;
}

template<typename T>
SparseCSR<T>::SparseCSR()
{
    mRowNum = mColNum = mNNZ = 0;
}


template<typename T> int 
SparseCSR<T>::constructCSRFromDiag(const std::vector<T>& vDiag)
{
    mRowNum = mColNum = mNNZ = vDiag.size();
    mColInd.resize(mRowNum);
    mRowPtr.resize(mRowNum);
    for (uint i = 0; i < mRowNum; ++i) {
        mVal[i] = vDiag[i];
        mColInd[i] = i + 1;
        mRowPtr[i] = i + 1;    
    }
    mRowPtr[mRowNum] = mNNZ + 1;
    return 0;
}

template<typename T> int 
SparseCSR<T>::constructCSRFromCOO(int nRow, int nCol, const std::vector<uint>& vII,
                         const std::vector<uint>& vJJ, const std::vector<T>& vVV) 
{
    mRowNum = nRow; mColNum = nCol;
    std::vector<TupleCOO<T>> vCOO;
    for (uint k = 0; k < vVV.size(); ++k) {
        vCOO.push_back(TupleCOO<T>(vII[k],vJJ[k],vVV[k]));
    }
    
    std::sort(vCOO.begin(), vCOO.end(), std::less<TupleCOO<T>>()); 

    std::vector<uint> vColInd, vRowPtr;
    std::vector<T> vVal;
    int prevRow = 0;
    for (typename std::vector<TupleCOO<T>>::iterator iter = vCOO.begin(); iter != vCOO.end(); ++iter) {
        if (iter->i > iter->j) continue;    //only upper triangle
        int curRow = iter->i;
        
        if (prevRow < curRow) {
            for (int i = prevRow + 1; i < curRow; ++i) {
                vVal.push_back(0.);
                vColInd.push_back(i);
                vRowPtr.push_back(vVal.size());
            }
            if (iter->i < iter->j) {
                vVal.push_back(0.);
                vColInd.push_back(curRow);
                vRowPtr.push_back(vVal.size());
            }
        }

        vVal.push_back(iter->v);
        vColInd.push_back(iter->j);
        if (iter->i == iter->j) vRowPtr.push_back(vVal.size());
        
        prevRow = curRow;
    }       
    vRowPtr.push_back(vVal.size()+1);
    

    mNNZ = vVal.size();
    mVal = vVal;
    mColInd = vColInd;
    mRowPtr = vRowPtr;

    return 0;
}

template<typename T> void
printSparseMatrix(const SparseCOO<T>& coo, std::ostream& out)
{
    for (uint i = 0; i < coo.mNNZ; ++i) {
        out << coo.mRowInd[i] << ' ' << coo.mColInd[i] << ' ' << coo.mVal[i] << std::endl;
    }
}

} //end of namespace ZMath

#endif
