#ifndef ZGEOM_DENSE_MATRIX_H
#define ZGEOM_DENSE_MATRIX_H
#include <fstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "common.h"
#include "VecN.h"

namespace ZGeom {

template<typename T>
class DenseMatrix
{
public:
	DenseMatrix() : mRow(0), mCol(0), mData(NULL) {}
	DenseMatrix(uint row, uint col) : mRow(row), mCol(col) { 
		mData = new T[row*col]; 
		for (uint i = 0; i < row*col; ++i) mData[i] = 0;
	}
    DenseMatrix(const DenseMatrix<T>& m2);
    DenseMatrix(DenseMatrix<T>&& m2);	
	DenseMatrix<T>& operator = (const DenseMatrix<T>& m2);
    DenseMatrix<T>& operator = (DenseMatrix<T>&& m2);
    DenseMatrix(const std::vector<VecN<T> >& vRowVecs);
    ~DenseMatrix() { delete []mData; }

	T& operator () (uint i, uint j) { return mData[i * mCol + j]; }	     // 0-based subscript access to matrix elements
	T operator () (uint i, uint j) const { return mData[i * mCol + j]; } // 0-based subscript access to matrix elements
		
	bool empty() const { return mData == NULL; }
	uint rowCount() const { return mRow; }
	uint colCount() const { return mCol; }
	bool symmetric() const { return mRow == mCol; }		
	T* raw_ptr() const { return mData; }
	T* raw_ptr_end() const { return mData + mRow*mCol; }
	void resize(uint row, uint col);
	void expand(uint newRow, uint newCol);
	void copyRows(const DenseMatrix<T>& m2, int startingRow = 0);
	void setRow(int rowIdx, const VecN<T>& vi);
	void print(const std::string& filepath) const;
	VecN<T> getRowVec(uint row) const;
	VecN<T> getColVec(uint col) const;
    std::vector<VecN<T> > toRowVecs() const;
    std::vector<VecN<T> > toColVecs() const;
	T rowSumNorm() const;
	T frobeniusNorm() const;
    DenseMatrix<T>& transpose();

	void write(const std::string& filename);
	void read(const std::string& filename);

private:
	T* mData;
	uint mRow, mCol;
};

typedef DenseMatrix<double> DenseMatrixd;
typedef DenseMatrix<float>  DenseMatrixs;

template<typename T>
inline void DenseMatrix<T>::print( const std::string& filepath ) const
{
	std::ofstream ofs(filepath.c_str());
	for (uint i = 0; i < mRow; ++i) {
		for (uint j = 0; j < mCol; ++j) {
			ofs << mData[i*mCol + j] << ' ';
		}
		ofs << '\n';
	}
	ofs.close();
}

template<typename T>
inline DenseMatrix<T>& DenseMatrix<T>::operator=( const DenseMatrix<T>& m2 )
{
	delete []mData;
	mRow = m2.mRow;
	mCol = m2.mCol;
	mData = new T[mRow*mCol];
	std::copy_n(m2.mData, mRow*mCol, mData);
	return *this;
}

template<typename T>
inline DenseMatrix<T>& DenseMatrix<T>::operator=(DenseMatrix<T>&& m2)
{
    delete []mData;
    mRow = m2.mRow;
    mCol = m2.mCol;
    mData = m2.mData;
    m2.mRow = m2.mCol = 0;
    m2.mData = nullptr;
    return *this;
}

template<typename T>
inline DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& m2) 
{
    *this = m2;
}

template<typename T>
inline DenseMatrix<T>::DenseMatrix(DenseMatrix<T>&& m2)
{
    *this = std::move(m2);
}

template<typename T>
inline DenseMatrix<T>::DenseMatrix(const std::vector<VecN<T> >& vRowVecs)
{
    mRow = vRowVecs.size();
    mCol = vRowVecs[0].size();
    mData = new T[mRow*mCol];
    for (int i = 0; i < mRow; ++i) {
        std::copy_n(vRowVecs[i].c_ptr(), mCol, mData + mCol*i);
    }
}

template<typename T>
inline DenseMatrix<T>& DenseMatrix<T>::transpose()
{
    T newData = new T[mRow*mCol];
    for (int i = 0; i < mRow; ++i)
        for (int j = 0; j < mCol; ++j)
            newData[i*mCol + j] = mData[j*mRow + i];
    delete[]mData;
    mData = newData;
    std::swap(mRow, mCol);
    return *this;
}

template<typename T>
inline void DenseMatrix<T>::write(const std::string& filename)
{
	std::ofstream ofs(filename.c_str(), std::ios::binary);
	int st = sizeof(T);
	ofs.write((const char*)&st, sizeof(st));
	ofs.write((const char*)&mRow, sizeof(mRow));
	ofs.write((const char*)&mCol, sizeof(mCol));		
	ofs.write((const char*)mData, sizeof(T)*mRow*mCol);
}

template<typename T>
inline void DenseMatrix<T>::read(const std::string& filename)
{
	std::ifstream ifs(filename.c_str(), std::ios::binary);
	int st;
	ifs.read((char*)&st, sizeof(st));
	if (sizeof(T) != st) throw std::runtime_error("DenseMatrix element type not compatible!");

	ifs.read((char*)&mRow, sizeof(mRow));
	ifs.read((char*)&mCol, sizeof(mCol));		
	mData = new T[mRow*mCol];
	ifs.read((char*)mData, sizeof(T)*mRow*mCol);
}

template<typename T>
inline T DenseMatrix<T>::rowSumNorm() const
{
	T rowSumMax(0);
	for (int i = 0; i < mRow; ++i) {
		T rowSum(0);
		for (int j = 0; j < mCol; ++j) rowSum += std::fabs(mData[i*mCol + j]);
		if (rowSum > rowSumMax) rowSumMax = rowSum;
	}
	return rowSumMax;
}

template<typename T>
VecN<T> DenseMatrix<T>::getRowVec(uint row) const
{
	assert(row < mRow);
	VecN<T> vec(mCol);
	for (int j = 0; j < mCol; ++j) vec[j] = mData[row*mCol + j];
	return vec;
}

template<typename T>
VecN<T> DenseMatrix<T>::getColVec(uint col) const
{
    assert(col < mCol);
    VecN<T> vec(mRow);
    for (int i = 0; i < mRow; ++i) vec[i] = mData[i*mCol + col];
    return vec;
}

template<typename T>
std::vector<VecN<T> > DenseMatrix<T>::toRowVecs() const
{
    std::vector<VecN<T> > vRows(mRow);
	for (int i = 0; i < mRow; ++i) vRows[i] = getRowVec(i);
    return vRows;
}

template<typename T>
std::vector<VecN<T> > DenseMatrix<T>::toColVecs() const
{
    std::vector<VecN<T> > vCols(mCol);
    for (int j = 0; j < mCol; ++j) vCols[j] = getColVec(j);
    return vCols;
}

template<typename T>
inline T DenseMatrix<T>::frobeniusNorm() const
{
	T sum(0);
	for (int i = 0; i < mRow*mCol; ++i)
		sum += mData[i] * mData[i];
	return std::sqrt(sum);
}

template<typename T>
inline void DenseMatrix<T>::copyRows(const DenseMatrix<T>& m2, int startingRow)
{
	if (this->mRow - startingRow < m2.mRow || this->mCol != m2.mCol) throw std::logic_error("DenseMatrix not compatible!");
	std::copy_n(m2.mData, m2.mRow * m2.mCol, mData + startingRow * this->mCol);
}

template<typename T>
void DenseMatrix<T>::setRow(int rowIdx, const VecN<T>& vi)
{
    assert(mData != nullptr);
	if (rowIdx < 0 || rowIdx >= (int)mRow || vi.size() != mCol)
		throw std::logic_error("Vector or index not compatible with matrix dimension!");
	std::copy_n(vi.c_ptr(), mCol, mData + rowIdx * mCol);
}

template<typename T>
inline void DenseMatrix<T>::resize(uint row, uint col) 
{ 
	delete []mData;
	mRow = row;
	mCol = col; 
	mData = new T[row*col]; 
}

template<typename T>
inline void DenseMatrix<T>::expand(uint newRow, uint newCol)
{
	if(newRow < mRow || newCol < mCol) throw std::runtime_error("DenseMatrix::expand error");
	if (newRow == mRow && newCol == mCol) return;
	double *newData = new double[newRow*newCol];
	for (int i = 0; i < newRow * newCol; ++i) newData[i] = 0;

	for (int i = 0; i < mRow; ++i) {
		for (int j = 0; j < mCol; ++j) {
			newData[i*newCol + j] = mData[i*mCol + j];
		}
	}

	mRow = newRow;
	mCol = newCol;
	delete []mData;
	mData = newData;
}

}	// end of namespace

#endif