#ifndef ZGEOM_DENSE_MATRIX_H
#define ZGEOM_DENSE_MATRIX_H
#include <fstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "common.h"

namespace ZGeom
{
	template<typename T>
	class DenseMatrix
	{
	public:
		DenseMatrix() : mRow(0), mCol(0), mData(NULL) {}
		DenseMatrix(uint row, uint col) : mRow(row), mCol(col) { 
			mData = new T[row*col]; 
			for (int i = 0; i < row*col; ++i) mData[i] = 0;
		}
		DenseMatrix(const DenseMatrix<T>& m2);
		DenseMatrix(DenseMatrix<T>&& m2);
		const DenseMatrix<T>& operator = (const DenseMatrix<T>& m2);
		~DenseMatrix() { delete []mData; }

		T& operator () (uint i, uint j) { return mData[i*mCol + j]; }	// 0-based subscript access to matrix elements
		T operator () (uint i, uint j) const { return mData[i*mCol + j]; } // 0-based subscript access to matrix elements

		bool empty() const { return mData == NULL; }
		uint rowCount() const { return mRow; }
		uint colCount() const { return mCol; }
		bool symmetric() const { return mRow == mCol; }		
		T* raw_ptr() const { return mData; }
		T* raw_ptr_end() const { return mData + mRow*mCol; }
		void resize(uint row, uint col);
		void expand(uint newRow, uint newCol);
		void copyRows(const DenseMatrix<T>& m2, int startingRow = 0);
		void print(const std::string& filepath) const;

	private:
		T* mData;
		uint mRow, mCol;
	};

	typedef DenseMatrix<double> DenseMatrixd;
	typedef DenseMatrix<float>  DenseMatrixs;

	template<typename T>
	inline void ZGeom::DenseMatrix<T>::print( const std::string& filepath ) const
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
	inline DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& m2) : mRow(m2.mRow), mCol(m2.mCol)
	{
		mData = new T[mRow*mCol];
		std::copy_n(m2.mData, mRow*mCol, mData);
	}

	template<typename T>
	inline DenseMatrix<T>::DenseMatrix(DenseMatrix<T>&& m2) : mRow(m2.mRow), mCol(m2.mCol)
	{
		mData = m2.mData;
		m2.mData = NULL;
		m2.mRow = m2.mCol = 0;
	}

	template<typename T>
	inline const DenseMatrix<T>& ZGeom::DenseMatrix<T>::operator=( const DenseMatrix<T>& m2 )
	{
		delete []mData;
		mRow = m2.mRow;
		mCol = m2.mCol;
		mData = new T[mRow*mCol];
		std::copy_n(m2.mData, mRow*mCol, mData);
	}

	template<typename T>
	inline void DenseMatrix<T>::copyRows(const DenseMatrix<T>& m2, int startingRow)
	{
		if (this->mRow - startingRow < m2.mRow || this->mCol != m2.mCol) throw std::logic_error("DenseMatrix not compatible!");
		std::copy_n(m2.mData, m2.mRow * m2.mCol, mData + startingRow * this->mCol);
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

} //end of namespace

#endif