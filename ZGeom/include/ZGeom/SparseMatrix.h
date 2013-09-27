#ifndef ZGEOM_SPARSE_MATRIX_H
#define ZGEOM_SPARSE_MATRIX_H

#include <iostream>
#include <vector>
#include "common.h"
#include "SparseMatrixCOO.h"
#include "SparseMatrixCSR.h"

namespace ZGeom
{
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
		bool operator< (const MatElem<T>& t2) const {
			return this->mRow < t2.mRow || (this->mRow == t2.mRow && this->mCol < t2.mCol);
		}

		MatElem() : mRow(0), mCol(0), mVal(0.) {}
		MatElem(uint ii, uint jj, T vv) : mRow(ii), mCol(jj), mVal(vv) {}
		MatElem(const MatElem<T>& elem) : mRow(elem.mRow), mCol(elem.mCol), mVal(elem.mVal) {}
		MatElem<T>& operator = (const MatElem<T>& elem) {
			mRow = elem.mRow;
			mCol = elem.mCol;
			mVal = elem.mVal;
			return *this;
		}
		uint row() const { return mRow; }
		uint col() const { return mCol; }
		T    val() const { return mVal; }

	private:   
		uint mRow, mCol;
		T mVal;
	};

	template<typename T>
	class SparseMatrix
	{
	public:
		SparseMatrix() : mRowCount(0), mColCount(0), mNonzeroCount(0) {}
		SparseMatrix(int row, int col) : mRowCount(row), mColCount(col), mNonzeroCount(0) {}

		uint rowCount() const { return mRowCount; }
		uint colCount() const { return mColCount; }
		uint nonzeroCount() const { return mElements.size(); }

		T getElemByIndex(uint index) const;
		T& getElemByIndex(uint index);
		T getElem(uint row, uint col) const;
		T& getElem(uint row, uint col);
		T operator() (uint row, uint col) const; 
		T& operator() (uint row, uint col);
		void insertElem(uint row, uint col, T val);
		void removeElem(uint row, uint col);   

		void copyElements(const SparseMatrix<T>& mat2);

		void computeSubMatrix(const std::vector<int>& vSelected, SparseMatrix<T>& subMat) const;

		template<typename F>
		void getDiagonal(std::vector<F>& diag) const;

		template<typename F>
		void convertFromDiagonal(const std::vector<F>& diag);    

		template<typename U, typename F>
		void convertFromCOO(uint rowCount, uint colCount, const std::vector<U>& rowInd, std::vector<U>& colInd, const std::vector<F>& val);

		template<typename U, typename F>
		void convertToCOO(std::vector<U>& rowInd, std::vector<U>& colInd, std::vector<F>& val, MatrixForm form) const;

		template<typename U, typename F>
		void convertToCOO(U* rowInd, U* colInd, F* val, int& nonzeroCount, MatrixForm form) const;    

		template<typename U, typename F> 
		void convertFromCSR(uint rowCount, uint colCount, F nzVal[], U colIdx[], U rowPtr[]);

		template<typename U, typename F> 
		void convertToCSR(std::vector<F>& nzVal, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form = MAT_UPPER) const;

		template<typename U, typename F> 
		void convertToCSR(SparseMatrixCSR<F,U>& MatCSR, MatrixForm form = MAT_UPPER) const;

		template<typename F> 
		void convertFromFull(F* fullMat, double sparse_eps = 1e-10);

		template<typename F> 
		void convertToFull(F* fullMat, MatrixForm form) const;

		void truncate(MatrixForm form = MAT_UPPER); // truncate matrix 
		void symmetrize();              // turn upper or lower matrix into symmetric one
		void fillEmptyDiagonal();       // fill empty diagonal elements with 0
		void setToIdentity(uint order);     // turn into an identity matrix
		void print(std::ostream& out) const;
		void print(const std::string& paht) const;
		void read(std::istream& in);

		const SparseMatrix<T>& operator *= (double coeff);
		friend VecN<T> mulMatVec(const SparseMatrix<T>& mat, const VecN<T>& vec, bool matIsSym);

	private:
		bool testNoEmptyRow() const;
		bool testSymmetric(double eps = 1e-7) const;

		std::vector< MatElem<T> > mElements;
		uint mRowCount;
		uint mColCount;
		uint mNonzeroCount;
	};
} // end of namespace ZGeom

#include "SparseMatrix.inl"

#endif
