#include "SparseMatrix.h"
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>

namespace ZGeom
{
	template<typename T>
	inline T SparseMatrix<T>::getElemByIndex(uint index) const
	{
		assert(index < mNonzeroCount);
		return mElements[index].mVal;
	}

	template<typename T>
	inline T& SparseMatrix<T>::getElemByIndex(uint index)
	{
		assert(index < mNonzeroCount);
		return mElements[index].mVal;
	}

	template<typename T>
	inline T SparseMatrix<T>::getElem(uint row, uint col) const
	{
		assert(row <= mRowCount && col <= mColCount);

		for (uint k = 0; k < mElements.size(); ++k) {
			if (mElements[k].mRow == row && mElements[k].mCol == col) 
				return mElements[k].mVal;
		}
		return 0;   
	}

	template<typename T>
	inline T& SparseMatrix<T>::getElem(uint row, uint col)
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
		return getElem(row, col);
	}

	template<typename T>
	inline T SparseMatrix<T>::operator() (uint row, uint col) const
	{
		return getElem(row, col);
	}

	template<typename T>
	inline void SparseMatrix<T>::insertElem( uint row, uint col, T val )
	{
		mElements.push_back(MatElem(row, col, val));
	}

	template<typename T>
	inline void SparseMatrix<T>::removeElem( uint row, uint col )
	{
		for (auto iter = mElements.begin(); iter != mElements.end(); ) {
			if (iter->row() == row && iter->col() == col) {
				iter = mElements.erase(iter);
				continue;
			}
			++iter;
		}
	}

	template<typename T>
	const SparseMatrix<T>& SparseMatrix<T>::operator *= (double coeff)
	{
		Concurrency::parallel_for_each(mElements.begin(), mElements.end(), [=](MatElem& e){
			e->mVal *= coeff;
		});
	}

	template<typename T>
	inline void SparseMatrix<T>::setToIdentity(uint order)
	{
		mRowCount = mColCount = mNonzeroCount = order;
		mElements.clear();
		for (uint k = 0; k < order; ++k) {
			mElements.push_back(MatElem<T>(k+1, k+1, 1.0));   
		}    
	}

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
		for (typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin(); iter != mElements.end(); ++iter) {
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
					if (fabs(mElements[l].val() - mElements[k].val()) < eps) symElemFound = true;
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
	template<typename F>
	void SparseMatrix<T>::getDiagonal(std::vector<F>& diag) const
	{
		assert(mRowCount == mColCount);

		diag.resize(mRowCount, 0.0);
		typename std::vector< MatElem<T> >::const_iterator iter = mElements.begin();
		for (; iter != mElements.end(); ++iter) {
			if (iter->row() == iter->col()) 
				diag[iter->row() - 1] = iter->val();
		}
	}

	template<typename T>
	template<typename F>
	void SparseMatrix<T>::convertFromDiagonal(const std::vector<F>& diag)
	{
		mRowCount = mColCount = mNonzeroCount = diag.size();
		mElements.resize(mRowCount);
		for (uint k = 0; k < mRowCount; ++k) {
			mElements[k].mRow = mElements[k].mCol = k + 1;
			mElements[k].mVal = diag[k];
		}
	}

	template<typename T>
	template<typename U, typename F>
	void SparseMatrix<T>::convertFromCOO( uint rowCount, uint colCount, const std::vector<U>& rowInd, std::vector<U>& colInd, const std::vector<F>& val )
	{
		assert(rowInd.size() == colInd.size() && rowInd.size() == val.size());

		mRowCount = rowCount; 
		mColCount = colCount;
		mNonzeroCount = val.size();
		
		mElements.clear();
		mElements.reserve(mNonzeroCount);
		for (uint k = 0; k < mNonzeroCount; ++k) {
			assert(rowInd[k] > 0 && colInd[k] > 0);
			assert(rowInd[k] <= rowCount && colInd[k] <= colCount);
			mElements.push_back(MatElem<T>(rowInd[k], colInd[k], val[k]));
		}
	}

	template<typename T>
	template<typename U, typename F>
	inline void SparseMatrix<T>::convertToCOO(std::vector<U>& rowInd, std::vector<U>& colInd,std::vector<F>& val, MatrixForm form) const
	{
		assert (mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR

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
		assert (mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR

		nonzeroCount = 0;
		typename std::vector< MatElem<T> >::const_iterator iter = sortedElements.begin();
		for (; iter != sortedElements.end(); ++iter) {
			if (form == MAT_UPPER && iter->row() > iter->col()) continue;   // only preserve upper triangle elements
			if (form == MAT_LOWER && iter->row() < iter->col()) continue;   // only preserve lower triangle elements

			rowInd[nonzeroCount] = static_cast<U>(iter->row());
			colInd[nonzeroCount] = static_cast<U>(iter->col());
			val[nonzeroCount]    = iter->val();

			nonzeroCount++;
		}
	}

	template<typename T>
	template<typename U, typename F>
	inline void SparseMatrix<T>::convertToCSR(std::vector<F>& nzVal, std::vector<U>& colIdx, std::vector<U>& rowPtr, MatrixForm form /*=MAT_UPPER*/) const
	{
		assert (mRowCount == mColCount || form == MAT_FULL); // only square matrix has upper or lower CSR        
		assert (testNoEmptyRow());          // CSR format requires at least one element in each row

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
	inline void SparseMatrix<T>::convertToCSR(SparseMatrixCSR<F,U>& matCSR, MatrixForm form) const
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

			fullMat[(iter->row()-1) * mColCount + (iter->col()-1)] = iter->val();    
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
			for (uint k = (uint)rowPtr[r-1]; k < (uint)rowPtr[r]; ++k) {
				mElements[k-1].mRow = r;
				mElements[k-1].mCol = colIdx[k-1];
				mElements[k-1].mVal = nzVal[k-1];
			}
		}
	}
}// end of namespace ZGeom