#ifndef ZGEOM_DENSE_MATRIX_H
#define ZGEOM_DENSE_MATRIX_H

#include "common.h"
#include <fstream>
#include <string>
#include <algorithm>

namespace ZGeom
{
    template<typename T>
    class DenseMatrix
    {
    public:
        DenseMatrix() : mRow(0), mCol(0), mData(NULL) {}
        DenseMatrix(uint row, uint col) : mRow(row), mCol(col) { mData = new T[row*col]; }
        DenseMatrix(const DenseMatrix<T>& m2);
        const DenseMatrix<T>& operator = (const DenseMatrix<T>& m2);
        ~DenseMatrix() { delete []mData; }

        T& operator () (uint i, uint j) { return mData[i*mCol + j]; }
        T operator () (uint i, uint j) const { return mData[i*mCol + j]; }

        bool empty() const { return mData == NULL; }
        void resize(uint row, uint col) { delete []mData; mRow = row; mCol = col; mData = new T[row*col]; }
        void print(const std::string& filepath) const;
        uint rowCount() const { return mRow; }
        uint colCount() const { return mCol; }
        bool symmetric() const { return mRow == mCol; }
        T* raw_ptr() const { return mData; }
        T* raw_ptr_end() const { return mData + mRow*mCol; }

    private:
        T* mData;
        uint mRow, mCol;
    };

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
    inline const DenseMatrix<T>& ZGeom::DenseMatrix<T>::operator=( const DenseMatrix<T>& m2 )
    {
        delete []mData;
        mRow = m2.mRow;
        mCol = m2.mCol;
        mData = new T[mRow*mCol];
        std::copy_n(m2.mData, mRow*mCol, mData);
    }
    
    typedef DenseMatrix<double> DenseMatrixd;
    typedef DenseMatrix<float>  DenseMatrixs;
} //end of namespace



#endif