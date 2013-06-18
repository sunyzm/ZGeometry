#ifndef ZNUMERICS_SPARSE_MATRIX_H
#define ZNUMERICS_SPARSE_MATRIX_H
#include <vector>

namespace ZNumerics {

template<typename T>
class SparseCOO
{
public:
	std::vector<T> mVal;
	std::vector<unsigned int> mRowCoo, mColCoo;

	unsigned int mRowNum, mColNum, mNNZ;
};

}

#endif