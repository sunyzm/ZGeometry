#include "SparseSolver.h"
#include "DenseMatrix.h"
#include "zassert.h"

namespace ZGeom{

void solveSparse(MatlabEngineWrapper& matlabEngine, const SparseMatrix<double>& matA, const VecNd& vec_b, VecNd& vec_x)
{
	/// solve sparse system Ax = b
	runtime_assert(matA.rowCount() == vec_b.size(), "Solve Ax=b - dimension of A and b not compatible");
	int dim1 = matA.rowCount(), dim2 = matA.colCount();	//A is a dim1*dim2 matrix

	matlabEngine.addColVec(vec_b.c_ptr(), dim1, "vec_b");
	matlabEngine.addSparseMat(matA, "matA");
	matlabEngine.eval("vec_x=matA\\vec_b;");
	double *px = matlabEngine.getDblVariablePtr("vec_x");
	vec_x = VecNd(px, dim2);
}

void solveSparseMultiColumn(MatlabEngineWrapper& matlabEngine, const SparseMatrix<double>& matA, const std::vector<VecNd>& vecB, std::vector<VecNd>& vecX)
{
	/// solve sparse system AX = B
	runtime_assert(matA.rowCount() == vecB[0].size(), "Solve AX=B : dimension of A and b not compatible");
	int dim1 = matA.rowCount(), dim2 = matA.colCount();	//A is a dim1*dim2 matrix, B is dim1*dim3, X is dim2*dim3
	int dim3 = (int)vecB.size();
	vecX.resize(dim3);

	matlabEngine.addSparseMat(matA, "matA");
	for (int i = 0; i < dim3; ++i) {
		matlabEngine.addColVec(vecB[i], "vec_b");
		matlabEngine.eval("vec_x=matA\\vec_b;");
		double *px = matlabEngine.getDblVariablePtr("vec_x");
		vecX[i] = VecNd(px, dim2);
	}
}

}