#ifndef ZGEOM_MATLAB_EINGINE_WRAPPER_H
#define ZGEOM_MATLAB_EINGINE_WRAPPER_H
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <engine.h>
#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "VecN.h"

namespace ZGeom {

class MatlabArrayWrapper
{
public:
	MatlabArrayWrapper(Engine* ep, double* data, int row, int col, bool rowMajor, const std::string& name) 
		: mName(name), mRow(row), mCol(col) 
	{
		mData = mxCreateDoubleMatrix(mRow, mCol, mxREAL);
		double *aa = mxGetPr(mData);
		if (rowMajor) {	
			// rowMajor needs transpose
			for (int j = 0; j < col; ++j)
				for (int i = 0; i < row; ++i)
					aa[j*row + i] = data[i*col + j];
		} else {	
			// colMajor just copy
			std::copy_n(data, row*col, aa);
		}
		
		engPutVariable(ep, mName.c_str(), mData);
	}

#if 0
	MatlabArrayWrapper(Engine *ep, int *data, int row, int col, bool rowMajor, const std::string& name) 
		: mName(name), mRow(row), mCol(col)
	{
		mData = engGetVariable(ep, name.c_str());
		if (mData) mxDestroyArray(mData);
	
		mwSize dims[] = {row, col};
		mData = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
		int *aa = (int*)mxGetData(mData);
		if (rowMajor) {
			for (int j = 0; j < col; ++j)
				for (int i = 0; i < row; ++i)
					aa[j*row + i] = data[i*col + j];
		} else {
			std::copy_n(data, row*col, aa);
		}
		engPutVariable(ep, mName.c_str(), mData);
	}
#endif

	~MatlabArrayWrapper() { mxDestroyArray(mData); }
	const std::string& name() const { return mName; }

private:
	mxArray* mData;
	std::string mName;
	int mRow, mCol;
};

class MatlabEngineWrapper
{
public:
	MatlabEngineWrapper(int bufSize = 256);
	~MatlabEngineWrapper();

	void open(const char* startcmd = NULL);
	bool isOpened() const { return (m_ep != NULL); }
	Engine* getEngine() const { return m_ep; }
	void eval(const std::string& str) const { engEvalString(m_ep, str.c_str()); }
	const char* getOutput() const { return mBuffer; }
	void resizeBuffer(int bufSize);
	mxArray* getVariable(const char *name) const; /* Get a variable with the specified name from MATLAB's workspace  */
	int putVariable(const char *var_name, const mxArray *ap) const; /* Put a variable into MATLAB's workspace with the specified name */
	
	void addArray(double* data, int row, int col, bool rowMajor, const std::string& name);
	void addArray(int *data, int row, int col, bool rowMajor, const std::string& name);
	void addDoubleScalar(double value, const std::string& name);
	void addColVec(double *data, int row, const std::string& name);
	void addColVec(int *data, int row, const std::string& name);
	void addColVec(const VecNd& data, const std::string& name);
	void addSparseMat(int *ii, int *jj, double *ss, int m, int n, int nnz, const std::string& name);
	void addSparseMat(const SparseMatrix<double>& mat, const std::string& varName);
	void addDenseMat(const DenseMatrix<double>& mat, const std::string& varName, bool asTranspose = false);

	double* getDblVariablePtr(const std::string& name);
	int* getIntVariablePtr(const std::string& name);
	void getSparseMat(const std::string& name, SparseMatrix<double>& mat);
    ZGeom::DenseMatrixd getDenseMat(const std::string& name, bool getTranspose = false);
	void removeVariable(const std::string& varName);
	
private:
	char *mBuffer;
	int  mBufSize;
	Engine* m_ep;
	std::vector<MatlabArrayWrapper*> mVariables;
};

}	// end of namespace
#endif