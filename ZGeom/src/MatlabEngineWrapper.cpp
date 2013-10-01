#include "MatlabEngineWrapper.h"
#include <iostream>


ZGeom::MatlabEngineWrapper::MatlabEngineWrapper( int bufSize )
{
	m_ep = NULL;
	mBuffer = NULL;
	resizeBuffer(bufSize);
}

ZGeom::MatlabEngineWrapper::~MatlabEngineWrapper()
{
	for (MatlabArrayWrapper* arr : mVariables) delete arr;
	engClose(m_ep);
	delete []mBuffer;
	std::cout << "MatlabEngineWrapper destroyed!" << std::endl;
}

void ZGeom::MatlabEngineWrapper::open( const char* startcmd/*="\0"*/ )
{
	m_ep = engOpen(startcmd);
	if (!m_ep) throw std::runtime_error("Fail to open Matlab engine");
	engOutputBuffer(m_ep, mBuffer, mBufSize);
}

void ZGeom::MatlabEngineWrapper::resizeBuffer( int bufSize )
{
	mBufSize = bufSize;
	delete []mBuffer;
	mBuffer = new char[bufSize + 1]; 
	mBuffer[bufSize] = '\0';
}

void ZGeom::MatlabEngineWrapper::removeVariable( const std::string& varName )
{
	for (auto iter = mVariables.begin(); iter != mVariables.end(); ++iter) {
		if ((*iter)->name() == varName) {
			delete *iter;
			mVariables.erase(iter);
			return;
		}
	}
}

void ZGeom::MatlabEngineWrapper::addSparseMat( int *ii, int *jj, double *ss, int m, int n, int nnz, const std::string& name )
{
	addColVec(ii, nnz, "ii");
	addColVec(jj, nnz, "jj");
	addColVec(ss, nnz, "ss");
	addDoubleScalar(m, "m");
	addDoubleScalar(n, "n");
	eval((name + "=sparse(ii,jj,ss,m,n)").c_str());
	removeVariable("ii");
	removeVariable("ii");
	removeVariable("ss");
	removeVariable("m");
	removeVariable("n");
}

void ZGeom::MatlabEngineWrapper::addSparseMat( const SparseMatrix<double>& mat, const std::string& varName )
{
	std::vector<double> ss;
	std::vector<int> ii, jj;
	mat.convertToCOO(ii, jj, ss, MAT_FULL);
	int m = mat.rowCount(), n = mat.colCount();
	int nnz = (int)ss.size();
	addSparseMat(&ii[0], &jj[0], &ss[0], m, n, nnz, varName);
}

mxArray* ZGeom::MatlabEngineWrapper::getVariable( const char *name ) const
{
	assert(m_ep);
	return engGetVariable(m_ep, name);
}

int ZGeom::MatlabEngineWrapper::putVariable( const char *var_name, const mxArray *ap ) const
{
	assert(m_ep);
	return engPutVariable(m_ep, var_name, ap);
}

void ZGeom::MatlabEngineWrapper::addVariable( double* data, int row, int col, bool rowMajor, const std::string& name )
{
	mVariables.push_back(new MatlabArrayWrapper(m_ep, data, row, col, rowMajor, name));
}

void ZGeom::MatlabEngineWrapper::addVariable( int *data, int row, int col, bool rowMajor, const std::string& name )
{
	double *buf = new double[row*col];
	for (int i = 0; i < row*col; ++i) buf[i] = data[i];
	addVariable(buf, row, col, rowMajor, name);
	delete []buf;
}

void ZGeom::MatlabEngineWrapper::addColVec( double *data, int row, const std::string& name )
{
	addVariable(data, row, 1, false, name);
}

void ZGeom::MatlabEngineWrapper::addColVec( int *data, int row, const std::string& name )
{
	addVariable(data, row, 1, false, name);
}

double* ZGeom::MatlabEngineWrapper::getDblVariablePtr( const std::string& name )
{
	mxArray *arr = getVariable(name.c_str());
	return mxGetPr(arr);
}

int * ZGeom::MatlabEngineWrapper::getIntVariablePtr( const std::string& name )
{
	mxArray *arr = getVariable(name.c_str());
	return (int*)mxGetData(arr);
}

void ZGeom::MatlabEngineWrapper::addDoubleScalar( double value, const std::string& name )
{
	addVariable(&value, 1, 1, false, name);
}
