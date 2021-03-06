#include "MatlabEngineWrapper.h"
#include <iostream>
#include <regex>

ZGeom::MatlabEngineWrapper::MatlabEngineWrapper( int bufSize )
{
	m_ep = NULL;
	mBuffer = NULL;
	resizeBuffer(bufSize);
}

ZGeom::MatlabEngineWrapper::~MatlabEngineWrapper()
{
    if (isOpened()) close();	
}

void ZGeom::MatlabEngineWrapper::open( const char* startcmd/*=NULL*/ )
{
    if (isOpened()) return;

	m_ep = engOpen(NULL);
	if (!m_ep) {
		std::cerr << "FAIL to open Matlab engine!" << std::endl;
		throw std::runtime_error("Fail to open Matlab engine");
	}
	engOutputBuffer(m_ep, mBuffer, mBufSize);
}

void ZGeom::MatlabEngineWrapper::close()
{
    for (MatlabArrayWrapper* arr : mVariables) delete arr;
    engClose(m_ep);
    delete[]mBuffer;
    m_ep = NULL;
}

void ZGeom::MatlabEngineWrapper::eval(const std::string& str) const
{
    engEvalString(m_ep, str.c_str());
}

void ZGeom::MatlabEngineWrapper::resizeBuffer( int bufSize )
{
    delete[]mBuffer;
    mBufSize = bufSize;	
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

void ZGeom::MatlabEngineWrapper::addArray( double* data, int row, int col, bool rowMajor, const std::string& name )
{
	removeVariable(name);
	mVariables.push_back(new MatlabArrayWrapper(m_ep, data, row, col, rowMajor, name));
}

void ZGeom::MatlabEngineWrapper::addArray( int *data, int row, int col, bool rowMajor, const std::string& name )
{
	double *buf = new double[row*col];
	for (int i = 0; i < row*col; ++i) buf[i] = data[i];
	addArray(buf, row, col, rowMajor, name);
	delete []buf;
}

void ZGeom::MatlabEngineWrapper::addColVec( double *data, int row, const std::string& name )
{
	addArray(data, row, 1, false, name);
}

void ZGeom::MatlabEngineWrapper::addColVec( int *data, int row, const std::string& name )
{
	addArray(data, row, 1, false, name);
}

void ZGeom::MatlabEngineWrapper::addColVec(const VecNd& data, const std::string& name )
{
	addColVec(data.c_ptr(), data.size(), name);
}

double* ZGeom::MatlabEngineWrapper::getDblVariablePtr( const std::string& name )
{
	mxArray *arr = getVariable(name.c_str());
	return mxGetPr(arr);
}

void ZGeom::MatlabEngineWrapper::addDoubleScalar( double value, const std::string& name )
{
	addArray(&value, 1, 1, false, name);
}

void ZGeom::MatlabEngineWrapper::getSparseMat( const std::string& name, SparseMatrix<double>& mat )
{
	eval(("[m,n]=size(" + name + ")").c_str());
	eval(("nonz=nnz("   + name + ")").c_str());
	eval(("[row,col,v] = find(" + name + ")").c_str());
	int m = (int)*getDblVariablePtr("m");
	int n = (int)*getDblVariablePtr("n");
	int nnz = (int)*getDblVariablePtr("nonz");
	double *row = getDblVariablePtr("row");
	double *col = getDblVariablePtr("col");
	double *v = getDblVariablePtr("v");

	std::vector<int> vRow(nnz), vCol(nnz);
	for (int i = 0; i < nnz; ++i) {
		vRow[i] = (int)row[i];
		vCol[i] = (int)col[i];
	}
	std::vector<double> vVal(nnz);
	std::copy_n(v, nnz, vVal.begin());

	mat.convertFromCOO(m, n, vRow, vCol, vVal);
}

void ZGeom::MatlabEngineWrapper::addDenseMat(const DenseMatrixd& mat, const std::string& varName, bool addAsTranspose)
{
    // if addAsTranspose, add as col-major; if not addAsTranspose, add as row-major
    if (addAsTranspose)
        addArray(mat.raw_ptr(), mat.colCount(), mat.rowCount(), false, varName);  
    else 
        addArray(mat.raw_ptr(), mat.rowCount(), mat.colCount(), true, varName);  
}

ZGeom::DenseMatrixd ZGeom::MatlabEngineWrapper::getDenseMat(const std::string& name, bool getTranspose /*= false*/)
{
    eval(("[m,n]=size(" + name + ")").c_str());
    int m = (int)*getDblVariablePtr("m");
    int n = (int)*getDblVariablePtr("n");
    double* pMat = getDblVariablePtr(name);
    
    if (getTranspose) 
    {
        ZGeom::DenseMatrixd mat(n, m);
        std::copy_n(pMat, m*n, mat.raw_ptr());
        return mat;
    }
    else 
    {
        ZGeom::DenseMatrixd mat(m, n);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                mat.raw_ptr()[i*n + j] = pMat[j*m + i];
        return mat;
    }    
}

std::string ZGeom::MatlabEngineWrapper::getMatlabVersion()
{
    eval("matlabver");
    std::string verstr = std::string(getOutput());
    std::regex ver_regex("R[0-9]{4}[ab]");
    std::smatch ver_match;
    std::regex_search(verstr, ver_match, ver_regex);
    if (ver_match.size() > 0) {
        std::ssub_match sub_match = ver_match[0];
        return sub_match.str();
    }
    else return "";
}
