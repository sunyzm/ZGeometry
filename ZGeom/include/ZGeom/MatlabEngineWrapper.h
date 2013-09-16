#ifndef ZGEOM_MATLAB_EINGINE_WRAPPER_H
#define ZGEOM_MATLAB_EINGINE_WRAPPER_H

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <engine.h>

namespace ZGeom
{

class MatlabArrayWrapper
{
public:
	MatlabArrayWrapper(Engine* ep, double* data, int row, int col, const std::string& name) 
		: mName(name), mRow(row), mCol(col) 
	{
		mData = engGetVariable(ep, name.c_str());
		if (mData) mxDestroyArray(mData);

		mData = mxCreateDoubleMatrix(mRow, mCol, mxREAL);
		double *aa = mxGetPr(mData);
		std::copy_n(data, row*col, aa);
		engPutVariable(ep, mName.c_str(), mData);
	}

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

	void open(const char* startcmd = "\0");
	bool isOpened() const { return (m_ep != NULL); }
	Engine* getEngine() const { return m_ep; }
	void eval(const std::string& str) const { engEvalString(m_ep, str.c_str()); }
	const char* getOutput() const { return mBuffer; }
	void resizeBuffer(int bufSize);
	/* Get a variable with the specified name from MATLAB's workspace  */
	mxArray* getVariable(const char *name) const;
	/* Put a variable into MATLAB's workspace with the specified name */
	int putVariable(const char *var_name, const mxArray *ap) const;    
	void addVariable(double* data, int row, int col, const std::string& name) {
		mVariables.push_back(new MatlabArrayWrapper(m_ep, data, row, col, name));
	}

	void addVariable(int *data, int row, int col, const std::string& name) {
		double *buf = new double[row*col];
		for (int i = 0; i < row*col; ++i) buf[i] = data[i];
		addVariable(buf, row, col, name);
		delete []buf;
	}

	double* retrieveVariable(const std::string& name) {
		mxArray *arr = getVariable(name.c_str());
		return mxGetPr(arr);
	}

	void removeVariable(const std::string& varName) {
		for (auto iter = mVariables.begin(); iter != mVariables.end(); ++iter) {
			if ((*iter)->name() == varName) {
				delete *iter;
				mVariables.erase(iter);
				return;
			}
		}
	}

private:
	char *mBuffer;
	int  mBufSize;
	Engine* m_ep;
	std::vector<MatlabArrayWrapper*> mVariables;

};

inline MatlabEngineWrapper::MatlabEngineWrapper(int bufSize) 
{
	m_ep = NULL;
	mBuffer = NULL;
	resizeBuffer(bufSize);
}

inline MatlabEngineWrapper::~MatlabEngineWrapper()
{   
	for (MatlabArrayWrapper* arr : mVariables) delete arr;
	engClose(m_ep);
	delete []mBuffer;
	std::cout << "MatlabEngineWrapper destroyed!" << std::endl;
}

inline void MatlabEngineWrapper::open(const char* startcmd/*="\0"*/)
{
	m_ep = engOpen(startcmd);
	if (!m_ep) throw std::runtime_error("Fail to open Matlab engine");
	engOutputBuffer(m_ep, mBuffer, mBufSize);
}

inline void MatlabEngineWrapper::resizeBuffer(int bufSize) 
{
	mBufSize = bufSize;
	delete []mBuffer;
	mBuffer = new char[bufSize + 1]; 
	mBuffer[bufSize] = '\0';
}

inline mxArray* MatlabEngineWrapper::getVariable(const char *name) const
{
	assert(m_ep);
	return engGetVariable(m_ep, name);
}

inline int MatlabEngineWrapper::putVariable( const char *var_name, const mxArray *ap ) const
{
	assert(m_ep);
	return engPutVariable(m_ep, var_name, ap);
}

}



#endif