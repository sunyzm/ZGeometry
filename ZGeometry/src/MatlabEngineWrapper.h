#pragma once

#include <stdexcept>
#include <engine.h>

class MatlabEngineWrapper
{
public:
	MatlabEngineWrapper(int bufSize = 256);
	~MatlabEngineWrapper();
	
	void open(const char* startcmd = "\0");
	bool isOpened() const { return (m_ep != NULL); }
	Engine* getEngine() const { return m_ep; }
	void eval(const std::string& str) { engEvalString(m_ep, str.c_str()); }
	char* getOutput() const { return mBuffer; }
	void resizeBuffer(int bufSize);

	/* Get a variable with the specified name from MATLAB's workspace  */
	mxArray* getVariable(const char *name);
	/* Put a variable into MATLAB's workspace with the specified name */
	int putVariable(const char *var_name, const mxArray *ap);

private:
	char *mBuffer;
	int  mBufSize;
	Engine* m_ep;
};

inline MatlabEngineWrapper::MatlabEngineWrapper(int bufSize) 
{
	m_ep = NULL;
	mBuffer = NULL;
	resizeBuffer(bufSize);
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

inline MatlabEngineWrapper::~MatlabEngineWrapper()
{
    delete []mBuffer;
	engClose(m_ep);
    std::cout << "MatlabEngineWrapper destroyed!" << std::endl;
}

inline mxArray* MatlabEngineWrapper::getVariable(const char *name)
{
	assert(m_ep);
	return engGetVariable(m_ep, name);
}

inline int MatlabEngineWrapper::putVariable(const char *var_name, const mxArray *ap)
{
	assert(m_ep);
	return engPutVariable(m_ep, var_name, ap);
}