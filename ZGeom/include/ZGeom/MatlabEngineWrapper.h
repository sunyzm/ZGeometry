#ifndef ZGEOM_MATLAB_EINGINE_WRAPPER_H
#define ZGEOM_MATLAB_EINGINE_WRAPPER_H

#include <stdexcept>
#include <iostream>
#include <engine.h>

namespace ZGeom
{
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

    inline MatlabEngineWrapper::~MatlabEngineWrapper()
    {        
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