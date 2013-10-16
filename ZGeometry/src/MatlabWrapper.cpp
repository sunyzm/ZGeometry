#include "MatlabWrapper.h"
#include <stdexcept>
using namespace std;

void matlab_cgls( Engine* ep, const vector<vector<double> >& SGW, const vector<double>& b, vector<double>& sol )
{
	//solve A*sol = B

	assert(SGW.size() == b.size());

	int sizeA = SGW.size(), sizeF = SGW[0].size();
	sol.resize(sizeF);

	mxArray *AA, *BB;
	AA = mxCreateDoubleMatrix(sizeA, sizeF, mxREAL);
	BB = mxCreateDoubleMatrix(sizeA, 1, mxREAL);
	double *pAA = mxGetPr(AA), *pBB = mxGetPr(BB);

	for (int i = 0; i < sizeA; ++i)
	{
		for (int j = 0; j < sizeF; ++j)
		{
			pAA[i + j*sizeA] = SGW[i][j];
		}
	}

	for (int j = 0; j < sizeA; ++j)
	{
		pBB[j] = b[j];
	}

	engPutVariable(ep, "A", AA);
	engPutVariable(ep, "b", BB);

	engEvalString(ep, "evals = cgls(A, b);");

	mxArray *SS = engGetVariable(ep, "evals");
	double *pSS = mxGetPr(SS);
	for (int j = 0; j < sizeF; ++j)
	{
		sol[j] = pSS[j];
	}

	mxDestroyArray(AA);
	mxDestroyArray(BB);
}

void matlab_scgls( Engine* ep, const vector<vector<double> >& SGW, const vector<double>& b, vector<double>& sol )
{
	//solve A*sol = B

	assert(SGW.size() == b.size());

	int sizeA = SGW.size(), sizeF = SGW[0].size();
	sol.resize(sizeF);

	const double sEpsilon = 1e-5;
	vector<double> vII, vJJ, vSS;
	for (int i = 0; i < sizeA; ++i)
	{
		for (int j = 0; j < sizeF; ++j)
		{
			if ( fabs(SGW[i][j]) > sEpsilon )
			{
				vII.push_back(i+1);
				vJJ.push_back(j+1);
				vSS.push_back(SGW[i][j]);
			}
		}
	}

	int sizeS = vSS.size();
	mxArray *II, *JJ, *SS, *BB, *DimM, *DimN;
	II = mxCreateDoubleMatrix(sizeS, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(sizeS, 1, mxREAL);
	SS = mxCreateDoubleMatrix(sizeS, 1, mxREAL);
	BB = mxCreateDoubleMatrix(sizeA, 1, mxREAL);
	DimM = mxCreateDoubleMatrix(1, 1, mxREAL);
	DimN = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *pII = mxGetPr(II), *pJJ = mxGetPr(JJ), 
		*pSS = mxGetPr(SS), *pBB = mxGetPr(BB),
		*pDimM = mxGetPr(DimM), *pDimN = mxGetPr(DimN);

	pDimM[0] = sizeA; pDimN[0] = sizeF;

	for (int n = 0; n < sizeS; ++n)
	{
		pII[n] = vII.at(n);
		pJJ[n] = vJJ.at(n);
		pSS[n] = vSS.at(n);
	}

	for (int j = 0; j < sizeA; ++j)
	{
		pBB[j] = b[j];
	}

	engPutVariable(ep, "II", II);
	engPutVariable(ep, "JJ", JJ);
	engPutVariable(ep, "SS", SS);
	engPutVariable(ep, "DimM", DimM);
	engPutVariable(ep, "DimN", DimN);
	engPutVariable(ep, "BB", BB);

	engEvalString(ep, "evals = scgls(II, JJ, SS, DimM, DimN, BB);");

	mxArray *evals = engGetVariable(ep, "evals");
	double *pEvals = mxGetPr(evals);
	for (int j = 0; j < sizeF; ++j)
	{
		sol[j] = pEvals[j];
	}

	mxDestroyArray(II);
	mxDestroyArray(JJ);
	mxDestroyArray(SS);
	mxDestroyArray(DimM);
	mxDestroyArray(DimN);
	mxDestroyArray(BB);
}

void MatlabWrapper::SparseBiConjugateGradient( int m_size, int n_size, const std::vector<int>& vI, const std::vector<int>& vJ, const std::vector<double>& vS, const std::vector<double>& vY, std::vector<double>& vX )
{
	if (vI.size() != vJ.size() || vI.size() != vS.size() || vY.size() != m_size)
		throw logic_error("Error: MatlabWrapper::SparseBiConjugateGradient; incompatible parameter");
	
	int sizeNZ = vS.size();
	mxArray *II, *JJ, *SS, *BB, *DimM, *DimN;
	II = mxCreateDoubleMatrix(sizeNZ, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(sizeNZ, 1, mxREAL);
	SS = mxCreateDoubleMatrix(sizeNZ, 1, mxREAL);
	BB = mxCreateDoubleMatrix(m_size, 1, mxREAL);
	DimM = mxCreateDoubleMatrix(1, 1, mxREAL);
	DimN = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *pII = mxGetPr(II), *pJJ = mxGetPr(JJ), 
		   *pSS = mxGetPr(SS), *pBB = mxGetPr(BB),
		   *pDimM = mxGetPr(DimM), *pDimN = mxGetPr(DimN);

	pDimM[0] = m_size; pDimN[0] = n_size;

	for (int n = 0; n < sizeNZ; ++n)
	{
		pII[n] = vI[n];
		pJJ[n] = vJ[n];
		pSS[n] = vS[n];
	}

	for (int j = 0; j < m_size; ++j)
	{
		pBB[j] = vY[j];
	}

	engPutVariable(m_ep, "II", II);
	engPutVariable(m_ep, "JJ", JJ);
	engPutVariable(m_ep, "SS", SS);
	engPutVariable(m_ep, "DimM", DimM);
	engPutVariable(m_ep, "DimN", DimN);
	engPutVariable(m_ep, "BB", BB);

	engEvalString(m_ep, "evals = sparse_bi_conj_grad(II, JJ, SS, DimM, DimN, BB);");

	mxArray *evals = engGetVariable(m_ep, "evals");
	double *pEvals = mxGetPr(evals);
	vX.resize(n_size, 0);
	for (int j = 0; j < n_size; ++j)
	{
		vX[j] = pEvals[j];
	}

	mxDestroyArray(II);
	mxDestroyArray(JJ);
	mxDestroyArray(SS);
	mxDestroyArray(DimM);
	mxDestroyArray(DimN);
	mxDestroyArray(BB);
}

void MatlabWrapper::SparseMatVecMul( int m_size, int n_size, const std::vector<int>& vI, const std::vector<int>& vJ, const std::vector<double>& vS, const std::vector<double>& vX, std::vector<double>& vY )
{
	if (vI.size() != vJ.size() || vI.size() != vS.size() || vX.size() != n_size)
		throw logic_error("Error: MatlabWrapper::SparseMatVecMul; incompatible parameter");

	int sizeNZ = vS.size();
	mxArray *II, *JJ, *SS, *XX, *DimM, *DimN;
	II = mxCreateDoubleMatrix(sizeNZ, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(sizeNZ, 1, mxREAL);
	SS = mxCreateDoubleMatrix(sizeNZ, 1, mxREAL);
	XX = mxCreateDoubleMatrix(n_size, 1, mxREAL);
	DimM = mxCreateDoubleMatrix(1, 1, mxREAL);
	DimN = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *pII = mxGetPr(II), *pJJ = mxGetPr(JJ), 
		   *pSS = mxGetPr(SS), *pXX = mxGetPr(XX),
		   *pDimM = mxGetPr(DimM), *pDimN = mxGetPr(DimN);

	pDimM[0] = m_size; pDimN[0] = n_size;

	for (int n = 0; n < sizeNZ; ++n)
	{
		pII[n] = vI[n];
		pJJ[n] = vJ[n];
		pSS[n] = vS[n];
	}

	for (int j = 0; j < n_size; ++j)
	{
		pXX[j] = vX[j];
	}

	engPutVariable(m_ep, "II", II);
	engPutVariable(m_ep, "JJ", JJ);
	engPutVariable(m_ep, "SS", SS);
	engPutVariable(m_ep, "DimM", DimM);
	engPutVariable(m_ep, "DimN", DimN);
	engPutVariable(m_ep, "XX", XX);

	engEvalString(m_ep, "evals = sparse(II, JJ, SS, DimM, DimN) * XX;");

	mxArray *evals = engGetVariable(m_ep, "evals");
	double *pEvals = mxGetPr(evals);
	vY.resize(m_size, 0);
	for (int j = 0; j < m_size; ++j)
		vY[j] = pEvals[j];
	
	mxDestroyArray(II);
	mxDestroyArray(JJ);
	mxDestroyArray(SS);
	mxDestroyArray(DimM);
	mxDestroyArray(DimN);
	mxDestroyArray(XX);
}

void MatlabWrapper::DenseConjugateGradient( const std::vector<std::vector<double> >& A, const std::vector<double>& b, std::vector<double>& sol )
{
	//solve A*sol = B

	assert(A.size() == b.size());

	int sizeA = A.size(), sizeF = A[0].size();
	sol.resize(sizeF);

	mxArray *AA, *BB;
	AA = mxCreateDoubleMatrix(sizeA, sizeF, mxREAL);
	BB = mxCreateDoubleMatrix(sizeA, 1, mxREAL);
	double *pAA = mxGetPr(AA), *pBB = mxGetPr(BB);

	for (int i = 0; i < sizeA; ++i)
	{
		for (int j = 0; j < sizeF; ++j)
		{
			pAA[i + j*sizeA] = A[i][j];
		}
	}

	for (int j = 0; j < sizeA; ++j)
	{
		pBB[j] = b[j];
	}

	engPutVariable(m_ep, "A", AA);
	engPutVariable(m_ep, "b", BB);

	engEvalString(m_ep, "evals = cgls(A, b);");

	mxArray *SS = engGetVariable(m_ep, "evals");
	double *pSS = mxGetPr(SS);
	for (int j = 0; j < sizeF; ++j)
	{
		sol[j] = pSS[j];
	}

	mxDestroyArray(AA);
	mxDestroyArray(BB);
}
