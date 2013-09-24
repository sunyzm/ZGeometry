#include "EigenCompute.h"

namespace ZGeom
{    
	void EigenCompute::solveGenSym( const SparseMatrix<double>& mLs, const SparseMatrix<double>& mW, int nEig, EigenSystem& eigSys )
	{
		std::vector<int> vII, vJJ;
		std::vector<double> vSS, vWeights;
		mLs.convertToCOO(vII, vJJ, vSS, ZGeom::MAT_FULL);
		mW.getDiagonal(vWeights);
		int nOrder = mLs.rowCount();
		int ns = (int) vII.size();
		double numv = nEig;

		m_ep->addVariable(&vWeights[0], nOrder, 1, false, "AA");
		m_ep->addVariable(&vII[0], ns, 1, false, "II");
		m_ep->addVariable(&vJJ[0], ns, 1, false, "JJ");
		m_ep->addVariable(&vSS[0], ns, 1, false, "SS");
		m_ep->addVariable(&numv, 1, 1, false, "Numv");

		m_ep->eval("[evecs,evals] = hspeigs(II,JJ,SS,AA,Numv);");
		
		double *evec = m_ep->getVariablePtr("evecs");
		double *eval = m_ep->getVariablePtr("evals");

		eigSys.setSize(nOrder, nEig);
		for (int i = 0; i < nEig; ++i) {
			eigSys.mEigVals[i] = std::fabs(eval[i]);
			std::copy_n(evec + i*nOrder, nOrder, eigSys.getEigVec(i).c_ptr());
		}

		m_ep->removeVariable("AA");
		m_ep->removeVariable("II");
		m_ep->removeVariable("JJ");
		m_ep->removeVariable("SS");
		m_ep->removeVariable("Numv");

#if 0
		mxArray *II, *JJ, *SS, *AA, *NUMV;
		mxArray *evecs, *evals;

		AA = mxCreateDoubleMatrix(nOrder, 1, mxREAL);
		double *aa = mxGetPr(AA);
		std::copy(vWeights.begin(), vWeights.end(), aa);

		NUMV = mxCreateDoubleMatrix(1, 1, mxREAL);
		double *numv = mxGetPr(NUMV);	
		numv[0] = nEig;			// number of eigen vectors to be computed

		int ns = (int) vII.size();
		II = mxCreateDoubleMatrix(ns, 1, mxREAL);
		JJ = mxCreateDoubleMatrix(ns, 1, mxREAL);
		SS = mxCreateDoubleMatrix(ns, 1, mxREAL);
		double *ii = mxGetPr(II);
		double *jj = mxGetPr(JJ);
		double *ss = mxGetPr(SS);
		std::copy(vII.begin(), vII.end(), ii);
		std::copy(vJJ.begin(), vJJ.end(), jj);
		std::copy(vSS.begin(), vSS.end(), ss);

		m_ep->putVariable("II", II);
		m_ep->putVariable("JJ", JJ);
		m_ep->putVariable("SS", SS);
		m_ep->putVariable("AA", AA);
		m_ep->putVariable("Numv", NUMV);

		m_ep->eval("[evecs,evals] = hspeigs(II,JJ,SS,AA,Numv);");
		
		evecs = m_ep->getVariable("evecs");
		evals = m_ep->getVariable("evals");
		double *evec = mxGetPr(evecs);				//eigenvectors
		double *eval = mxGetPr(evals);				//eigenvalues

		eigSys.setSize(nOrder, nEig);

		for (int i = 0; i < nEig; ++i) {
			eigSys.mEigVals[i] = std::fabs(eval[i]);
			std::copy(evec + i*nOrder, evec + (i+1)*nOrder, eigSys.mEigVecs[i].c_ptr());
		}

		mxDestroyArray(II);
		mxDestroyArray(JJ);
		mxDestroyArray(SS);
		mxDestroyArray(AA);
		mxDestroyArray(NUMV);
#endif
	}

	double computeHeatKernel( const EigenSystem& eigSys, uint x, uint y, double t )
	{
		double sum(0);
		for (uint k = 0; k < eigSys.eigVecCount(); ++k) {
			double lambda = eigSys.getEigVal(k);
			double* phi = eigSys.getEigVec(k).c_ptr();
			sum += std::exp(- lambda * t) * phi[x] * phi[y];
		}
		return sum;
	}



}//end of namespace