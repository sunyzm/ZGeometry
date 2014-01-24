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

		m_ep->addArray(&vII[0], ns, 1, false, "II");
		m_ep->addArray(&vJJ[0], ns, 1, false, "JJ");
		m_ep->addArray(&vSS[0], ns, 1, false, "SS");
		m_ep->addArray(&vWeights[0], nOrder, 1, false, "AA");
		m_ep->addArray(&numv, 1, 1, false, "Numv");

		m_ep->eval("[evecs,evals] = genspeigs(II,JJ,SS,AA,Numv);");
		
		double *evec = m_ep->getDblVariablePtr("evecs");
		double *eval = m_ep->getDblVariablePtr("evals");

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
	}

	void EigenCompute::solveStdSym( const SparseMatrix<double>& mLs, int nEig, EigenSystem& eigSys )
	{
		std::vector<int> vII, vJJ;
		std::vector<double> vSS;
		mLs.convertToCOO(vII, vJJ, vSS, ZGeom::MAT_FULL);
		int nOrder = mLs.rowCount();
		int ns = (int)vII.size();
		double numv = nEig;
		double order = (double)nOrder;
		
		m_ep->addArray(&vII[0], ns, 1, false, "II");
		m_ep->addArray(&vJJ[0], ns, 1, false, "JJ");
		m_ep->addArray(&vSS[0], ns, 1, false, "SS");
		m_ep->addArray(&order, 1, 1, false, "order");
		m_ep->addArray(&numv, 1, 1, false, "Numv");

		m_ep->eval("[evecs,evals] = stdspeigs(II,JJ,SS,order,Numv);");

		double *evec = m_ep->getDblVariablePtr("evecs");
		double *eval = m_ep->getDblVariablePtr("evals");

		eigSys.setSize(nOrder, nEig);
		for (int i = 0; i < nEig; ++i) {
			eigSys.mEigVals[i] = std::fabs(eval[i]);
			std::copy_n(evec + i*nOrder, nOrder, eigSys.getEigVec(i).c_ptr());
		}

		m_ep->removeVariable("II");
		m_ep->removeVariable("JJ");
		m_ep->removeVariable("SS");
		m_ep->removeVariable("order");
		m_ep->removeVariable("Numv");
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