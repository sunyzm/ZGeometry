#ifndef ZGEOM_SPARSE_SYMMETRIC_SOLVER_H
#define ZGEOM_SPARSE_SYMMETRIC_SOLVER_H
#include <ctime>
#include <stdexcept>
#include <mkl.h>
#include <mkl_service.h>
#include <boost/type_traits/is_floating_point.hpp>
#include "VecN.h"

namespace ZGeom
{
	//
	// solver for positive-definite symmetric (PDS) or indefinite symmetric linear system CX = U.
	// (mc,mjc,mic) is the CSR representation of C;
	// U and X should have the size of ne*nrhs.
	// ne: number of equations; nrhs: number of right-hand sides
	//
	template<typename T>
	class SparseSymSolver
	{
	public:
		static const MKL_INT PHASE_FACTORIZATION = 12;
		static const MKL_INT PHASE_SOLVE         = 33;
		static const MKL_INT PHASE_RELEASE       = -1;

		SparseSymSolver() : mInitialized(false) {}
		SparseSymSolver(const SparseMatrix<T>& mat, bool positive_definite, bool verbose = false) { initialize(mat, positive_definite, verbose); }
		SparseSymSolver(MKL_INT ne, T* c, MKL_INT* ic, MKL_INT *jc, bool positive_definite = true, bool verbose = false) { initialize(ne, c, ic, jc, positive_definite, verbose); }
		~SparseSymSolver() { clear(); }

		uint equationCount() const { return mEquationCount; }
		void initialize(const SparseMatrix<T>& mat, bool positive_definite, bool verbose = false);
		void initialize(MKL_INT ne, T* c, MKL_INT* ic, MKL_INT *jc, bool positive_definite = true, bool verbose = false);
		int solve(MKL_INT rhsCount, T* U, T* X);    // U: right-hand side; X: solution
		int solve(const VecN<T>& U, VecN<T>& X);

	protected:
		SparseMatrixCSR<T, MKL_INT> mCsrMat;

	private:
		void clear();

		bool	 mInitialized;
		bool     mIsSingle;
		bool     mVerbose;
		T*       mc;
		MKL_INT* mic;
		MKL_INT* mjc;
		MKL_INT  mEquationCount;  // number of equations
		MKL_INT  mIparm[64];      // parameter array for PARDISO
		MKL_INT  mMtype;          // real and symmetric (semi-)positive definite
		MKL_INT  mMaxfct;         // maximum number of numerical factorization
		MKL_INT  mMnum;           // which matrix to factorize
		MKL_INT  mMsglvl;         // print statistical information in file
		void*    mPt[64];        // pointer to solver internal data
	};

	template<typename T>
	void ZGeom::SparseSymSolver<T>::initialize( const SparseMatrix<T>& mat, bool positive_definite, bool verbose /*= false*/ )
	{
		mat.convertToCSR(mCsrMat, MAT_UPPER);
		this->initialize(mCsrMat.rowCount(), mCsrMat.nzVal(), mCsrMat.rowPtr(), mCsrMat.colIdx(), positive_definite, verbose);
	}

	template<typename T>
	inline void SparseSymSolver<T>::initialize(MKL_INT ne, T *c, MKL_INT *ic, MKL_INT *jc, bool positive_definite, bool verbose)
	{
		assert(boost::is_floating_point<T>::value);
		mIsSingle = (sizeof(T) == sizeof(float));
		mVerbose = verbose;

		// set parameters of PARDISO
		//
		mc = c; mic = ic; mjc = jc;
		mEquationCount  = ne;
		mMtype          = (positive_definite? 2 : -2);
		mMaxfct         = 1;
		mMnum           = 1;
		mMsglvl         = 0; //(verbose ? 1 : 0);

		for (int i = 0; i < 64; ++i) {
			mIparm[i] = 0;
			mPt[i] = NULL;
		}

		mIparm[0]  = 1;                      // not use default parameters
		//mIparm[1]  = 3;                    // parallel nested disection in
		//mIparm[7]  = 2;                    // max numbers of iterative refinement steps
		//mIparm[23] = 1;                    // parallel factorization control
		mIparm[24] = 0;                      // parallel solve
		mIparm[27] = (mIsSingle ? 1 : 0);    // single or double precision data

		MKL_INT idum;
		T       ddum;
		MKL_INT error = 0;              // error flag

		if (verbose) std::cout << "MKL max thread num: " << mkl_get_max_threads() << ", ";
		if (verbose) std::cout << "PARDISO max thread num: " << mkl_domain_get_max_threads(MKL_DOMAIN_PARDISO) << std::endl;

		// Reordering and symbolic factorization
		//
		MKL_INT phase = PHASE_FACTORIZATION;

		PARDISO(mPt, &mMaxfct, &mMnum, &mMtype, &phase, &mEquationCount, mc, mic, mjc,
			&idum, &idum, mIparm, &mMsglvl, &ddum, &ddum, &error);
		if (error != 0) {
			throw std::runtime_error("PARDISO error during factorization!");
		}

		mInitialized = true;	//important for clearing up
	}

	template<typename T>
	inline int SparseSymSolver<T>::solve(MKL_INT rhsCount, T* U, T* X)
	{
		MKL_INT error = 0;
		MKL_INT idum;

		// Solving linear equations
		//
		MKL_INT phase = PHASE_SOLVE;

		PARDISO(mPt, &mMaxfct, &mMnum, &mMtype, &phase, &mEquationCount, mc, mic, mjc,
			&idum, &rhsCount, mIparm, &mMsglvl, U, X, &error);
		if (error != 0) {        
			throw std::runtime_error("PARDISO error during solution!");
		}

		return 0;
	}
	
	template<typename T>
	inline int SparseSymSolver<T>::solve( const VecN<T>& U, VecN<T>& X )
	{
		assert(U.size() == mEquationCount);
		X.resize(mEquationCount);
		return solve(1, U.c_ptr(), X.c_ptr());
	}


	template<typename T>
	inline void SparseSymSolver<T>::clear()
	{
		if (!mInitialized) return;

		MKL_INT idum;
		T       ddum;
		MKL_INT error = 0;              // error flag
		MKL_INT phase = PHASE_RELEASE;
		PARDISO(mPt, &mMaxfct, &mMnum, &mMtype, &phase, &mEquationCount, mc, mic, mjc, &idum, &idum, mIparm, &mMsglvl, &ddum, &ddum, &error);
	}

} //end of namespace ZGeom

#endif