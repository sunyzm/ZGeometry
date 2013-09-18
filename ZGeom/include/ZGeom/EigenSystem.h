#ifndef ZGEOM_EIGEN_SYSTEM_H
#define ZGEOM_EIGEN_SYSTEM_H

#include <vector>
#include <functional>
#include "VecN.h"
#include "MatVecFunctor.h"

namespace ZGeom
{
	class EigenSystem
	{
	public:
		friend class EigenCompute;
		EigenSystem() : mOrder(0), mEvCount(0) {}

		void setSize(int order, int nev);
		bool empty() const { return mEvCount == 0 ; }
		uint eigVecSize() const { return mOrder; }
		uint eigVecCount() const { return mEvCount; }
		void setValues(int index, double eigVal, double eigVec[]);
		void setValues(int index, double eigVal, const std::vector<double>& eigVec);
		const std::vector<double>& getEigVals() const { return mEigVals; }
		double getEigVal(int index) const { return mEigVals[index]; }
		const VecNd& getEigVec(int index) const { return mEigVecs[index]; }
		VecNd& getEigVec(int index) { return mEigVecs[index]; }

		void print(const std::string& file1, const std::string& file2) const;
		void printEigVals(const std::string& file1) const;
		void save(const std::string& file) const;
		void load(const std::string& file);

		void evalError(MatVecFunctor* A, std::vector<double>& vErrors) const;
		void evalGeneralError(MatVecFunctor* A, MatVecFunctor* M, std::vector<double>& vErrors) const;
		void inverseEigVals();

		double heatKernel(uint x, uint y, double t) const;
		double kernel(uint x, uint y, std::function<double(double)> transferFunc) const;
		double kernel(uint x, uint y, double t, std::function<double(double, double)> transferFunc) const;

	protected:
		uint mOrder;
		uint mEvCount;
		std::vector<double> mEigVals;
		std::vector<VecNd> mEigVecs;
	};
}



#endif