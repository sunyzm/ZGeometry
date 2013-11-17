#ifndef ZGEOM_APPROXIMATION_H
#define ZGEOM_APPROXIMATION_H

#include <vector>
#include <tuple>
#include <ZUtil/ZUtil.h>
#include "VecN.h"

namespace ZGeom
{

class  ApproxItem
{
public:
	ApproxItem() : mRes(0), mBasisIdx(-1), mCoeff(0) {}
	ApproxItem(double r, int i, double c) : mRes(r), mBasisIdx(i), mCoeff(c) {}

	double res() const { return mRes; }
	double& res() { return mRes; }
	int index() const { return mBasisIdx; }
	int& index() { return mBasisIdx; }
	double coeff() const { return mCoeff; }
	double& coeff() { return mCoeff; }

private:
	double mRes;
	int mBasisIdx;
	double mCoeff;
};

class ApproximationBasis
{
	friend class FunctionApproximation;

public:


private:
	std::vector<VecNd> mBasis; 
};

class FunctionApproximation
{
public:
	FunctionApproximation() : mpBasis(NULL) {}

	void clear() { mApproxItems.clear(); mpBasis = NULL; }
	void addItem(double r, int i, double c)
	{
		mApproxItems.push_back(ApproxItem(r, i, c));
	}

	const std::vector<ApproxItem>& getApproxItems() const { return mApproxItems; }
	ApproxItem& operator [] (int i) { return mApproxItems[i]; }
	ApproxItem& back() { return mApproxItems.back(); }
	void resize(int n) { mApproxItems.resize(n); }
	int size() const { return (int)mApproxItems.size(); }

private:
	ApproximationBasis *mpBasis;
	std::vector<ApproxItem> mApproxItems;
};

} // end of namespace

#endif