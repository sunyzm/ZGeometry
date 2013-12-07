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

class SignalAtom : public VecNd
{
public:
	SignalAtom() : VecNd(), mScale(-1), mPosition(-1) {}
	SignalAtom(const VecNd& v, int scale = -1, int pos = -1) : VecNd(v), mScale(scale), mPosition(pos) {}
	SignalAtom(double* p, uint dim, int scale = -1, int pos = -1) : VecNd(p, dim), mScale(scale), mPosition(pos) {}
	SignalAtom(const std::vector<double>& v, int scale = -1, int pos = -1) : VecNd(v), mScale(scale), mPosition(pos) {}
	
protected:
	int mScale;
	int mPosition;
};

typedef std::vector<SignalAtom*> ApproximationAtoms;

class FunctionApproximation
{
public:
	FunctionApproximation() : mAtoms(NULL) {}

	void addItem(double r, int i, double c)
	{
		mApproxItems.push_back(ApproxItem(r, i, c));
	}

	void clear() { mApproxItems.clear(); mAtoms = NULL; }
	const std::vector<ApproxItem>& getApproxItems() const { return mApproxItems; }
	ApproxItem& operator [] (int i) { return mApproxItems[i]; }
	const ApproxItem& operator [] (int i) const { return mApproxItems[i]; }
	ApproxItem& back() { return mApproxItems.back(); }
	void resize(int n) { mApproxItems.resize(n); }
	int size() const { return (int)mApproxItems.size(); }

	std::vector<int> getAllAtomIndex() const 
	{
		std::vector<int> vIdx;
		for (auto ai : mApproxItems) vIdx.push_back(ai.index());
		return vIdx;
	}

private:
	ApproximationAtoms *mAtoms;
	std::vector<ApproxItem> mApproxItems;
};

} // end of namespace

#endif
