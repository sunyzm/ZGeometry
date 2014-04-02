#ifndef ZGEOM_APPROXIMATION_H
#define ZGEOM_APPROXIMATION_H
#include <vector>
#include <tuple>
#include <ZUtil/ZUtil.h>
#include "VecN.h"

namespace ZGeom
{

class Dictionary
{
public:
	Dictionary() : mDim(0) {}
	void setDimension(int m) { mDim = m; }
	void clear() { mAtoms.clear(); }
	void resize(int N, int m) 
	{ 
		mDim = m;
		mAtoms.resize(N); 
		for (VecNd& v : mAtoms) v.resize(m);
	}
	void resize(int N) { mAtoms.resize(N); }
	const VecNd& operator [] (int i) const { return mAtoms[i]; }
	VecNd& operator[] (int i) { return mAtoms[i]; }
	int atomCount() const { return (int)mAtoms.size(); }
	const std::vector<VecNd>& getAtoms() const { return mAtoms; }
	void expandTo(int N) {
		if (N <= (int)mAtoms.size()) return;
		mAtoms.resize(N);
	}

private:
	std::vector<VecNd> mAtoms;
	int mDim;
};

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

class FunctionApproximation
{
public:
	FunctionApproximation() {}

	void addItem(double r, int i, double c)
	{
		mApproxItems.push_back(ApproxItem(r, i, c));
	}

	void clear() { mApproxItems.clear(); }
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
	std::vector<ApproxItem> mApproxItems;
};

} // end of namespace

#endif
