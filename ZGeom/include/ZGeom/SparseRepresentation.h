#ifndef ZGEOM_DICTIONARY_H
#define ZGEOM_DICTIONARY_H
#include <algorithm>
#include <functional>
#include "VecN.h"

namespace ZGeom {

class Dictionary
{
public:
	Dictionary() : mDim(0) {}
	Dictionary(int m) : mDim(m) {}
    Dictionary(const Dictionary& dict) = default;
    Dictionary& operator=(const Dictionary& dict) = default;
    Dictionary(Dictionary&& dict)
    {
        mAtoms = std::move(dict.mAtoms);
        mDim = dict.mDim;
        dict.mDim = 0;
    }
    Dictionary& operator=(Dictionary&& dict)
    {
        mAtoms = std::move(dict.mAtoms);
        mDim = dict.mDim;
        dict.mDim = 0;
        return *this;
    }
	const VecNd& operator [] (int i) const { return mAtoms[i]; }
	VecNd& operator[] (int i) { return mAtoms[i]; }

    const std::vector<VecNd>& getAtoms() const { return mAtoms; }
    const std::vector<VecNd>& operator() () const { return mAtoms; }
    
	int atomDim() const { return mDim; }
	int atomCount() const { return (int)mAtoms.size(); }
	int size() const { return (int)mAtoms.size(); }
	void setDimension(int m) { mDim = m; }
	void resize(int ct) { mAtoms.resize(ct); }
	void resize(int ct, int dim)
	{
		mAtoms.resize(ct);
		mDim = dim;
		for (VecNd& v : mAtoms) v.resize(mDim);
	}

	void clear() { mAtoms.clear(); }
	bool empty() const { return mAtoms.size() == 0; }
    double calCoherence() const;
    double calCoherence2() const;
    double mutualCoherence(const Dictionary& dict2) const;

	void expandTo(int N) 
	{
		if (N <= (int)mAtoms.size()) return;
		mAtoms.resize(N);
	}

	void mergeDict(const Dictionary& d2);
	friend void combineDictionary(const Dictionary& d1, const Dictionary& v2, Dictionary& v3);

private:
	std::vector<VecNd> mAtoms;
	int mDim;
};

class SparseCodingItem
{
public:
	friend class SparseCoding;

	SparseCodingItem() : mBasisIdx(-1), mCoeff(0) {}
	SparseCodingItem(int i, double c) : mBasisIdx(i), mCoeff(c) {}

	int index() const { return mBasisIdx; }
	int& index() { return mBasisIdx; }
	double coeff() const { return mCoeff; }
	double& coeff() { return mCoeff; }

private:
	int mBasisIdx;
	double mCoeff;
};

class SparseCoding
{
public:
	SparseCoding() {}
	SparseCoding(const VecNd& vDense, double epsilon = 1e-6) { fromDense(vDense, epsilon); }

	SparseCodingItem& operator [] (int i) { return mApproxItems[i]; }
	const SparseCodingItem& operator [] (int i) const { return mApproxItems[i]; }

	void fromDense(const VecNd& vDense, double epsilon = 1e-6)
	{
		clear();
		for (int i = 0; i < (int)vDense.size(); ++i) {
			if (std::fabs(vDense[i]) > epsilon) {
				addItem(i, vDense[i]);
			}
		}
	}

    VecNd toDense(int totalSize) const
    {
        VecNd vDense(totalSize, 0);
        for (auto sci : mApproxItems) vDense[sci.index()] = sci.coeff();
        return vDense;
    }

	void clear() { mApproxItems.clear(); }
	void addItem(const SparseCodingItem& vi) { mApproxItems.push_back(vi);  }
	void addItem(int i, double c) {	mApproxItems.push_back(SparseCodingItem(i, c));	}		
	const std::vector<SparseCodingItem>& getApproxItems() const { return mApproxItems; }
	SparseCodingItem& back() { return mApproxItems.back(); }
	void resize(int n) { mApproxItems.resize(n); }
	int size() const { return (int)mApproxItems.size(); }
	
	std::vector<int> getAllAtomIndex() const
	{
		std::vector<int> vIdx;
		for (auto ai : mApproxItems) vIdx.push_back(ai.index());
		return vIdx;
	}

	friend std::ostream& operator << (std::ostream &os, const SparseCoding& fa)
	{
		for (auto t : fa.mApproxItems) {
			os << t.index() << ", " << t.coeff() << '\n';
		}
		return os;
	}

	void sortByIndex() 
	{
		std::sort(mApproxItems.begin(), mApproxItems.end(),	
			[](const SparseCodingItem& c1, const SparseCodingItem& c2) {
				return c1.index() < c2.index();
		});
	}

	const SparseCoding& sortByCoeff() 
	{
		std::sort(mApproxItems.begin(), mApproxItems.end(), 
			[](const SparseCodingItem& c1, const SparseCodingItem& c2) {
				return fabs(c1.coeff()) > fabs(c2.coeff());
		});
		return *this;
	}

	SparseCoding extract(std::function<bool(const SparseCodingItem&)> funcSelect) 
	{
		SparseCoding newSC;
		for (auto sci : mApproxItems) {
			if (funcSelect(sci)) newSC.addItem(sci);
		}
		return newSC;
	}

	SparseCoding extractFront(int n) 
	{
		SparseCoding newSC;
		for (auto sci : mApproxItems) {
			if (sci.index() < n) newSC.addItem(sci);
		}
		return newSC;
	}

	SparseCoding extractBack(int n)
	{
		SparseCoding newSC;
		for (auto sci : mApproxItems) {
			if (sci.index() >= n) newSC.addItem(sci);
		}
		return newSC;
	}

private:
	std::vector<SparseCodingItem> mApproxItems;
};

void filterSparseCoding(const SparseCoding& coding1, SparseCoding& coding2, std::function<bool(const SparseCodingItem&)> funcFilter);
void splitSparseCoding(const SparseCoding& coding1, SparseCoding& coding2, SparseCoding& coding3, std::function<bool(const SparseCodingItem&)> funcSplit);
void multiChannelSplitSparseCoding(const std::vector<SparseCoding>& vCodings1, std::vector<SparseCoding>& vCodings2, std::vector<SparseCoding>& vCodings3, std::function<bool(const SparseCodingItem&)> funcSplit);
VecNd dictVecInnerProducts(const Dictionary& dict, const VecNd& vec);
VecNd singleChannelDenseReconstruct(const Dictionary& dict, const VecNd& vCoeff);
std::vector<VecNd> multiChannelDenseReconstruct(const Dictionary& dict, const std::vector<VecNd>& vCoeffs);
VecNd singleChannelSparseReconstruct(const Dictionary& dict, const SparseCoding& approx, int nCode = -1);
VecNd singleChannelSparseReconstruct(const std::vector<VecNd>& vDictAtoms, const SparseCoding& approx, int nCode = -1);
void multiChannelSparseReconstruct(const Dictionary& dict, const std::vector<SparseCoding>& vCodings, std::vector<VecNd>& vReconstructed, int nCode = -1);
void multiChannelSparseReconstruct(const std::vector<VecNd>& vAtoms, const std::vector<SparseCoding>& vApprox, std::vector<VecNd>& vReconstructed, int nCode = -1);
//std::vector<VecNd> SparseReconstructMultiChannel(const Dictionary& dict, const std::vector<SparseCoding>& vCodings);

}
#endif