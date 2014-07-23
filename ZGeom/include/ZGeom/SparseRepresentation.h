#ifndef ZGEOM_DICTIONARY_H
#define ZGEOM_DICTIONARY_H
#include "VecN.h"

namespace ZGeom {

class Dictionary
{
public:
	Dictionary() : mDim(0) {}
	Dictionary(int m) : mDim(m) {}

	const VecNd& operator [] (int i) const { return mAtoms[i]; }
	VecNd& operator[] (int i) { return mAtoms[i]; }

	int atomDim() const { return mDim; }
	int atomCount() const { return (int)mAtoms.size(); }
	int size() const { return (int)mAtoms.size(); }

	void setDimension(int m) { mDim = m; }
	void resize(int N) { mAtoms.resize(N); }
	void resize(int N, int m)
	{
		mDim = m;
		mAtoms.resize(N);
		for (VecNd& v : mAtoms) v.resize(m);
	}

	const std::vector<VecNd>& getAtoms() const { return mAtoms; }
	void clear() { mAtoms.clear(); }

	void expandTo(int N) 
	{
		if (N <= (int)mAtoms.size()) return;
		mAtoms.resize(N);
	}

	const std::vector<VecNd>& operator() () const { return mAtoms; }

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
	SparseCoding(const VecNd& vDense, double epsilon = 1e-6) 
	{
		fromDense(vDense, epsilon);
	}

	void fromDense(const VecNd& vDense, double epsilon = 1e-6)
	{
		clear();
		for (int i = 0; i < (int)vDense.size(); ++i) {
			if (std::fabs(vDense[i]) >= epsilon) {
				addItem(i, vDense[i]);
			}
		}
	}

	void clear() { mApproxItems.clear(); }

	void addItem(int i, double c)
	{
		mApproxItems.push_back(SparseCodingItem(i, c));
	}
		
	const std::vector<SparseCodingItem>& getApproxItems() const { return mApproxItems; }
	SparseCodingItem& operator [] (int i) { return mApproxItems[i]; }
	const SparseCodingItem& operator [] (int i) const { return mApproxItems[i]; }
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

private:
	std::vector<SparseCodingItem> mApproxItems;
};

VecNd DictVecInnerProducts(const Dictionary& dict, const VecNd& vec);
VecNd DenseReconstructSingleChannel(const Dictionary& dict, const VecNd& vCoeff);
std::vector<VecNd> DenseReconstructMultiChannel(const Dictionary& dict, const std::vector<VecNd>& vCoeffs);
VecNd SparseReconstructSingleChannel(const Dictionary& dict, const SparseCoding& approx);
VecNd SparseReconstructSingleChannel(const std::vector<VecNd>& vAtoms, const SparseCoding& approx);
std::vector<VecNd> SparseReconstructMultiChannel(const Dictionary& dict, const std::vector<SparseCoding>& vCodings);
void ReconstructApproximationMultiChannel(const std::vector<VecNd>& vAtoms, const std::vector<SparseCoding>& vApprox, std::vector<VecNd>& vReconstructed);

}
#endif