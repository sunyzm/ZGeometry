#ifndef ZGEOM_DICTIONARY_H
#define ZGEOM_DICTIONARY_H
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

	void expandTo(int N) {
		if (N <= (int)mAtoms.size()) return;
		mAtoms.resize(N);
	}

	const std::vector<VecNd>& operator() () const { return mAtoms; }

private:
	std::vector<VecNd> mAtoms;
	int mDim;
};

}
#endif