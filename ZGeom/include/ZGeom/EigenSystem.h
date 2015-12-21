#ifndef ZGEOM_EIGEN_SYSTEM_H
#define ZGEOM_EIGEN_SYSTEM_H
#include <vector>
#include <functional>
#include "VecN.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"

namespace ZGeom {

class EigenSystem
{
public:
	friend class EigenCompute;

    void clear();
	void setSize(int order, int nev);
	bool empty() const { return mEigVecs.empty(); }
    int eigVecCount() const { return (int)mEigVecs.size(); }
    int eigVecSize() const { return empty() ? 0 : (int)mEigVecs[0].size(); }
	const std::vector<double>& getAllEigVals() const { return mEigVals; }
    const double* getAllEigValsPtr() const { return mEigVals.data(); }
    const std::vector<VecNd>& getAllEigVecs() const { return mEigVecs; }
	double getEigVal(int index) const { return mEigVals[index]; }
	const VecNd& getEigVec(int index) const { return mEigVecs[index]; }
	VecNd& getEigVec(int index) { return mEigVecs[index]; }
    ZGeom::DenseMatrixd toDenseMatrix() const;
	void print(const std::string& file1, const std::string& file2) const;
	void printEigVals(const std::string& file1) const;
	void save(const std::string& file) const;
	void load(const std::string& file);
    void resize(int m);
    bool hasInducingMat() const { return !mInducingMat.empty(); }
    const ZGeom::SparseMatrixd getInducingMat() const { return mInducingMat; }
    bool hasEigenMat() const { return !mEigenMat.empty(); }
    void computeEigenMat();
    const ZGeom::DenseMatrixd& getEigenMat();

    void validate();

public:
	std::vector<double> mEigVals;
	std::vector<VecNd> mEigVecs;

private:    
    ZGeom::SparseMatrixd mInducingMat;
    ZGeom::DenseMatrixd mEigenMat;
};

}	// end of namespace

#endif
