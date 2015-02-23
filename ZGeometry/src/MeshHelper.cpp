#include "MeshHelper.h"
#include <fstream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <functional>
#include <ppl.h>
#include <mkl.h>
#include <ZGeom/util.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include "differential_geometry.h"

using namespace std;
using ZGeom::MatlabEngineWrapper;
using ZGeom::VectorPointwiseProduct;

MeshHelper::MeshHelper()
{
	mMesh = nullptr;
	mRefVert = 0;
	mActiveHandle = -1;	
}

MeshHelper::MeshHelper(CMesh* tm, CMesh* originalMesh)
{
	mMesh = nullptr;
	mRefVert = 0;
	mActiveHandle = -1;

	init_lite(tm, originalMesh);
}

void MeshHelper::init(CMesh* tm)
{
	mMesh = tm;
	mOriMesh = tm;
	mRefVert = g_configMgr.getConfigValueInt("INITIAL_REF_POINT");
	mRefPos = mMesh->getVertex(mRefVert)->pos();
}

void MeshHelper::init_lite( CMesh* tm, CMesh* originalMesh )
{
	mMesh = tm;
	mOriMesh = originalMesh;
	mRefVert = 0;
    mRefPos = mMesh->getVertex(mRefVert)->pos();
}

void MeshHelper::constructLaplacian( LaplacianType laplacianType /*= CotFormula*/ )
{
	if (hasLaplacian(laplacianType)) return;
	
	MeshLaplacian& laplacian = mMeshLaplacians[laplacianType];
	switch(laplacianType)
	{
	case Tutte:
	case Umbrella:
	case NormalizedUmbrella:
	case CotFormula:
	case SymCot:
	case Anisotropic1:
	case Anisotropic2:
		(laplacian.*(laplacian.getConstructFunc(laplacianType)))(mMesh);
		break;

	case IsoApproximate:
		laplacian.constructIsoApprox(mMesh);                
		break;

	default: throw std::logic_error("Unrecognized Laplacian type");
	}       
}


void MeshHelper::decomposeLaplacian( int nEigFunc, LaplacianType laplacianType /*= CotFormula*/ )
{
	ZGeom::logic_assert(hasLaplacian(laplacianType), "laplacian is not available for decomposition");	
	mMeshLaplacians[laplacianType].meshEigenDecompose(nEigFunc, &g_engineWrapper, mMHBs[laplacianType]);
}

void MeshHelper::loadMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	mMHBs[laplacianType].load(path);
}

void MeshHelper::saveMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	mMHBs[laplacianType].save(path);	
}

std::string MeshHelper::generateMHBPath( const std::string& prefix, LaplacianType laplacianType )
{
	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = prefix + mMesh->getMeshName() + ".mhb." + s_idx;
	return pathMHB;
}

bool MeshHelper::isMHBCacheValid( const std::string& pathMHB, int eigenCount )
{
	if (!fileExist(pathMHB)) return false;

	int nEig, nSize;
	ifstream ifs(pathMHB.c_str(), ios::binary);
	ifs.read((char*)&nEig, sizeof(int));
	ifs.read((char*)&nSize, sizeof(int));
	ifs.close();

	if (nEig != eigenCount || nSize != mMesh->vertCount()) return false;

	return true;
}

void MeshHelper::addNewHandle( int hIdx )
{
	auto iter = mHandles.find(hIdx);
	if (iter != mHandles.end()) mHandles.erase(iter);
	else mHandles.insert(std::make_pair(hIdx, mMesh->getVertex(hIdx)->pos()));	 
}

double MeshHelper::calHK( int v1, int v2, double timescale ) const
{
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k)	{
		double lambda = mhb.getEigVal(k);
		const ZGeom::VecNd& phi = mhb.getEigVec(k);
		sum += std::exp(-lambda * timescale) * phi[v1] * phi[v2];
	}
	return sum;
}

void MeshHelper::calHeat( int vSrc, double tMultiplier, std::vector<double>& vHeat )
{
	const int vertCount = mMesh->vertCount();
	vHeat.resize(vertCount);
	ZGeom::VecNd vInitHeat(vertCount, 0);
	vInitHeat[vSrc] = 1.0;
	ZGeom::VecNd vSolvedHeat(vertCount, 0);

	if (mHeatDiffuseMat.empty()) computeHeatDiffuseMat(tMultiplier);
	
	mHeatDiffuseSolver.solve(1, vInitHeat.c_ptr(), vSolvedHeat.c_ptr());
	std::copy_n(vSolvedHeat.c_ptr(), vertCount, vHeat.begin());
}

void MeshHelper::computeHeatDiffuseMat( double tMultiplier )
{
	const int vertCount = mMesh->vertCount();
	const double t = std::pow(mMesh->getAvgEdgeLength(), 2) * tMultiplier;

	const MeshLaplacian& laplacian = getMeshLaplacian(CotFormula);
	const ZGeom::SparseMatrix<double>& matW = laplacian.getW();
	const ZGeom::SparseMatrix<double>& matLc = laplacian.getLS();	// negative

	ZGeom::addMatMat(matW, matLc, -t, mHeatDiffuseMat);	//A = W - t*Lc
	mHeatDiffuseSolver.initialize(mHeatDiffuseMat, true, true);
}

const ZGeom::EigenSystem& MeshHelper::prepareEigenSystem(const MeshLaplacian& laplaceMat, int eigenCount)
{
	LaplacianType laplaceType = laplaceMat.laplacianType();
    if (!mMHBs[laplaceType].empty()) return mMHBs[laplaceType];

	std::string pathMHB = generateMHBPath("cache/", laplaceType);
	if (eigenCount == -1) eigenCount = mMesh->vertCount() - 1;

	int useCache = gSettings.LOAD_MHB_CACHE;
	if (useCache != 0 && isMHBCacheValid(pathMHB, eigenCount)) {
		mMHBs[laplaceType].load(pathMHB);
	}
	else {
        std::cout << "==== Do Eigendecomposition ====\n";

		laplaceMat.meshEigenDecompose(eigenCount, &g_engineWrapper, mMHBs[laplaceType]);
		mMHBs[laplaceType].save(pathMHB);
	}
	
	return mMHBs[laplaceType];
}
