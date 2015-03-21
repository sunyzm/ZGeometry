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
	mRefVert = 0;
	mActiveHandle = -1;	
    currentMeshIdx = -1;
}

 MeshHelper::MeshHelper(MeshHelper && mh)
 {
     mMeshHistory = std::move(mMeshHistory);
     mRefVert = mh.mRefVert;
     mRefPos = mh.mRefPos;
     mHandles = std::move(mh.mHandles);
     mActiveHandle = std::move(mh.mActiveHandle);
 }

void MeshHelper::init(CMesh* tm)
{
    assert(tm != nullptr);
    mMeshHistory.resize(1);
    mMeshHistory[0].reset(tm);
    tm->addAttr<std::string>("original_mesh", CMesh::StrAttrMeshDescription, AR_UNIFORM, AT_STRING);
    currentMeshIdx = 0;
	mRefVert = g_configMgr.getConfigValueInt("INITIAL_REF_POINT");
	mRefPos = getMesh()->vert(mRefVert)->pos();
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
        (laplacian.*(laplacian.getConstructFunc(laplacianType)))(getMesh());
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
    std::string pathMHB = prefix + getMesh()->getMeshName() + ".mhb." + s_idx;
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

    if (nEig != eigenCount || nSize != getMesh()->vertCount()) return false;

	return true;
}

void MeshHelper::addNewHandle( int hIdx )
{
	auto iter = mHandles.find(hIdx);
	if (iter != mHandles.end()) mHandles.erase(iter);
    else mHandles.insert(std::make_pair(hIdx, getMesh()->vert(hIdx)->pos()));
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

const ZGeom::EigenSystem& MeshHelper::prepareEigenSystem(const MeshLaplacian& laplaceMat, int eigenCount)
{
	LaplacianType laplaceType = laplaceMat.laplacianType();
    if (!mMHBs[laplaceType].empty()) return mMHBs[laplaceType];

	std::string pathMHB = generateMHBPath("cache/", laplaceType);
    if (eigenCount == -1) eigenCount = getMesh()->vertCount() - 1;

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

void MeshHelper::revertOriginal()
{
    currentMeshIdx = 0;
    clearMeshRelated();
}

void MeshHelper::addMesh(std::unique_ptr<CMesh> && newMesh, const std::string description)
{
    mMeshHistory.push_back(std::move(newMesh));
    CMesh *mesh = mMeshHistory.back().get();
    mesh->setMeshDescription(description);
    currentMeshIdx = (int)mMeshHistory.size() - 1;
    clearMeshRelated();
}

void MeshHelper::nextMesh()
{
    if (mMeshHistory.size() <= 1) return;
    currentMeshIdx = (currentMeshIdx + 1) % mMeshHistory.size();    

    std::cout << "Switch to: " << mMeshHistory[currentMeshIdx]->getMeshDescription() << std::endl;
    clearMeshRelated();
}

void MeshHelper::prevMesh()
{
    if (mMeshHistory.size() <= 1) return;
    currentMeshIdx = (currentMeshIdx + mMeshHistory.size() - 1) % mMeshHistory.size();

    std::cout << "Switch to: " << mMeshHistory[currentMeshIdx]->getMeshDescription() << std::endl;
    clearMeshRelated();
}

void MeshHelper::clearMeshRelated()
{
    mRefVert = 0;
    mRefPos = getMesh()->vertPos(mRefVert);
    mHandles.clear();
    mActiveHandle = -1;
    for (MeshLaplacian& lm : mMeshLaplacians) lm.clear();
    for (ZGeom::EigenSystem & es : mMHBs) es.clear();
}

CMesh* MeshHelper::getMeshByName(const std::string mesh_descript)
{
    for (int k = 0; k < mMeshHistory.size(); ++k) {
        if (mMeshHistory[k]->getMeshDescription() == mesh_descript)
            return mMeshHistory[k].get();
    }

    return nullptr; 
}

void MeshHelper::removeCurrentMesh()
{
    if (currentMeshIdx == 0) {
        cout << "Original mesh cannot be removed!" << std::endl;
        return;
    }

    mMeshHistory.erase(mMeshHistory.begin() + currentMeshIdx);

    currentMeshIdx = currentMeshIdx % mMeshHistory.size();
    std::cout << "Switch to: " << mMeshHistory[currentMeshIdx]->getMeshDescription() << std::endl;
    clearMeshRelated();
}

