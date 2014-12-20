#include "DifferentialMeshProcessor.h"
#include <fstream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <functional>
#include <ppl.h>
#include <amp.h>
#include <amp_math.h>
#include <mkl.h>
#include <ZGeom/util.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include "global.h"
#include "SpectralGeometry.h"

using namespace std;
using ZGeom::MatlabEngineWrapper;
using ZGeom::VectorPointwiseProduct;
using ZGeom::VectorDotProduct;

std::function<double(double)> transferScaling1 = [](double lambda) { return std::exp(-std::pow(lambda, 2)); };
std::function<double(double, double)> transferHeatKernel = [](double lambda, double t) { return std::exp(-lambda*t); };
std::function<double(double, double)> transferMHW = [](double lambda, double t) { return lambda * std::exp(-lambda*t); };


double transferFunc2(double lambda, double t)
{
	double coeff = lambda * t;
	return coeff * std::exp(-coeff);
}

double heatKernelTransferFunc( double lambda, double t )
{
	return std::exp(-lambda * t);
}

double mhwTransferFunc1( double lambda, double t )
{
	return lambda * std::exp(-lambda * t);
}

DifferentialMeshProcessor::DifferentialMeshProcessor()
{
	mMesh = nullptr;
	mRefVert = 0;
	mActiveHandle = -1;	
}

DifferentialMeshProcessor::DifferentialMeshProcessor(CMesh* tm, CMesh* originalMesh)
{
	mMesh = nullptr;
	mRefVert = 0;
	mActiveHandle = -1;

	init_lite(tm, originalMesh);
}

DifferentialMeshProcessor::~DifferentialMeshProcessor()
{
	std::cout << "DifferentialMeshProcessor destroyed!" << std::endl;
}

void DifferentialMeshProcessor::init(CMesh* tm)
{
	mMesh = tm;
	mOriMesh = tm;
	mRefVert = g_configMgr.getConfigValueInt("INITIAL_REF_POINT");
	mRefPos = mMesh->getVertex(mRefVert)->pos();
}

void DifferentialMeshProcessor::init_lite( CMesh* tm, CMesh* originalMesh )
{
	mMesh = tm;
	mOriMesh = originalMesh;
	mRefVert = 0;
    mRefPos = mMesh->getVertex(mRefVert)->pos();
}

void DifferentialMeshProcessor::constructLaplacian( LaplacianType laplacianType /*= CotFormula*/ )
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


void DifferentialMeshProcessor::decomposeLaplacian( int nEigFunc, LaplacianType laplacianType /*= CotFormula*/ )
{
	ZGeom::logic_assert(hasLaplacian(laplacianType), "laplacian is not available for decomposition");	
	mMeshLaplacians[laplacianType].meshEigenDecompose(nEigFunc, &g_engineWrapper, mMHBs[laplacianType]);
}

void DifferentialMeshProcessor::loadMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	mMHBs[laplacianType].load(path);
}

void DifferentialMeshProcessor::saveMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	mMHBs[laplacianType].save(path);	
}

std::string DifferentialMeshProcessor::generateMHBPath( const std::string& prefix, LaplacianType laplacianType )
{
	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
	std::string pathMHB = prefix + mMesh->getMeshName() + ".mhb." + s_idx;
	return pathMHB;
}

bool DifferentialMeshProcessor::isMHBCacheValid( const std::string& pathMHB, int eigenCount )
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

void DifferentialMeshProcessor::addNewHandle( int hIdx )
{
	auto iter = mHandles.find(hIdx);
	if (iter != mHandles.end()) mHandles.erase(iter);
	else mHandles.insert(std::make_pair(hIdx, mMesh->getVertex(hIdx)->pos()));	 
}

void DifferentialMeshProcessor::calKernelSignature( double scale, KernelType kernelType, std::vector<double>& values ) const
{
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	const int vertCount = mMesh->vertCount();
	values.resize(vertCount);

	TransferFunc pTF = &heatKernelTransferFunc;	// default as heat kernel

	switch(kernelType)
	{
	case HEAT_KERNEL:
		pTF = &heatKernelTransferFunc;
		break;
	case MHW_KERNEL:
		pTF = &mhwTransferFunc1;
		break;
	case SGW_KERNEL:
		pTF = &transferFunc2;
		break;
	}

	Concurrency::parallel_for(0, vertCount, [&](int i)	{
		double sum = 0;
		for (int k = 0; k < mhb.eigVecCount(); ++k)
		{
			sum += (*pTF)(mhb.getEigVal(k), scale) * mhb.getEigVec(k)[i] * mhb.getEigVec(k)[i];
		}
		values[i] = sum;
	});
}

void DifferentialMeshProcessor::computeKernelSignatureFeatures( const std::vector<double>& timescales, KernelType kernelType )
{
	MeshFeatureList mfl;
	int nScales = timescales.size();

	for (int s = 0; s < nScales; ++s) {
		vector<double> vSig;
		vector<int> vFeatures;
		calKernelSignature(timescales[s], kernelType, vSig);
		mMesh->extractExtrema(vSig, 2, 1e-5, vFeatures);

		for (vector<int>::iterator iter = vFeatures.begin(); iter != vFeatures.end(); ++iter) {
			mfl.addFeature(new MeshFeature(*iter, s));			
		}
	}

	std::string featureStr = StrAttrFeatureUnnamed;

	switch(kernelType)
	{
	case HEAT_KERNEL:
		featureStr = StrAttrFeatureHKS;
		break;
	case MHW_KERNEL:
		featureStr = StrAttrFeatureMHWS;
		break;
	}

	mMesh->addAttrMeshFeatures(mfl, featureStr);
}

double DifferentialMeshProcessor::calHK( int v1, int v2, double timescale ) const
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

double DifferentialMeshProcessor::calMHW( int v1, int v2, double timescale ) const
{
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k)	{
		double lambda = mhb.getEigVal(k);
		const ZGeom::VecNd& phi = mhb.getEigVec(k);
		sum += lambda * std::exp(-lambda * timescale) * phi[v1] * phi[v2];
	}
	return sum;
}

double DifferentialMeshProcessor::calHeatTrace(double timescale) const
{
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k) {
		sum += std::exp(-mhb.getEigVal(k) * timescale);
	}
	return sum;
}

double DifferentialMeshProcessor::calBiharmonic(int v1, int v2) const
{
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k) {
		sum += pow( (mhb.getEigVec(k)[v1] - mhb.getEigVec(k)[v2]) / mhb.getEigVal(k), 2 );
	}
	return std::sqrt(sum);
}

void DifferentialMeshProcessor::calHeat( int vSrc, double tMultiplier, std::vector<double>& vHeat )
{
	const int vertCount = mMesh->vertCount();
	vHeat.resize(vertCount);
	ZGeom::VecNd vInitHeat(vertCount, 0);
	vInitHeat[vSrc] = 1.0;
	ZGeom::VecNd vSolvedHeat(vertCount, 0);

	if (mHeatDiffuseMat.empty()) {
		computeHeatDiffuseMat(tMultiplier);
	}
	
	mHeatDiffuseSolver.solve(1, vInitHeat.c_ptr(), vSolvedHeat.c_ptr());
	std::copy_n(vSolvedHeat.c_ptr(), vertCount, vHeat.begin());
}

void DifferentialMeshProcessor::computeHeatDiffuseMat( double tMultiplier )
{
	const int vertCount = mMesh->vertCount();
	const double t = std::pow(mMesh->getAvgEdgeLength(), 2) * tMultiplier;

	const MeshLaplacian& laplacian = getMeshLaplacian(CotFormula);
	const ZGeom::SparseMatrix<double>& matW = laplacian.getW();
	const ZGeom::SparseMatrix<double>& matLc = laplacian.getLS();	// negative

	ZGeom::addMatMat(matW, matLc, -t, mHeatDiffuseMat);	//A = W - t*Lc
	mHeatDiffuseSolver.initialize(mHeatDiffuseMat, true, true);
}

const ZGeom::EigenSystem& DifferentialMeshProcessor::prepareEigenSystem(const MeshLaplacian& laplaceMat, int eigenCount)
{
	LaplacianType laplaceType = laplaceMat.laplacianType();
	std::string pathMHB = generateMHBPath("cache/", laplaceType);
	if (eigenCount == -1) eigenCount = mMesh->vertCount() - 1;

	int useCache = gSettings.LOAD_MHB_CACHE;
	if (useCache != 0 && isMHBCacheValid(pathMHB, eigenCount)) {
		mMHBs[laplaceType].load(pathMHB);
	}
	else {
		laplaceMat.meshEigenDecompose(eigenCount, &g_engineWrapper, mMHBs[laplaceType]);
		mMHBs[laplaceType].save(pathMHB);
	}
	
	return mMHBs[laplaceType];
}
