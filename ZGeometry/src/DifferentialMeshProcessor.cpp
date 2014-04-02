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
#include <ZUtil/SimpleConfigLoader.h>
#include <ZUtil/zassert.h>
#include <ZUtil/timer.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include "global.h"

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
	mRefPos = mMesh->getVertex(mRefVert)->getPosition();
}

void DifferentialMeshProcessor::init_lite( CMesh* tm, CMesh* originalMesh )
{
	mMesh = tm;
	mOriMesh = originalMesh;
	mRefVert = 0;
	mRefPos = mMesh->getVertex(0)->getPosition();
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
		(laplacian.*(laplacian.getConstructFunc(laplacianType)))(mMesh);
		break;

	case Anisotropic1:
		{
			double para1 = 2 * mMesh->getAvgEdgeLength() * mMesh->getAvgEdgeLength();
			double para2 = para1;
			int nRing = 1;
			laplacian.constructAnisotropic1(mMesh, nRing, para1, para2);
		}
		break;

	case Anisotropic2:
		{
			double para1 = 2 * mMesh->getAvgEdgeLength() * mMesh->getAvgEdgeLength();
			double para2 = mMesh->getAvgEdgeLength() / 2;
			int nRing = 1;
			laplacian.constructAnisotropic2(mMesh, nRing, para1, para2);
		}
		break;

	case IsoApproximate:
		laplacian.constructFromMesh5(mMesh);                
		break;

	default: throw std::logic_error("Unrecognized Laplacian type");
	}       
}


void DifferentialMeshProcessor::decomposeLaplacian( int nEigFunc, LaplacianType laplacianType /*= CotFormula*/ )
{
	ZUtil::logic_assert(hasLaplacian(laplacianType), "laplacian is not available for decomposition");
	if (!g_engineWrapper.isOpened())
		throw std::logic_error("Matlab engine not opened for Laplacian decomposition!");
	
	if (laplacianType == Tutte || laplacianType == CotFormula)
		mMeshLaplacians[laplacianType].decompose(nEigFunc, &g_engineWrapper, mMHBs[laplacianType], true);
	else 
		mMeshLaplacians[laplacianType].decompose(nEigFunc, &g_engineWrapper, mMHBs[laplacianType], false);
}

void DifferentialMeshProcessor::loadMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	mMHBs[laplacianType].load(path);
}

void DifferentialMeshProcessor::saveMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	mMHBs[laplacianType].save(path);	
}

void DifferentialMeshProcessor::computeCurvature( std::vector<double>& vCurvature, int curvatureType /*= 0*/ )
{
	const int vertCount = mMesh->vertCount();
	vCurvature.resize(vertCount);
	if (curvatureType == 0)	vCurvature = mMesh->getMeanCurvature();
	else if (curvatureType == 1) vCurvature = mMesh->getGaussCurvature();
}

void DifferentialMeshProcessor::addNewHandle( int hIdx )
{
	auto iter = mHandles.find(hIdx);
	if (iter != mHandles.end()) mHandles.erase(iter);
	else mHandles.insert(std::make_pair(hIdx, mMesh->getVertex(hIdx)->getPosition()));	 
}

void DifferentialMeshProcessor::computeSGW1( LaplacianType laplacianType /*= CotFormula*/ )
{
	CStopWatch timer;	
	timer.startTimer();
	const ZGeom::EigenSystem& mhb = getMHB(laplacianType);
	computeSGWMat(mhb, 5, mMatAtoms);
	timer.stopTimer("Time to compute SGW1: ");
}

void DifferentialMeshProcessor::computeSGW2(LaplacianType laplacianType /*= CotFormula*/)
{
	CStopWatch timer;	
	timer.startTimer();
	const ZGeom::EigenSystem& mhb = getMHB(laplacianType);
	computeSGWMat2(mhb, 5, mMatAtoms);
	timer.stopTimer("Time to compute SGW2: ");
}

void DifferentialMeshProcessor::computeMixedAtoms1( LaplacianType laplacianType /*= CotFormula*/ )
{
	CStopWatch timer;	
	timer.startTimer();

	const ZGeom::EigenSystem& mhb = getMHB(laplacianType);
	const int vertCount = mhb.eigVecSize();
	const int eigCount = mhb.eigVecCount();

	std::function<double(double)> generator1 = [](double x) {
		if (x < 1) return x*x;
		else if (x <= 2) return (-5. + 11.*x - 6.*x*x + x*x*x);
		else return 4.0/x/x;
	};

	std::cout << "Ev(1) = " << mhb.getEigVal(1) 
			  << ", Ev(" << mhb.eigVecCount()-1 << ") = " << mhb.getEigVal(mhb.eigVecCount()-1)
			  << "\n";

	const int nWaveletScales = 5;
	std::vector<double> vWaveletScales(nWaveletScales);	
	double maxEigVal = mhb.getEigVal(mhb.eigVecCount()-1);
	const double K = 20.0;
	double minEigVal = maxEigVal / K;
	double minT = 1./maxEigVal, maxT = 2./minEigVal;		
	const double tMultiplier = std::pow(maxT/minT, 1.0 / double(nWaveletScales - 1));	
	for (int s = 0; s < nWaveletScales; ++s) 
		vWaveletScales[s] = minT * std::pow(tMultiplier, s);

	double scalingScales[] = {/*0.5, 1.0,*/ 2.0 /*, 4.0, 8.0*/};
	std::vector<double> vScalingScales(std::begin(scalingScales), std::end(scalingScales));
	const int mScalingScales = (int)vScalingScales.size();;
	
	// allocate storage for mMatAtoms
	mMatAtoms.resize(vertCount * (nWaveletScales + mScalingScales), vertCount);

	ZGeom::DenseMatrixd matEigVecs(eigCount, vertCount);	
	const double *pEigVals = &(mhb.getEigVals()[0]);
	double *pEigVec = matEigVecs.raw_ptr();
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(mhb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);

	std::vector<double> vDiag(eigCount);
	//////////////////////////////////////////////////////////////////////////
	// compute SGW with AMP
	for (int s = 0; s < nWaveletScales; ++s) {
		for (int i = 0; i < eigCount; ++i) 
			vDiag[i] = generator1(vWaveletScales[s] * pEigVals[i]);
		double *pResult = mMatAtoms.raw_ptr() + vertCount * vertCount * s;
		ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
	}

	for (int s = 0; s < mScalingScales; ++s) {
		for (int i = 0; i < eigCount; ++i) 
			vDiag[i] = std::exp(-vScalingScales[s] * 8/*7.716*/ * std::pow(pEigVals[i]/minEigVal, 4));
		double *pResult = mMatAtoms.raw_ptr() + vertCount * vertCount * (nWaveletScales + s);
		ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
	}

	timer.stopTimer("Time to compute mixed atoms 1: ");
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

	std::string featureStr = StrFeatureUnnamed;

	switch(kernelType)
	{
	case HEAT_KERNEL:
		featureStr = StrFeatureHKS;
		break;
	case MHW_KERNEL:
		featureStr = StrFeatureMHWS;
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

double DifferentialMeshProcessor::calSGW( int v1, int v2, double timescale, LaplacianType laplacianType ) const
{
	const ZGeom::EigenSystem& mhb = getMHB(laplacianType);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k)	{
		double lambda = mhb.getEigVal(k);
		const ZGeom::VecNd& phi = mhb.getEigVec(k);
		sum += std::pow(lambda * timescale, 2) * std::exp(-std::pow(lambda * timescale, 2)) * phi[v1] * phi[v2];
		//sum += std::exp(-lambda*lambda) * phi[v1] * phi[v2];
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

void DifferentialMeshProcessor::computeSimilarityMap1( int refPoint )
{
	const int vertCount = mMesh->vertCount();
	const CVertex* pvi = mMesh->getVertex(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);
	vFaces.resize(mMesh->faceCount());
	for (int f = 0; f < mMesh->faceCount(); ++f)
		vFaces[f] = f;

	vector<double> vSimilarities;
	vSimilarities.resize(mMesh->vertCount(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mMesh->vertCount(), 0.);

	double hPara1 = std::pow(mMesh->getAvgEdgeLength() * 5, 2);
	double hPara2 = std::pow(mMesh->getAvgEdgeLength(), 2);
	const std::vector<Vector3D>& vVertNormals = mMesh->getVertNormals();

	for (int fi = 0; fi < vFaces.size(); ++fi) {
		const CFace* pfi = mMesh->getFace(vFaces[fi]);
		double face_area = pfi->computeArea();
		for (int k = 0; k < 3; ++k) {
			int vki = pfi->getVertexIndex(k);
//			if (vki == refPoint) continue;
			const CVertex* pvk = pfi->getVertex(k);

			double w1 = 1., w2 = 1.;
//			double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
			w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
//			w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
			w2 = std::exp(-std::pow(dotProduct3D(vVertNormals[refPoint], pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
//			w2 = std::exp((dotProduct3D(pvi->getNormal(), pfi->getNormal()) - 1) / 1.0);

			double svalue = w1 * w2;
			
			vSimilarities[vki] += svalue * face_area;
			vAreas[vki] += face_area;
		}
	}

	std::transform(vSimilarities.begin(), vSimilarities.end(), vAreas.begin(), vSimilarities.begin(), [](double v1, double v2) { return v1 / v2; } );

	mMesh->addAttrVertScalars(vSimilarities, "Simialrity_Map_Signature");
}

void DifferentialMeshProcessor::computeSimilarityMap2( int refPoint )
{
	const int vertCount = mMesh->vertCount();
	const CVertex* pvi = mMesh->getVertex(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
	//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);
	vFaces.resize(mMesh->faceCount());
	for (int f = 0; f < mMesh->faceCount(); ++f) vFaces[f] = f;

	vector<double> vSimilarities;
	vSimilarities.resize(mMesh->vertCount(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mMesh->vertCount(), 0.);

	double hPara1 = std::pow(mMesh->getAvgEdgeLength() * 5, 2);
	double hPara2 = std::pow(mMesh->getAvgEdgeLength(), 2);
	const std::vector<Vector3D>& vVertNormals = mMesh->getVertNormals();
	const std::vector<Vector3D>& vFaceNormals = mMesh->getFaceNormals();
	const int faceNum = mMesh->faceCount();

	for (int fIndex = 0; fIndex < faceNum; ++fIndex) {
		const CFace* pfi = mMesh->getFace(vFaces[fIndex]);
		double face_area = pfi->computeArea();
		for (int k = 0; k < 3; ++k) {
			int vki = pfi->getVertexIndex(k);
			//			if (vki == refPoint) continue;
			const CVertex* pvk = pfi->getVertex(k);

			double w1 = 1., w2 = 1.;
//			double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
			w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
//			w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
//			w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
			w2 = std::exp((dotProduct3D(vVertNormals[refPoint], vFaceNormals[fIndex]) - 1) / 1.0);

			double svalue = w1 * w2;

			vSimilarities[vki] += svalue * face_area;
			vAreas[vki] += face_area;
		}
	}

	std::transform(vSimilarities.begin(), vSimilarities.end(), vAreas.begin(), vSimilarities.begin(), [](double v1, double v2) { return v1 / v2; } );

	mMesh->addAttrVertScalars(vSimilarities, "Simialrity_Map_Signature");
}

void DifferentialMeshProcessor::computeSimilarityMap3( int refPoint )
{
	const int vertCount = mMesh->vertCount();
	const CVertex* pvi = mMesh->getVertex(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
	vFaces.resize(mMesh->faceCount());
	for (int f = 0; f < mMesh->faceCount(); ++f)
		vFaces[f] = f;
//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);

	vector<double> vSimilarities;
	vSimilarities.resize(mMesh->vertCount(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mMesh->vertCount(), 0.);

	double hPara1 = std::pow(mMesh->getAvgEdgeLength() * 5, 2);
//	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);
	double hPara2 = mMesh->getAvgEdgeLength();
	const std::vector<Vector3D>& vVertNormals = mMesh->getVertNormals();

	for (int fi = 0; fi < vFaces.size(); ++fi)
	{
		const CFace* pfi = mMesh->getFace(vFaces[fi]);
		double face_area = pfi->computeArea();
		for (int k = 0; k < 3; ++k)
		{
			int vki = pfi->getVertexIndex(k);
			//			if (vki == refPoint) continue;
			const CVertex* pvk = pfi->getVertex(k);

			double w1 = 1., w2 = 1.;
			//			double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
			w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
			//			w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
			w2 = std::exp(-dotProduct3D(vVertNormals[refPoint], pvk->getPosition() - pvi->getPosition()) / hPara2);
//			w2 = std::exp((dotProduct3D(pvi->getNormal(), pfi->getNormal()) - 1) / 1.0);

			double svalue = w1 * w2;

			vSimilarities[vki] += svalue * face_area;
			vAreas[vki] += face_area;
		}
	}

	std::transform(vSimilarities.begin(), vSimilarities.end(), vAreas.begin(), vSimilarities.begin(), [](double v1, double v2) { return v1 / v2; } );

	mMesh->addAttrVertScalars(vSimilarities, "Simialrity_Map_Signature");
}

void DifferentialMeshProcessor::calHeat( int vSrc, double tMultiplier, std::vector<double>& vHeat )
{
	const int vertCount = mMesh->vertCount();
	vHeat.resize(vertCount);

	if (mHeatDiffuseMat.empty()) computeHeatDiffuseMat(tMultiplier);

	ZGeom::VecNd vInitHeat(vertCount, 0);
	vInitHeat[vSrc] = 1.0;
	ZGeom::VecNd vSolvedHeat(vertCount, 0);

	mHeatDiffuseSolver.solve(1, vInitHeat.c_ptr(), vSolvedHeat.c_ptr());

// 	mpEngineWrapper->addColVec(vInitHeat, "vInit");
// 	mpEngineWrapper->addSparseMat(mHeatDiffuseSolver, "matA");
// 	mpEngineWrapper->addColVec(vSolvedHeat, "vSolved");

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

void DifferentialMeshProcessor::computeHeatKernelMat( double t, ZGeom::DenseMatrix<double>& hkmat )
{
	const int vertCount = mMesh->vertCount();
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	const int eigCount = mhb.eigVecCount();
	hkmat.resize(vertCount, vertCount);
	double *pResult = hkmat.raw_ptr();

	double *pEigVec = new double[vertCount*eigCount];
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(mhb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);
	
	const double *pEigVals = &(mhb.getEigVals()[0]);
	std::vector<double> vDiag(eigCount);
	for (int i = 0; i < eigCount; ++i) vDiag[i] = std::exp(-pEigVals[i]*t);

	CStopWatch timer;
	timer.startTimer();
	Concurrency::parallel_for(0, vertCount, [&](int i) {
		for (int j = 0; j <= i; ++j) {
			double sum(0);
			for (int k = 0; k < eigCount; ++k) {
				sum += vDiag[k] * pEigVec[k*vertCount + i] * pEigVec[k*vertCount + j];
			}
			pResult[i*vertCount + j] = pResult[j*vertCount + i] = sum;
		}
	});
	timer.stopTimer("HK computation time with PPL: ");

	delete []pEigVec;
}

void DifferentialMeshProcessor::computeHeatKernelMat_AMP( double t, ZGeom::DenseMatrix<double>& hkmat )
{
	const int vertCount = mMesh->vertCount();
	const ZGeom::EigenSystem& mhb = getMHB(CotFormula);
	const int eigCount = mhb.eigVecCount();
	hkmat.resize(vertCount, vertCount);
	double *pResult = hkmat.raw_ptr();

	double *pEigVec = new double[eigCount*vertCount];
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(mhb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);

	const double *pEigVals = &(mhb.getEigVals()[0]);
	std::vector<double> vDiag(eigCount);
	for (int i = 0; i < eigCount; ++i) vDiag[i] = std::exp(-pEigVals[i]*t);

	CStopWatch timer;
	timer.startTimer();
	ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
	timer.stopTimer("HK computation time with AMP: ");
	delete []pEigVec;
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
	if (!ZUtil::fileExist(pathMHB)) return false;

	int nEig, nSize;
	ifstream ifs(pathMHB.c_str(), ios::binary);
	ifs.read((char*)&nEig, sizeof(int));
	ifs.read((char*)&nSize, sizeof(int));
	ifs.close();

	if (nEig != eigenCount || nSize != mMesh->vertCount()) return false;

	return true;
}

void DifferentialMeshProcessor::computeGeometricLaplacianCoordinate( const CMesh& mesh, const MeshCoordinates& eCoord, MeshCoordinates& lCoord )
{
	using ZGeom::Vec3d;
	assert(mesh.vertCount() == eCoord.size());
	int vertCount = mesh.vertCount();
	lCoord.resize(vertCount);

	for (int i = 0; i < vertCount; ++i)
	{
		Vec3d vi = eCoord.getVertCoordinate(i);
		std::vector<int> neighbors = mesh.getVertNeighborVerts(i, 1, false);
		Vec3d vavg(0,0,0);
		double weightSum(0);
		for (int j : neighbors) {
			Vec3d vj = eCoord.getVertCoordinate(j);
			double weight = 1.0 / (vi - vj).length();
			vavg += vj * weight;
			weightSum += weight;
		}
		vavg /= weightSum;
		lCoord[i] = vi - vavg;
	}
}

void DifferentialMeshProcessor::computeSGWMat( const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW )
{
	const int vertCount = mhb.eigVecSize();
	const int eigCount = mhb.eigVecCount();
	const int nWaveletScales = waveletScaleNum;
	const int nScalingScales = 1;
	const int totalAtomCount = vertCount*(nWaveletScales+1);
	matSGW.resize(totalAtomCount, vertCount);

	double maxEigVal = mhb.getEigVal(mhb.eigVecCount()-1);
	const double K = 20.0;
	double minEigVal = maxEigVal / K;
	double minT = 2./maxEigVal, maxT = 2./minEigVal;		
	const double tMultiplier = std::pow(maxT/minT, 1.0 / double(nWaveletScales - 1));
	const double gamma = 1.3849;

	std::function<double(double)> genG = [](double x) {
		if (x < 1) return x*x;
		else if (x <= 2) return (-5. + 11.*x - 6.*x*x + x*x*x);
		else return 4.0/x/x;
	};
	std::function<double(double)> genH = [=](double x) {
		return gamma * std::exp(-std::pow(x/(0.6*minEigVal), 4));
	};

	std::vector<double> vWaveletScales(nWaveletScales);	
	for (int s = 0; s < nWaveletScales; ++s)
		vWaveletScales[s] = minT * std::pow(tMultiplier, s);

	ZGeom::DenseMatrixd matEigVecs(eigCount, vertCount);	
	const double *pEigVals = &(mhb.getEigVals()[0]);
	double *pEigVec = matEigVecs.raw_ptr();
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(mhb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);

	//////////////////////////////////////////////////////////////////////////
	// compute SGW with AMP
	std::vector<double> vDiag(eigCount);
	for (int s = 0; s < nWaveletScales; ++s) 
	{
		for (int i = 0; i < eigCount; ++i) 
			vDiag[i] = genG(vWaveletScales[s] * pEigVals[i]);
		double *pResult = matSGW.raw_ptr() + vertCount * vertCount * s;
		ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
	}

	{
		for (int i = 0; i < eigCount; ++i) 
			vDiag[i] = genH(pEigVals[i]);
		double *pResult = matSGW.raw_ptr() + vertCount * vertCount * nWaveletScales;
		ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
	}
}

void DifferentialMeshProcessor::computeSGWMat2( const ZGeom::EigenSystem& mhb, int nWaveletScales, ZGeom::DenseMatrixd& matSGW )
{
	const int vertCount = mhb.eigVecSize();
	const int eigCount = mhb.eigVecCount();

	std::function<double(double)> generator1 = [](double x) { return x * std::exp(-x); };
	std::function<double(double)> generator2 = [](double x) { return std::exp(-x-1); };

	const int nScales = 5;
	const int totalAtomCount = vertCount*(nWaveletScales+1);
	double maxEigVal = mhb.getEigVal(mhb.eigVecCount()-1);
	double minEigVal = maxEigVal / 20.;//mhb.getEigVal(1);
	double maxT = 1.0 / minEigVal;
	double minT = 1.0 / maxEigVal;	
	const double tMultiplier = std::pow(maxT/minT, 1.0 / double(nScales - 1));

	std::vector<double> waveletScales(nScales);	
	for (int s = 0; s < nScales; ++s) waveletScales[s] = minT * std::pow(tMultiplier, s);

	matSGW.resize(totalAtomCount, vertCount);

	ZGeom::DenseMatrixd matEigVecs(eigCount, vertCount);	
	const double *pEigVals = &(mhb.getEigVals()[0]);
	double *pEigVec = matEigVecs.raw_ptr();
	for (int i = 0; i < eigCount; ++i) 
		std::copy_n(mhb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);


	//////////////////////////////////////////////////////////////////////////
	// compute SGW with AMP
	//
	std::vector<double> vDiag(eigCount);
	for (int s = 0; s < nScales; ++s) {
		for (int i = 0; i < eigCount; ++i) {
			vDiag[i] = generator1(waveletScales[s] * pEigVals[i]);
		}
		double *pResult = matSGW.raw_ptr() + vertCount * vertCount * s;
		ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
	}

	for (int i = 0; i < eigCount; ++i) {
		vDiag[i] = generator2(pEigVals[i]);
	}
	double *pResult = matSGW.raw_ptr() + vertCount * vertCount * nScales;
	ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
}

void DifferentialMeshProcessor::setMesh( CMesh* newMesh )
{
	mMesh = newMesh;
}

