#include "DifferentialMeshProcessor.h"
#include <fstream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <functional>
#include <ppl.h>
#include <mkl.h>
#include <ZUtil/SimpleConfigLoader.h>
#include <ZUtil/zassert.h>
#include <ZGeom/arithmetic.h>
#include <ZGeom/EigenSystem.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/MatVecArithmetic.h>
#include "global.h"

using namespace std;
using ZGeom::MatlabEngineWrapper;
using ZGeom::VectorPointwiseProduct;
using ZGeom::VectorDotProduct;

std::function<double(double)> transferScaling1 = [](double lambda) { return std::exp(-std::pow(lambda, 2)); };
std::function<double(double, double)> transferHeatKernel = [](double lambda, double t) { return std::exp(-lambda*t); };
std::function<double(double, double)> transferMHW = [](double lambda, double t) { return lambda * std::exp(-lambda*t); };

double transferScalingFunc1( double lambda )
{
	return std::exp(-std::pow(lambda, 2.0));
}

double transferFunc1( double lambda, double t )
{
	double coeff = (lambda*t) * (lambda*t);
	return coeff * std::exp(-coeff);
}

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
	mpEngineWrapper = NULL;
	mesh = NULL;
	mActiveFeature = FEATURE_ID;
	mRefVert = 0;
	m_size = 0;
	mActiveHandle = -1;	
}

DifferentialMeshProcessor::DifferentialMeshProcessor(CMesh* tm, CMesh* originalMesh)
{
	mpEngineWrapper = nullptr;
	mesh = nullptr;
	mActiveFeature = FEATURE_ID;
	mRefVert = 0;
	m_size = 0;
	mActiveHandle = -1;

	init_lite(tm, originalMesh);
}

DifferentialMeshProcessor::~DifferentialMeshProcessor()
{
	std::cout << "DifferentialMeshProcessor destroyed!" << std::endl;
}

void DifferentialMeshProcessor::init(CMesh* tm, MatlabEngineWrapper* e)
{
	mesh = tm;
	ori_mesh = tm;
	mpEngineWrapper = e;
	//matlabWrapper.setEngine(e);
	m_size = mesh->vertCount();
	mRefVert = g_configMgr.getConfigValueInt("INITIAL_REF_POINT");
	mRefPos = mesh->getVertex(mRefVert)->getPosition();
}

void DifferentialMeshProcessor::init_lite( CMesh* tm, CMesh* originalMesh )
{
	mesh = tm;
	ori_mesh = originalMesh;
	m_size = mesh->vertCount();
	mRefVert = 0;
	mRefPos = mesh->getVertex(0)->getPosition();
}

void DifferentialMeshProcessor::constructLaplacian( MeshLaplacian::LaplacianType laplacianType /*= CotFormula*/ )
{
	if (hasLaplacian(laplacianType)) return;
	
	MeshLaplacian& laplacian = mMeshLaplacians[laplacianType];
	switch(laplacianType)
	{
	case MeshLaplacian::Umbrella:
	case MeshLaplacian::Tutte:
	case MeshLaplacian::CotFormula:
	case MeshLaplacian::SymCot:
		(laplacian.*(laplacian.getConstructFunc(laplacianType)))(mesh);
		break;

	case MeshLaplacian::Anisotropic1:
		{
			double para1 = 2 * mesh->getAvgEdgeLength() * mesh->getAvgEdgeLength();
			double para2 = para1;
			mMeshLaplacians[laplacianType].constructFromMesh3(mesh, 1, para1, para2);
		}
		break;

	case MeshLaplacian::Anisotropic2:
		{
			double para1 = 2 * mesh->getAvgEdgeLength() * mesh->getAvgEdgeLength();
			double para2 = mesh->getAvgEdgeLength() / 2;
			mMeshLaplacians[laplacianType].constructFromMesh4(mesh, 1, para1, para2);
		}
		break;

	case MeshLaplacian::IsoApproximate:
		mMeshLaplacians[laplacianType].constructFromMesh5(mesh);                
		break;

	default: throw std::logic_error("Unrecognized Laplacian type");
	}       
}


void DifferentialMeshProcessor::decomposeLaplacian( int nEigFunc, MeshLaplacian::LaplacianType laplacianType /*= CotFormula*/ )
{
	ZUtil::logic_assert(hasLaplacian(laplacianType), "laplacian is not available for decomposition");
	if (!mpEngineWrapper->isOpened())
		throw std::logic_error("Matlab engine not opened for Laplacian decomposition!");
	
	Engine *ep = mpEngineWrapper->getEngine();
	mMeshLaplacians[laplacianType].decompose(nEigFunc, mpEngineWrapper, mMHBs[laplacianType]);
}

void DifferentialMeshProcessor::loadMHB( const std::string& path, MeshLaplacian::LaplacianType laplacianType /*= MeshLaplacian::CotFormula*/ )
{
	mMHBs[laplacianType].load(path);
}

void DifferentialMeshProcessor::saveMHB( const std::string& path, MeshLaplacian::LaplacianType laplacianType /*= MeshLaplacian::CotFormula*/ )
{
	mMHBs[laplacianType].save(path);	
}

void DifferentialMeshProcessor::computeCurvature( std::vector<double>& vCurvature, int curvatureType /*= 0*/ )
{
	vCurvature.resize(m_size);
	if (curvatureType == 0)	
		vCurvature = mesh->getMeanCurvature();
	else if (curvatureType == 1) 
		vCurvature = mesh->getGaussCurvature();
}

void DifferentialMeshProcessor::addNewHandle( int hIdx )
{
	auto iter = mHandles.find(hIdx);
	if (iter != mHandles.end()) mHandles.erase(iter);
	else mHandles.insert(std::make_pair(hIdx, mesh->getVertex(hIdx)->getPosition()));	 
}

void DifferentialMeshProcessor::computeSGW()
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);
	const int vertCount = mhb.eigVecSize();
	const int eigCount = mhb.eigVecCount();
	double waveletScales[] = {20, 40};
	const int scales = 1;

	mMatWavelet.resize(vertCount*(scales+1), vertCount);
	Concurrency::parallel_for(0, vertCount, [&](int i) {
		for (int s = 0; s < scales; ++s) {
			for (int j = 0; j <= i; ++j) {
				double elemVal(0);
				for (int k = 0; k < eigCount; ++k) {
					double lambda = mhb.getEigVal(k);
					const ZGeom::VecNd& phi = mhb.getEigVec(k);
					double lt = std::pow(lambda * waveletScales[s], 2);
					elemVal += lt * std::exp(-lt) * phi[i] * phi[j];
				}
				mMatWavelet(vertCount * s + i, j) = elemVal;
				mMatWavelet(vertCount * s + j, i) = elemVal;
			}
		}
#if 1
		for (int j = 0; j <= i; ++j) {
			double elemVal(0);
			for (int k = 0; k < eigCount; ++k) {	
				double lambda = mhb.getEigVal(k);
				const ZGeom::VecNd& phi = mhb.getEigVec(k);
				elemVal += std::exp(-lambda*lambda) * phi[i] * phi[j];
			}
			mMatWavelet(vertCount * scales + i, j) = elemVal;
			mMatWavelet(vertCount * scales + j, i) = elemVal;
		}
#endif
	});
}

void DifferentialMeshProcessor::calKernelSignature( double scale, KernelType kernelType, std::vector<double>& values ) const
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);

	values.resize(m_size);

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

	Concurrency::parallel_for(0, m_size, [&](int i)	{
		double sum = 0;
		for (int k = 0; k < mhb.eigVecCount(); ++k)
		{
			sum += (*pTF)(mhb.getEigVal(k), scale) * mhb.getEigVec(k)[i] * mhb.getEigVec(k)[i];
		}
		values[i] = sum;
	});
}

void DifferentialMeshProcessor::calNormalizedKernelSignature(double scale, KernelType kernelType, std::vector<double>& normalized_values) const
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);

	calKernelSignature(scale, kernelType, normalized_values);

	double normalize_factor = 0.;
	for (auto iter = normalized_values.begin(); iter != normalized_values.end(); ++iter)
	{
		normalize_factor += *iter;
	}

	double normalize_factor2 = 0.;
	for (auto iter = mhb.getEigVals().begin(); iter != mhb.getEigVals().end(); ++iter)
	{
		normalize_factor2 += std::exp(-*iter * scale);
	}
	cout << "factor1 = " << normalize_factor << ", " << "factor2 = " << normalize_factor2 << endl;

	std::transform( normalized_values.begin(), normalized_values.end(), normalized_values.begin(), [=](double v){return v / normalize_factor;} );
}

void DifferentialMeshProcessor::computeKernelSignatureFeatures( const std::vector<double>& timescales, KernelType kernelType )
{
	MeshFeatureList *mfl = new MeshFeatureList;
	int nScales = timescales.size();

	for (int s = 0; s < nScales; ++s)
	{
		vector<double> vSig;
		vector<int> vFeatures;
		calKernelSignature(timescales[s], kernelType, vSig);
		mesh->extractExtrema(vSig, 2, 1e-5, vFeatures);

		for (vector<int>::iterator iter = vFeatures.begin(); iter != vFeatures.end(); ++iter)
		{
			mfl->addFeature(new MeshFeature(*iter, s));			
		}
	}

	switch(kernelType)
	{
	case HEAT_KERNEL:
		removePropertyByID(FEATURE_HKS);
		mfl->setIDandName(FEATURE_HKS, "Feature_HKS");
		mActiveFeature = FEATURE_HKS;
		break;
	case MHW_KERNEL:
		removePropertyByID(FEATURE_MHWS);
		mfl->setIDandName(FEATURE_MHWS, "Feature_MHWS");
		mActiveFeature = FEATURE_MHWS;
		break;
	case SGW_KERNEL:
		removePropertyByID(FEATURE_SGWS);
		mfl->setIDandName(FEATURE_SGWS, "Feature_SGWS");
		mActiveFeature = FEATURE_SGWS;
		break;
	}

	addProperty(mfl);
}

double DifferentialMeshProcessor::calHK( int v1, int v2, double timescale ) const
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);
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
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k)	{
		double lambda = mhb.getEigVal(k);
		const ZGeom::VecNd& phi = mhb.getEigVec(k);
		sum += lambda * std::exp(-lambda * timescale) * phi[v1] * phi[v2];
	}
	return sum;
}

double DifferentialMeshProcessor::calSGW( int v1, int v2, double timescale ) const
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k)	{
		double lambda = mhb.getEigVal(k);
		const ZGeom::VecNd& phi = mhb.getEigVec(k);
		sum += std::pow(lambda * timescale, 2) * std::exp(-std::pow(lambda * timescale, 2)) * phi[v1] * phi[v2];
	}
	return sum;
}

double DifferentialMeshProcessor::calHeatTrace(double timescale) const
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k) {
		sum += std::exp(-mhb.getEigVal(k) * timescale);
	}
	return sum;
}

double DifferentialMeshProcessor::calBiharmonic(int v1, int v2) const
{
	const ManifoldHarmonics& mhb = getMHB(MeshLaplacian::CotFormula);
	double sum = 0;
	for (int k = 0; k < mhb.eigVecCount(); ++k) {
		sum += pow( (mhb.getEigVec(k)[v1] - mhb.getEigVec(k)[v2]) / mhb.getEigVal(k), 2 );
	}
	return std::sqrt(sum);
}

void DifferentialMeshProcessor::computeSimilarityMap1( int refPoint )
{
	const CVertex* pvi = mesh->getVertex(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);
	vFaces.resize(mesh->faceCount());
	for (int f = 0; f < mesh->faceCount(); ++f)
		vFaces[f] = f;

	vector<double> vSimilarities;
	vSimilarities.resize(mesh->vertCount(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mesh->vertCount(), 0.);

	double hPara1 = std::pow(mesh->getAvgEdgeLength() * 5, 2);
	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);
	const std::vector<Vector3D>& vVertNormals = mesh->getVertNormals();

	for (int fi = 0; fi < vFaces.size(); ++fi) {
		const CFace* pfi = mesh->getFace(vFaces[fi]);
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

	MeshFunction *mf = new MeshFunction(m_size);

	for (int i = 0; i < m_size; ++i) {
		mf->setValue(i, vSimilarities[i]);
	}

	removePropertyByID(SIGNATURE_SIMILARITY_MAP);
	mf->setIDandName(SIGNATURE_SIMILARITY_MAP, "Simialrity_Map_Signature");
	addProperty(mf);
}

void DifferentialMeshProcessor::computeSimilarityMap2( int refPoint )
{
	const CVertex* pvi = mesh->getVertex(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
	//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);
	vFaces.resize(mesh->faceCount());
	for (int f = 0; f < mesh->faceCount(); ++f) vFaces[f] = f;

	vector<double> vSimilarities;
	vSimilarities.resize(mesh->vertCount(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mesh->vertCount(), 0.);

	double hPara1 = std::pow(mesh->getAvgEdgeLength() * 5, 2);
	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);
	const std::vector<Vector3D>& vVertNormals = mesh->getVertNormals();
	const std::vector<Vector3D>& vFaceNormals = mesh->getFaceNormals();
	const int faceNum = mesh->faceCount();

	for (int fIndex = 0; fIndex < faceNum; ++fIndex) {
		const CFace* pfi = mesh->getFace(vFaces[fIndex]);
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

	MeshFunction *mf = new MeshFunction(m_size);

	for (int i = 0; i < m_size; ++i) {
		mf->setValue(i, vSimilarities[i]);
	}

	removePropertyByID(SIGNATURE_SIMILARITY_MAP);
	mf->setIDandName(SIGNATURE_SIMILARITY_MAP, "Simialrity_Map_Signature");
	addProperty(mf);
}

void DifferentialMeshProcessor::computeSimilarityMap3( int refPoint )
{
	const CVertex* pvi = mesh->getVertex(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
	vFaces.resize(mesh->faceCount());
	for (int f = 0; f < mesh->faceCount(); ++f)
		vFaces[f] = f;
//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);

	vector<double> vSimilarities;
	vSimilarities.resize(mesh->vertCount(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mesh->vertCount(), 0.);

	double hPara1 = std::pow(mesh->getAvgEdgeLength() * 5, 2);
//	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);
	double hPara2 = mesh->getAvgEdgeLength();
	const std::vector<Vector3D>& vVertNormals = mesh->getVertNormals();

	for (int fi = 0; fi < vFaces.size(); ++fi)
	{
		const CFace* pfi = mesh->getFace(vFaces[fi]);
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

	MeshFunction *mf = new MeshFunction(m_size);

	for (int i = 0; i < m_size; ++i)
	{
		mf->setValue(i, vSimilarities[i]);
	}

	removePropertyByID(SIGNATURE_SIMILARITY_MAP);
	mf->setIDandName(SIGNATURE_SIMILARITY_MAP, "Simialrity_Map_Signature");
	addProperty(mf);
}

const MeshFeatureList* DifferentialMeshProcessor::getActiveFeatures() const
{
	const MeshProperty* feat = retrievePropertyByID(mActiveFeature);
	if (feat == NULL) return NULL;
	else return dynamic_cast<const MeshFeatureList*>(feat);
}

void DifferentialMeshProcessor::calHeat( int vSrc, double tMultiplier, std::vector<double>& vHeat )
{
	const int vertCount = mesh->vertCount();
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
	const int vertCount = mesh->vertCount();
	const double t = std::pow(mesh->getAvgEdgeLength(), 2) * tMultiplier;

	const MeshLaplacian& laplacian = getMeshLaplacian(MeshLaplacian::CotFormula);
	const ZGeom::SparseMatrix<double>& matW = laplacian.getW();
	const ZGeom::SparseMatrix<double>& matLc = laplacian.getLS();	// negative

	ZGeom::addMatMat(matW, matLc, -t, mHeatDiffuseMat);	//A = W - t*Lc
	mHeatDiffuseSolver.initialize(mHeatDiffuseMat, true, true);
}
