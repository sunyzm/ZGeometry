#include "DifferentialMeshProcessor.h"
#include <fstream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <ppl.h>
#include <SimpleConfigLoader.h>

using namespace std;

extern SimpleConfigLoader g_configMgr;

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

double transferFunc3( double lambda, double t )
{
	return std::exp(-lambda * t);
}

double transferFunc4( double lambda, double t )
{
	return lambda * std::exp(-lambda * t);
}

TransferFunc HeatKernel_TF = &transferFunc3, MHW_TF = &transferFunc4;

DifferentialMeshProcessor::DifferentialMeshProcessor(void)
{
	m_ep = NULL;
	mesh = NULL;
	pvActiveFeatures = NULL;
	m_bLaplacianDecomposed = false;
	pRef = 0;
	m_size = 0;
	active_handle = -1;
	
	constrain_weight = 0.1;
	m_bSGWComputed = false;
}

DifferentialMeshProcessor::DifferentialMeshProcessor(CMesh* tm)
{
	DifferentialMeshProcessor();
	mesh = tm;

	init_lite(tm);
}

DifferentialMeshProcessor::~DifferentialMeshProcessor(void)
{
}

void DifferentialMeshProcessor::init(CMesh* tm, Engine* e)
{
	mesh = tm;
	m_ep = e;
	matlabWrapper.setEngine(e);
	m_size = mesh->getVerticesNum();
	pRef = g_configMgr.getConfigValueInt("INITIAL_REF_POINT");
	posRef = mesh->getVertex(pRef)->getPosition();
}

void DifferentialMeshProcessor::init_lite(CMesh* tm)
{
	mesh = tm;
	m_size = mesh->getVerticesNum();
	pRef = 0;
	posRef = mesh->getVertex(0)->getPosition();
}

void DifferentialMeshProcessor::decomposeLaplacian( int nEigFunc, LaplacianType laplacianType /*= CotFormula*/ )
{
	
//	mLaplacian.decompose(mhb, nEigFunc, m_ep);
//	meshKernel.decompose(mhb, nEigFunc, m_ep);

	vMeshLaplacian[laplacianType].decompose(vMHB[laplacianType], nEigFunc, m_ep);
	mhb = vMHB[laplacianType];

	this->m_bLaplacianDecomposed = true;
}

void DifferentialMeshProcessor::readMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	ifstream ifs(path, ios::binary);
	ifs.seekg (0, ios::end);
	int length = (int)ifs.tellg();
	ifs.seekg (0, ios::beg);
	char *buffer = new char[length];
	ifs.read(buffer, length);
	ifs.close();
	
	istringstream iss (buffer, istringstream::in);
	iss >> mhb.m_nEigFunc;
	iss >> mhb.m_size;
	mhb.m_func.resize(mhb.m_nEigFunc);
	for (int i = 0; i < mhb.m_nEigFunc; ++i)
	{
		mhb.m_func[i].m_vec.resize(mhb.m_size);
		iss >> mhb.m_func[i].m_val;
		for (int j = 0; j < mhb.m_size; ++j)
			iss >> mhb.m_func[i].m_vec[j];
	}
	delete []buffer;

	vMHB[laplacianType] = mhb;
	this->m_bLaplacianDecomposed = true;
}

void DifferentialMeshProcessor::writeMHB( const std::string& path, LaplacianType laplacianType /*= CotFormula*/ )
{
	ofstream ofs(path.c_str(), ios::trunc);
	ofs << mhb.m_nEigFunc << endl;
	ofs << mhb.m_size << endl;
	for (int i = 0; i < mhb.m_nEigFunc; ++i)
	{
		ofs << mhb.m_func[i].m_val;
		for (int j = 0; j < mhb.m_size; ++j)
			ofs << ' ' << mhb.m_func[i].m_vec[j];
		ofs << endl;
	}
	ofs.close();	
}

void DifferentialMeshProcessor::computeCurvature( std::vector<double>& vCurvature, int curvatureType /*= 0*/ )
{
	vCurvature.resize(m_size);
	if (curvatureType == 0)
	{
		for (int i = 0; i < m_size; ++i)
		{
//			mesh->calVertexCurvature(i);
			vCurvature[i] = mesh->getVertex(i)->getMeanCurvature();
		}
	}
	else if (curvatureType == 1)
	{
		for (int i = 0; i < m_size; ++i)
			vCurvature[i] = mesh->getVertex(i)->getGaussCurvature();
	}
	
}

void DifferentialMeshProcessor::addNewHandle( int hIdx )
{
	auto iter = mHandles.find(hIdx);
	if (iter != mHandles.end())
		mHandles.erase(iter);
	else
		mHandles[hIdx] = mesh->getVertex_const(hIdx)->getPosition();
	 
}

void DifferentialMeshProcessor::computeMexicanHatWavelet( std::vector<double>& vMHW, double scale, int wtype /*= 1*/ )
{
	vMHW.resize(m_size);

	if (wtype == 1)
	{
		for (int i = 0; i < m_size; ++i)
		{
			double sum = 0;
			for (int k = 0; k < mhb.m_nEigFunc; ++k)
			{
				double coef = mhb.m_func[k].m_val * scale;
				sum += mhb.m_func[k].m_val  * exp(-mhb.m_func[k].m_val * scale) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[pRef];
			}
			vMHW[i] = sum;
		}
	}
	else if (wtype == 2)
	{
		for (int i = 0; i < m_size; ++i)
		{
			double sum = 0;
			for (int k = 0; k < mhb.m_nEigFunc; ++k)
			{
				double coef = pow(mhb.m_func[k].m_val * scale, 2.0);
				sum += coef * exp(-coef) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[pRef];
			}
			vMHW[i] = sum;
		}
	}
}

void DifferentialMeshProcessor::computeExperimentalWavelet( std::vector<double>& vExp, double scale )
{
	vExp.resize(m_size);
	
	for (int i = 0; i < m_size; ++i)
	{
		double sum = 0;
		for (int k = 0; k < mhb.m_nEigFunc; ++k)
		{
			double coef = mhb.m_func[k].m_val * scale;
			sum += mhb.m_func[k].m_val * scale * exp(-mhb.m_func[k].m_val * scale) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[pRef];
		}
		vExp[i] = sum;
	}

	ofstream fout("output/frequency.txt", ios::trunc);
	for (int k = 0; k < mhb.m_nEigFunc; ++k)
	{
		double coef = pow(mhb.m_func[k].m_val * scale, 1.0);
		fout << coef * exp(-coef) /* * mhb.m_func[k].m_vec[pRef]*/ << " ";
	}
	fout << "\n";
	fout.close();
}

void DifferentialMeshProcessor::calGeometryDWT()
{
// output: 1. coordinates of all vertices
//         2. wavelets coefficients of the x-coordinates

	ofstream ofs("output/coord.dat");
	for (int i = 0; i < m_size; ++i)
		ofs << mesh->getVertex(i)->getPosition().x << ' ';
	ofs << endl;
	for (int i = 0; i < m_size; ++i)
		ofs << mesh->getVertex(i)->getPosition().y << ' ';
	ofs << endl;
	for (int i = 0; i < m_size; ++i)
		ofs << mesh->getVertex(i)->getPosition().z << ' ';
	ofs.close();

	int totalScales = 4;
	double lambdaMax = mhb.m_func.back().m_val;
	double lambdaMin = lambdaMax / pow(3.0, totalScales-1);
	vector<double> vScales;
	vScales.push_back(2.0/lambdaMin);
	for (int l = 1; l < totalScales; ++l)
	{
		vScales.push_back(vScales.back()/3);
	}

	ofs.open("output/wavelets.txt");

	for (int x = 0; x < m_size; ++x)
	{
		for (int y = 0; y < m_size; ++y)
		{
			double itemSum = 0;
			for (int k = 0; k < mhb.m_nEigFunc; ++k)
			{
//				itemSum += exp(-pow(mhb.m_func[k].m_val / (0.6*lambdaMin), 4)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
				itemSum += exp(-1.0) * exp(-pow(mhb.m_func[k].m_val / (0.6*lambdaMin), 4)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
			}
			ofs << itemSum << ' ';
		}
		ofs << endl;

		for (int l = 0; l < totalScales; ++l)
		{
			for (int y = 0; y < m_size; ++y)
			{
				double itemSum = 0;
				for (int k = 0; k < mhb.m_nEigFunc; ++k)
				{
//					itemSum += TransferFunc1(mhb.m_func[k].m_val * vScales[l]) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
					itemSum += pow(mhb.m_func[k].m_val * vScales[l], 2) * exp(-pow(mhb.m_func[k].m_val * vScales[l],2)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
				}
				ofs << itemSum << ' ';
			}
			ofs << endl;
		}
	}
}

void DifferentialMeshProcessor::reconstructExperimental1( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint /*= false*/ ) const
{
	const MeshLaplacian& mLaplacian = vMeshLaplacian[CotFormula];

	vx.resize(m_size);
	vy.resize(m_size);
	vz.resize(m_size);

	vector<double> vxcoord0, vycoord0, vzcoord0;
	mesh->getCoordinateFunction(0, vxcoord0);
	mesh->getCoordinateFunction(1, vycoord0);
	mesh->getCoordinateFunction(2, vzcoord0);
	
	vector<double> vWeight = mLaplacian.getVerticesWeight();
	vector<double> xWeightedCoord, yWeightedCoord, zWeightedCoord;
	VectorPointwiseProduct(vxcoord0, vWeight, xWeightedCoord);
	VectorPointwiseProduct(vycoord0, vWeight, yWeightedCoord);
	VectorPointwiseProduct(vzcoord0, vWeight, zWeightedCoord);

 	int scales = 1;
 	double t_scales[4] = {80, 40, 20, 10};
 	vector<vector<double> > SGW;
 	
//	scaling function
 	for (int x = 0; x < m_size; ++x)
 	{
 		SGW.push_back(vector<double>());
		SGW.back().resize(m_size);
 
 		for (int y = 0; y < m_size; ++y)
 		{
 			double itemSum = 0;
 			for (int k = 0; k < mhb.m_nEigFunc; ++k)
 			{
 //				itemSum += exp(-pow(mhb.m_func[k].m_val / (0.6*lambdaMin), 4)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
 //				itemSum += exp(-1.0) * exp(-pow(mhb.m_func[k].m_val / (0.6*lambdaMin), 4)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
 				itemSum += exp(-1.0) * exp(-mhb.m_func[k].m_val) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
 			}
 			SGW.back().at(y) = itemSum;
 		}
 	}
 
//	wavelet functions
 	for (int s = 0; s < scales; ++s)
 	{
 		for (int x = 0; x < m_size; ++x)
 		{
			SGW.push_back(vector<double>());
			SGW.back().resize(m_size);
 
 			for (int y = 0; y < m_size; ++y)
 			{
 				double itemSum = 0;
 				for (int k = 0; k < mhb.m_nEigFunc; ++k)
 				{
					double coef = mhb.m_func[k].m_val * t_scales[s];
 					itemSum += coef * exp(-coef) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
 				}
 				SGW.back().at(y) = itemSum;
 			}
 		}
 	}
 
	int sizeCoeff = (scales + 1) * m_size;
 	vector<double> vxCoeff, vyCoeff, vzCoeff;
	vxCoeff.resize(sizeCoeff);
	vyCoeff.resize(sizeCoeff);
	vzCoeff.resize(sizeCoeff);

	assert(SGW.size() == sizeCoeff);

	for (int i = 0; i < sizeCoeff; ++i)
	{
		double itemSumX = 0, itemSumY = 0, itemSumZ = 0;
		for (int j = 0; j < m_size; ++j)
		{
// 				itemSumX += SGW[i][j] * xWeightedCoord[j];
// 				itemSumY += SGW[i][j] * yWeightedCoord[j];
// 				itemSumZ += SGW[i][j] * zWeightedCoord[j];

				itemSumX += SGW[i][j] * vxcoord0[j];
				itemSumY += SGW[i][j] * vycoord0[j];
				itemSumZ += SGW[i][j] * vzcoord0[j];
		}
		vxCoeff[i] = itemSumX;
		vyCoeff[i] = itemSumY;
		vzCoeff[i] = itemSumZ;
	}

	double weightI = 0.1;
	if (withConstraint)
	{
		SGW.push_back(vector<double>());
		SGW.back().resize(m_size, 0.0);
		SGW.back().at(pRef) = weightI;
		
		vxCoeff.push_back(posRef.x * weightI);
		vyCoeff.push_back(posRef.y * weightI);
		vzCoeff.push_back(posRef.z * weightI);
	}
	
/*
	ofstream of1("output/sgw.dat"), of2("output/coeff.dat");
	for (int i = 0; i < sizeCoeff; ++i)
	{
		for (int j = 0; j < m_size; ++j)
		{
			of1 << SGW[i].at[j] << ' ';
		}
		of1 << endl;
	}

	for (int i = 0; i < sizeCoeff; ++i)
	{
		of2 << vxCoeff[i] << endl;
	}

	VectorPointwiseDivide(xWeightedCoord, vWeight, vx);
	VectorPointwiseDivide(yWeightedCoord, vWeight, vy);
	VectorPointwiseDivide(zWeightedCoord, vWeight, vz);
*/

// 	matlab_cgls(m_ep, SGW, vxCoeff, vx);
// 	matlab_cgls(m_ep, SGW, vyCoeff, vy);
// 	matlab_cgls(m_ep, SGW, vzCoeff, vz);

	matlab_scgls(m_ep, SGW, vxCoeff, vx);
	matlab_scgls(m_ep, SGW, vyCoeff, vy);
	matlab_scgls(m_ep, SGW, vzCoeff, vz);	
 
//  VectorPointwiseDivide(vx, vWeight, vx);
//  VectorPointwiseDivide(vy, vWeight, vy);
//  VectorPointwiseDivide(vz, vWeight, vz);
}

void DifferentialMeshProcessor::reconstructByMHB( int aN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz ) const
{
	const MeshLaplacian& mLaplacian = vMeshLaplacian[CotFormula];

	vx.resize(m_size);
	vy.resize(m_size);
	vz.resize(m_size);

	vector<double> vxcoord0, vycoord0, vzcoord0;
	mesh->getCoordinateFunction(0, vxcoord0);
	mesh->getCoordinateFunction(1, vycoord0);
	mesh->getCoordinateFunction(2, vzcoord0);

	vector<double> vWeight = mLaplacian.getVerticesWeight();
	vector<double> xWeightedCoord, yWeightedCoord, zWeightedCoord;
	VectorPointwiseProduct(vWeight, vxcoord0, xWeightedCoord);
	VectorPointwiseProduct(vWeight, vycoord0, yWeightedCoord);
	VectorPointwiseProduct(vWeight, vzcoord0, zWeightedCoord);

	vector<double> xCoeff, yCoeff, zCoeff;
	int approxN = min(aN, mhb.m_nEigFunc);
	for (int k = 0; k < approxN; ++k)
	{
		xCoeff.push_back(VectorDotProduct(xWeightedCoord, mhb.m_func[k].m_vec));
		yCoeff.push_back(VectorDotProduct(yWeightedCoord, mhb.m_func[k].m_vec));
		zCoeff.push_back(VectorDotProduct(zWeightedCoord, mhb.m_func[k].m_vec));
	}

	for (int i = 0; i < m_size; ++i)
	{
		double sumX(0), sumY(0), sumZ(0);
		for (int k = 0; k < approxN; ++k)
		{
			sumX += xCoeff[k] * mhb.m_func[k].m_vec[i];
			sumY += yCoeff[k] * mhb.m_func[k].m_vec[i];
			sumZ += zCoeff[k] * mhb.m_func[k].m_vec[i];
		}
		vx[i] = sumX; 
		vy[i] = sumY;
		vz[i] = sumZ;
	}
	
}

void DifferentialMeshProcessor::reconstructByDifferential( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint /* =false */ ) const
{
	const MeshLaplacian& mLaplacian = vMeshLaplacian[CotFormula];

	vx.resize(m_size);
	vy.resize(m_size);
	vz.resize(m_size);

	vector<double> vxcoord0, vycoord0, vzcoord0;
	mesh->getCoordinateFunction(0, vxcoord0);
	mesh->getCoordinateFunction(1, vycoord0);
	mesh->getCoordinateFunction(2, vzcoord0);
	
 	vector<vector<double> > sysM;
 	
	for (int i = 0; i < m_size; ++i)
	{
		sysM.push_back(vector<double>());
		sysM.back().resize(m_size, 0.0);
	}

	vector<int> vII, vJJ;
	vector<double> vSS;
	mLaplacian.getSparseMatrix(vII, vJJ, vSS);

	int sparseLapSize = vII.size();
	for (int n = 0; n < sparseLapSize; ++n)
	{
		sysM.at(vII[n]-1).at(vJJ[n]-1) += vSS[n];
	}

	vector<double> vxCoeff, vyCoeff, vzCoeff;
	vxCoeff.resize(m_size);
	vyCoeff.resize(m_size);
	vzCoeff.resize(m_size);

	for (int i = 0; i < m_size; ++i)
	{
		double itemSumX = 0, itemSumY = 0, itemSumZ = 0;
		for (int j = 0; j < m_size; ++j)
		{
				itemSumX += sysM[i][j] * vxcoord0[j];
				itemSumY += sysM[i][j] * vycoord0[j];
				itemSumZ += sysM[i][j] * vzcoord0[j];
		}
		vxCoeff[i] = itemSumX;
		vyCoeff[i] = itemSumY;
		vzCoeff[i] = itemSumZ;
	}

	if (withConstraint)
	{
		double weightI = 1.0;
		sysM.push_back(vector<double>());
		sysM.back().resize(m_size, 0.0);
		sysM.back().at(pRef) = weightI;
		
		vxCoeff.push_back(posRef.x * weightI);
		vyCoeff.push_back(posRef.y * weightI);
		vzCoeff.push_back(posRef.z * weightI);
	}
	
	matlab_scgls(m_ep, sysM, vxCoeff, vx);
	matlab_scgls(m_ep, sysM, vyCoeff, vy);
	matlab_scgls(m_ep, sysM, vzCoeff, vz);	
}

void DifferentialMeshProcessor::reconstructBySGW( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint /*= false*/ )
{
	vx.resize(m_size);
	vy.resize(m_size);
	vz.resize(m_size);

	vector<double> vxcoord0, vycoord0, vzcoord0;
	mesh->getCoordinateFunction(0, vxcoord0);
	mesh->getCoordinateFunction(1, vycoord0);
	mesh->getCoordinateFunction(2, vzcoord0);
	
	if (!m_bSGWComputed)
		throw logic_error("SGW not computed!");
	 	
 	vector<vector<double> > SGW = m_vSGW;
	int scales = m_vTimescales.size();

	int sizeCoeff = SGW.size();

	vector<double> vxCoeff, vyCoeff, vzCoeff;
	vxCoeff.resize(sizeCoeff);
	vyCoeff.resize(sizeCoeff);
	vzCoeff.resize(sizeCoeff);
	
	for (int i = 0; i < sizeCoeff; ++i)
	{
		double itemSumX = 0, itemSumY = 0, itemSumZ = 0;
		for (int j = 0; j < m_size; ++j)
		{
				itemSumX += SGW[i][j] * vxcoord0[j];
				itemSumY += SGW[i][j] * vycoord0[j];
				itemSumZ += SGW[i][j] * vzcoord0[j];
		}
		vxCoeff[i] = itemSumX;
		vyCoeff[i] = itemSumY;
		vzCoeff[i] = itemSumZ;
	}

	double weightI = constrain_weight;
	if (withConstraint)
	{
		for (auto iter = mHandles.begin(); iter != mHandles.end(); ++iter)
		{
			SGW.push_back(vector<double>());
			SGW.back().resize(m_size, 0.0);
			SGW.back().at(iter->first) = weightI;
			vxCoeff.push_back(iter->second.x * weightI);
			vyCoeff.push_back(iter->second.y * weightI);
			vzCoeff.push_back(iter->second.z * weightI);
		}
	}

 	matlabWrapper.DenseConjugateGradient(SGW, vxCoeff, vx);
	matlabWrapper.DenseConjugateGradient(SGW, vyCoeff, vy);
	matlabWrapper.DenseConjugateGradient(SGW, vzCoeff, vz);

//	matlab_scgls(m_ep, SGW, vxCoeff, vx);
//	matlab_scgls(m_ep, SGW, vyCoeff, vy);
//	matlab_scgls(m_ep, SGW, vzCoeff, vz);	
}

void DifferentialMeshProcessor::filterBySGW( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz )
{
	vx.resize(m_size);
	vy.resize(m_size);
	vz.resize(m_size);

	vector<double> vxcoord0, vycoord0, vzcoord0;
	mesh->getCoordinateFunction(0, vxcoord0);
	mesh->getCoordinateFunction(1, vycoord0);
	mesh->getCoordinateFunction(2, vzcoord0);

	if (m_vSGW.empty())
	{
		vector<double> t_scales;
		t_scales.push_back(80);
//		t_scales.push_back(40);
//		t_scales.push_back(20);
//		t_scales.push_back(10);
		computeSGW(t_scales);
	}

	vector<vector<double> > SGW = m_vSGW;
	int nScales = m_vTimescales.size();
	int sizeCoeff = SGW.size();

	vector<double> vxCoeff, vyCoeff, vzCoeff;
	vxCoeff.resize(sizeCoeff);
	vyCoeff.resize(sizeCoeff);
	vzCoeff.resize(sizeCoeff);

	for (int i = 0; i < sizeCoeff; ++i)
	{
		double itemSumX = 0, itemSumY = 0, itemSumZ = 0;
		for (int j = 0; j < m_size; ++j)
		{
			itemSumX += SGW[i][j] * vxcoord0[j];
			itemSumY += SGW[i][j] * vycoord0[j];
			itemSumZ += SGW[i][j] * vzcoord0[j];
		}
		vxCoeff[i] = itemSumX;
		vyCoeff[i] = itemSumY;
		vzCoeff[i] = itemSumZ;

//		if ( i >= m_size && i < m_size * 2)
 		if ( i >= sizeCoeff - m_size )
 		{
 			vxCoeff[i] *= 0.5;
 			vyCoeff[i] *= 0.5;
 			vzCoeff[i] *= 0.5;
 		}
 	}

	matlab_scgls(m_ep, SGW, vxCoeff, vx);
	matlab_scgls(m_ep, SGW, vyCoeff, vy);
	matlab_scgls(m_ep, SGW, vzCoeff, vz);
}

void DifferentialMeshProcessor::deform( const std::vector<int>& vHandleIdx, const std::vector<Vector3D>& vHandlePos, const std::vector<int>& vFreeIdx, std::vector<Vector3D>& vDeformedPos, DeformType dfType )
{
	const MeshLaplacian& mLaplacian = vMeshLaplacian[CotFormula];

	if (vHandleIdx.size() != vHandlePos.size())
		throw logic_error("Error: DifferentialMeshProcessor::deform; parameters size incompatible!");
	
	vDeformedPos.clear();
	
	if (dfType == Simple)
	{
		int hIdx = vHandleIdx[0];		// only use the first handle to deform
		Vector3D handleTrans = vHandlePos[0] - mesh->getVertex_const(hIdx)->getPosition();
		
		int nFreeVertices = vFreeIdx.size();
		vector<double> vDist2Handle;
		double distMax = 0;
		for (int i = 0; i < nFreeVertices; ++i)
		{
			double dist = mesh->getGeodesic(hIdx, vFreeIdx[i]);
			vDist2Handle.push_back(dist);
			if (dist > distMax) distMax = dist;
		}
				
		for (int i = 0; i < nFreeVertices; ++i)
		{
			Vector3D newPos = mesh->getVertex_const(vFreeIdx[i])->getPosition() + handleTrans * (1.0 - vDist2Handle[i]/distMax);
			vDeformedPos.push_back(newPos);
		}
		
	}
	else if (dfType == Laplace)
	{
		vector<int> vI, vJ;
		vector<double> vS, vBX, vBY, vBZ;
		vector<double> vRX, vRY, vRZ;
		int nFree = vFreeIdx.size();
		int n = this->m_size, //m = n * 2 - nFree;
			m = n + vHandleIdx.size();
				
		int nz_laplacian = vI.size();
		vI = mLaplacian.vII;
		vJ = mLaplacian.vJJ;
		vS = mLaplacian.vSS;
		
		vector<double> xCoord, yCoord, zCoord;
		mesh->getCoordinateFunction(0, xCoord);
		mesh->getCoordinateFunction(1, yCoord);
		mesh->getCoordinateFunction(2, zCoord);
		matlabWrapper.SparseMatVecMul(n, n, vI, vJ, vS, xCoord, vBX);
		matlabWrapper.SparseMatVecMul(n, n, vI, vJ, vS, yCoord, vBY);
		matlabWrapper.SparseMatVecMul(n, n, vI, vJ, vS, zCoord, vBZ);

		int posII = n;
		for (int h = 0; h < vHandleIdx.size(); ++h)
		{
			vI.push_back(posII);
			vJ.push_back(vHandleIdx[h]);
			vS.push_back(1);
			vBX.push_back(vHandlePos[h].x);
			vBY.push_back(vHandlePos[h].y);
			vBZ.push_back(vHandlePos[h].z);
			posII += 1;
		}

		matlabWrapper.SparseBiConjugateGradient(m, n, vI, vJ, vS, vBX, vRX);		//A*R=B
		matlabWrapper.SparseBiConjugateGradient(m, n, vI, vJ, vS, vBY, vRY);
		matlabWrapper.SparseBiConjugateGradient(m, n, vI, vJ, vS, vBZ, vRZ);

		for (int i = 0; i < vFreeIdx.size(); ++i)
		{
			Vector3D newPos = Vector3D(vRX[vFreeIdx[i]], vRY[vFreeIdx[i]], vRZ[vFreeIdx[i]]);
			vDeformedPos.push_back(newPos);
		}

	}
	else if (dfType == SGW)
	{
		if (!m_bSGWComputed)
			throw logic_error("Error: SGW not computed!");

		vector<double> vRX, vRY, vRZ;

		vRX.resize(m_size);
		vRY.resize(m_size);
		vRZ.resize(m_size);

		vector<double> vxcoord0, vycoord0, vzcoord0;
		mesh->getCoordinateFunction(0, vxcoord0);
		mesh->getCoordinateFunction(1, vycoord0);
		mesh->getCoordinateFunction(2, vzcoord0);

		vector<vector<double> > SGW = m_vSGW;
		int scales = m_vTimescales.size();

		int sizeCoeff = SGW.size();

		vector<double> vxCoeff, vyCoeff, vzCoeff;
		vxCoeff.resize(sizeCoeff);
		vyCoeff.resize(sizeCoeff);
		vzCoeff.resize(sizeCoeff);

		for (int i = 0; i < sizeCoeff; ++i)
		{
			double itemSumX = 0, itemSumY = 0, itemSumZ = 0;
			for (int j = 0; j < m_size; ++j)
			{
				itemSumX += SGW[i][j] * vxcoord0[j];
				itemSumY += SGW[i][j] * vycoord0[j];
				itemSumZ += SGW[i][j] * vzcoord0[j];
			}
			vxCoeff[i] = itemSumX;
			vyCoeff[i] = itemSumY;
			vzCoeff[i] = itemSumZ;
		}

		double weightI = constrain_weight;
		int nHandleSize = vHandleIdx.size();
		for (int i = 0; i < nHandleSize; ++i)
		{
			SGW.push_back(vector<double>());
			SGW.back().resize(m_size, 0.0);
			SGW.back().at(vHandleIdx[i]) = weightI;

			vxCoeff.push_back(vHandlePos[i].x * weightI);
			vyCoeff.push_back(vHandlePos[i].y * weightI);
			vzCoeff.push_back(vHandlePos[i].z * weightI);
		}

		matlabWrapper.DenseConjugateGradient(SGW, vxCoeff, vRX);
		matlabWrapper.DenseConjugateGradient(SGW, vyCoeff, vRY);
		matlabWrapper.DenseConjugateGradient(SGW, vzCoeff, vRZ);

		for (int i = 0; i < vFreeIdx.size(); ++i)
		{
			Vector3D newPos = Vector3D(vRX[vFreeIdx[i]], vRY[vFreeIdx[i]], vRZ[vFreeIdx[i]]);
			vDeformedPos.push_back(newPos);
		}
	}
}

void DifferentialMeshProcessor::computeSGW( const std::vector<double>& timescales, double (*transferWavelet)(double, double) /*= &transferFunc1*/, bool withScaling /*= false*/, double (*transferScaling)(double) /*= &transferScalingFunc1*/ )
{
	m_vTimescales = timescales;
	int nScales = m_vTimescales.size();

	m_vSGW.clear();

	//	scaling function
	if (withScaling)	
	{
		for (int x = 0; x < m_size; ++x)
		{
			m_vSGW.push_back(vector<double>());
			m_vSGW.back().resize(m_size);
		}
		for (int x = 0; x < m_size; ++x)
		{
			for (int y = 0; y <= x; ++y)
			{
				double itemSum = 0;
				for (int k = 0; k < mhb.m_nEigFunc; ++k)
				{
					double transfer = (*transferScaling)(mhb.m_func[k].m_val);
					itemSum +=  transfer * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
				}
				m_vSGW[x][y] = m_vSGW[y][x] = itemSum;
			}
		}
	}	

	//	wavelet functions: nm * m
	for (int s = 0; s < nScales; ++s)
	{
		int offset = withScaling ? m_size * (s+1) : m_size * s;
		
		for (int i = 0; i < m_size; ++i)
		{
			m_vSGW.push_back(vector<double>());
			m_vSGW.back().resize(m_size);
		}

		for (int x = 0; x < m_size; ++x)
		{
			for (int y = 0; y <= x; ++y)
			{
				double itemSum = 0;
				for (int k = 0; k < mhb.m_nEigFunc; ++k)
				{
					double transfer = (*transferWavelet)(mhb.m_func[k].m_val, m_vTimescales[s]);
					itemSum += transfer * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
				}
				m_vSGW[offset + x][y] = m_vSGW[offset + y][x] = itemSum;
			}
		}
	}

	m_bSGWComputed = true;
}

void DifferentialMeshProcessor::calKernelSignature( double timescale, KernelType kernelType, std::vector<double>& values ) const
{
	if (!m_bLaplacianDecomposed) return;
	values.resize(m_size);

	TransferFunc pTF = &transferFunc3;	// default as heat kernel

	switch(kernelType)
	{
	case HEAT_KERNEL:
		pTF = &transferFunc3;
		break;
	case MHW_KERNEL:
		pTF = &transferFunc4;
		break;
	case SGW_KERNEL:
		pTF = &transferFunc2;
		break;
	}

	Concurrency::parallel_for(0, m_size, [&](int i)
//	for (int i = 0; i < m_size; ++i)
	{
		double sum = 0;
		for (int k = 0; k < mhb.m_nEigFunc; ++k)
		{
			sum += (*pTF)(mhb.m_func[k].m_val, timescale) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[i];
		}
		values[i] = sum;
	}
	);
}

void DifferentialMeshProcessor::computeKernelSignature( double timescale, KernelType kernelType )
{
	if (!m_bLaplacianDecomposed) return;

	MeshFunction *mf = new MeshFunction(m_size);

	switch(kernelType)
	{
	case HEAT_KERNEL:
		mf->setIDandName(SIGNATURE_HKS, "HKS");		
		break;

	case MHW_KERNEL:
		mf->setIDandName(SIGNATURE_MHWS, "Mexican_Hat_Signature");
		break;
	
	case SGW_KERNEL:
		mf->setIDandName(SIGNATURE_SGWS, "Spectral_Graph_Wavelet_Signature");
		break;
	}
	
	calKernelSignature(timescale, kernelType, mf->getMeshFunction());
	replaceProperty(mf);	
}

void DifferentialMeshProcessor::computeKernelDistanceSignature( double timescale, KernelType kernelType, int refPoint )
{
	if (!m_bLaplacianDecomposed) return;
	if (refPoint < 0 || refPoint >= m_size)
		throw runtime_error("Error computeKernelDistanceSignature: invalid reference point");

	MeshFunction *mf = new MeshFunction(m_size);
	TransferFunc pTF = &transferFunc3;

	switch(kernelType)
	{
	case HEAT_KERNEL:
		pTF = &transferFunc3;
		mf->setIDandName(SIGNATURE_HK, "HK");		
		break;
	case MHW_KERNEL:
		pTF = &transferFunc4;
		mf->setIDandName(SIGNATURE_MHW, "MHW");		
		break;
	case SGW_KERNEL:
		pTF = &transferFunc4;
		mf->setIDandName(SIGNATURE_SGW, "SGW");		
		break;
	}
	
	for (int i = 0; i < m_size; ++i)
	{
		double sum = 0;
		for (int k = 0; k < mhb.m_nEigFunc; ++k)
		{
			sum += (*pTF)(mhb.m_func[k].m_val, timescale) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[refPoint];
		}
		mf->setValue(i, sum);
	}

	replaceProperty(mf);
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
		break;
	case MHW_KERNEL:
		removePropertyByID(FEATURE_MHWS);
		mfl->setIDandName(FEATURE_MHWS, "Feature_MHWS");
		break;
	case SGW_KERNEL:
		removePropertyByID(FEATURE_SGWS);
		mfl->setIDandName(FEATURE_SGWS, "Feature_SGWS");
		break;
	}

	addProperty(mfl);

	pvActiveFeatures = mfl->getFeatureVector();
}

void DifferentialMeshProcessor::setActiveFeaturesByID( int feature_id )
{
	if (retrievePropertyByID(feature_id))
	{
		MeshFeatureList* mfl = dynamic_cast<MeshFeatureList*>(retrievePropertyByID(feature_id));
		pvActiveFeatures = mfl->getFeatureVector();
	}
}

double DifferentialMeshProcessor::calHK( int v1, int v2, double timescale ) const
{
	double sum = 0;
	for (int k = 0; k < mhb.m_nEigFunc; ++k)
	{
		sum += std::exp(-mhb.m_func[k].m_val * timescale) * mhb.m_func[k].m_vec[v1] * mhb.m_func[k].m_vec[v2];
	}
	return sum;
}

double DifferentialMeshProcessor::calBiharmonic(int v1, int v2) const
{
	double sum = 0;
	for (int k = 0; k < mhb.m_nEigFunc; ++k)
	{
		sum += pow( (mhb.m_func[k].m_vec[v1] - mhb.m_func[k].m_vec[v2]) / mhb.m_func[k].m_val, 2 );
	}
	return sum;
}

void DifferentialMeshProcessor::computeBiharmonicDistanceSignature( int refPoint )
{
	if (!m_bLaplacianDecomposed) return;
	if (refPoint < 0 || refPoint >= m_size)
		throw runtime_error("Error computeBiharmonicDistanceSignature: invalid reference point");

	MeshFunction *mf = new MeshFunction(m_size);
	
	for (int i = 0; i < m_size; ++i)
	{
		double sum = 0;
		for (int k = 0; k < mhb.m_nEigFunc; ++k)
		{
			sum += pow((mhb.m_func[k].m_vec[i] - mhb.m_func[k].m_vec[refPoint]) / mhb.m_func[k].m_val, 2);
		}
		mf->setValue(i, sum);
	}

	removePropertyByID(SIGNATURE_BIHARMONIC_DISTANCE);
	mf->setIDandName(SIGNATURE_BIHARMONIC_DISTANCE, "Biharmonic_Distance_signature");
	addProperty(mf);
}

void DifferentialMeshProcessor::computeSimilarityMap1( int refPoint )
{
	const CVertex* pvi = mesh->getVertex_const(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);
	vFaces.resize(mesh->getFaceNum());
	for (int f = 0; f < mesh->getFaceNum(); ++f)
		vFaces[f] = f;

	vector<double> vSimilarities;
	vSimilarities.resize(mesh->getVerticesNum(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mesh->getVerticesNum(), 0.);

	double hPara1 = std::pow(mesh->getAvgEdgeLength() * 5, 2);
	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);

	for (int fi = 0; fi < vFaces.size(); ++fi)
	{
		const CFace* pfi = mesh->getFace_const(vFaces[fi]);
		double face_area = pfi->getArea();
		for (int k = 0; k < 3; ++k)
		{
			int vki = pfi->getVertexIndex(k);
//			if (vki == refPoint) continue;
			const CVertex* pvk = pfi->getVertex_const(k);

			double w1 = 1., w2 = 1.;
//			double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
			w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
//			w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
			w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
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

void DifferentialMeshProcessor::computeSimilarityMap2( int refPoint )
{
	const CVertex* pvi = mesh->getVertex_const(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
	//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);
	vFaces.resize(mesh->getFaceNum());
	for (int f = 0; f < mesh->getFaceNum(); ++f)
		vFaces[f] = f;

	vector<double> vSimilarities;
	vSimilarities.resize(mesh->getVerticesNum(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mesh->getVerticesNum(), 0.);

	double hPara1 = std::pow(mesh->getAvgEdgeLength() * 5, 2);
	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);

	for (int fi = 0; fi < vFaces.size(); ++fi)
	{
		const CFace* pfi = mesh->getFace_const(vFaces[fi]);
		double face_area = pfi->getArea();
		for (int k = 0; k < 3; ++k)
		{
			int vki = pfi->getVertexIndex(k);
			//			if (vki == refPoint) continue;
			const CVertex* pvk = pfi->getVertex_const(k);

			double w1 = 1., w2 = 1.;
//			double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
			w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
//			w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
//			w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
			w2 = std::exp((dotProduct3D(pvi->getNormal(), pfi->getNormal()) - 1) / 1.0);

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

void DifferentialMeshProcessor::computeSimilarityMap3( int refPoint )
{
	const CVertex* pvi = mesh->getVertex_const(refPoint); 
	int ringT = 10;
	vector<int> vFaces;
	vFaces.resize(mesh->getFaceNum());
	for (int f = 0; f < mesh->getFaceNum(); ++f)
		vFaces[f] = f;
//	vFaces = mesh->getVertexAdjacentFacesIndex(refPoint, ringT);

	vector<double> vSimilarities;
	vSimilarities.resize(mesh->getVerticesNum(), 1.0);

	vector<double> vAreas;
	vAreas.resize(mesh->getVerticesNum(), 0.);

	double hPara1 = std::pow(mesh->getAvgEdgeLength() * 5, 2);
//	double hPara2 = std::pow(mesh->getAvgEdgeLength(), 2);
	double hPara2 = mesh->getAvgEdgeLength();

	for (int fi = 0; fi < vFaces.size(); ++fi)
	{
		const CFace* pfi = mesh->getFace_const(vFaces[fi]);
		double face_area = pfi->getArea();
		for (int k = 0; k < 3; ++k)
		{
			int vki = pfi->getVertexIndex(k);
			//			if (vki == refPoint) continue;
			const CVertex* pvk = pfi->getVertex_const(k);

			double w1 = 1., w2 = 1.;
			//			double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
			w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
			//			w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
			w2 = std::exp(-dotProduct3D(pvi->getNormal(), pvk->getPosition() - pvi->getPosition()) / hPara2);
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

