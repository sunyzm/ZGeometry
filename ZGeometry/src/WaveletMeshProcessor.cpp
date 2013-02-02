#include "WaveletMeshProcessor.h"
#include <fstream>

using namespace std;

void matlab_cgls(Engine* ep, const vector<vector<double> >& SGW, const vector<double>& b, vector<double>& sol)
{
	//solve A*sol = B

 	assert(SGW.size() == b.size());

	int sizeA = SGW.size(), sizeF = SGW[0].size();
	sol.resize(sizeF);

	mxArray *AA, *BB;
	AA = mxCreateDoubleMatrix(sizeA, sizeF, mxREAL);
	BB = mxCreateDoubleMatrix(sizeA, 1, mxREAL);
	double *pAA = mxGetPr(AA), *pBB = mxGetPr(BB);

	for (int i = 0; i < sizeA; ++i)
	{
		for (int j = 0; j < sizeF; ++j)
		{
			pAA[i + j*sizeA] = SGW[i][j];
		}
	}

	for (int j = 0; j < sizeA; ++j)
	{
		pBB[j] = b[j];
	}

	engPutVariable(ep, "A", AA);
	engPutVariable(ep, "b", BB);
	
	engEvalString(ep, "evals = cgls(A, b);");

	mxArray *SS = engGetVariable(ep, "evals");
	double *pSS = mxGetPr(SS);
	for (int j = 0; j < sizeF; ++j)
	{
		sol[j] = pSS[j];
	}

	mxDestroyArray(AA);
	mxDestroyArray(BB);
}

void matlab_scgls(Engine* ep, const vector<vector<double> >& SGW, const vector<double>& b, vector<double>& sol)
{
	//solve A*sol = B

	assert(SGW.size() == b.size());

	int sizeA = SGW.size(), sizeF = SGW[0].size();
	sol.resize(sizeF);

	const double sEpsilon = 1e-5;
	vector<double> vII, vJJ, vSS;
	for (int i = 0; i < sizeA; ++i)
	{
		for (int j = 0; j < sizeF; ++j)
		{
			if ( fabs(SGW[i][j]) > sEpsilon )
			{
				vII.push_back(i+1);
				vJJ.push_back(j+1);
				vSS.push_back(SGW[i][j]);
			}
		}
	}

	int sizeS = vSS.size();
	mxArray *II, *JJ, *SS, *BB, *DimM, *DimN;
	II = mxCreateDoubleMatrix(sizeS, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(sizeS, 1, mxREAL);
	SS = mxCreateDoubleMatrix(sizeS, 1, mxREAL);
	BB = mxCreateDoubleMatrix(sizeA, 1, mxREAL);
	DimM = mxCreateDoubleMatrix(1, 1, mxREAL);
	DimN = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *pII = mxGetPr(II), *pJJ = mxGetPr(JJ), 
		   *pSS = mxGetPr(SS), *pBB = mxGetPr(BB),
		   *pDimM = mxGetPr(DimM), *pDimN = mxGetPr(DimN);
	
	pDimM[0] = sizeA; pDimN[0] = sizeF;

	for (int n = 0; n < sizeS; ++n)
	{
		pII[n] = vII.at(n);
		pJJ[n] = vJJ.at(n);
		pSS[n] = vSS.at(n);
	}

	for (int j = 0; j < sizeA; ++j)
	{
		pBB[j] = b[j];
	}

	engPutVariable(ep, "II", II);
	engPutVariable(ep, "JJ", JJ);
	engPutVariable(ep, "SS", SS);
	engPutVariable(ep, "DimM", DimM);
	engPutVariable(ep, "DimN", DimN);
	engPutVariable(ep, "BB", BB);

	engEvalString(ep, "evals = scgls(II, JJ, SS, DimM, DimN, BB);");

	mxArray *evals = engGetVariable(ep, "evals");
	double *pEvals = mxGetPr(evals);
	for (int j = 0; j < sizeF; ++j)
	{
		sol[j] = pEvals[j];
	}

	mxDestroyArray(II);
	mxDestroyArray(JJ);
	mxDestroyArray(SS);
	mxDestroyArray(DimM);
	mxDestroyArray(DimN);
	mxDestroyArray(BB);
}

void WaveletMeshProcessor::computeMexicanHatWavelet( std::vector<double>& vMHW, double scale, int wtype /*= 1*/ )
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



void WaveletMeshProcessor::computeExperimentalWavelet( std::vector<double>& vExp, double scale )
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


inline double TransferFunc1(double x)
{
	if (x < 1) return x;
	else if (1 <= x && x <= 2) return (-5 + 11*x - 6*x*x + x*x*x);
	else return 2/x;
}

void WaveletMeshProcessor::calGeometryDWT()
{
// output: 1. coordinates of all vertices
//         2. wavelets coefficients of the x-coordinates

	ofstream ofs("output/coord.dat");
	for (int i = 0; i < m_size; ++i)
		ofs << mesh->getVertex(i)->getPos().x << ' ';
	ofs << endl;
	for (int i = 0; i < m_size; ++i)
		ofs << mesh->getVertex(i)->getPos().y << ' ';
	ofs << endl;
	for (int i = 0; i < m_size; ++i)
		ofs << mesh->getVertex(i)->getPos().z << ' ';
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

void WaveletMeshProcessor::reconstructExperimental1( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint /*= false*/ ) const
{
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

void WaveletMeshProcessor::reconstructByMHB( int aN, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz ) const
{
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
	int approxN = std::min(aN, mhb.m_nEigFunc);
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

void WaveletMeshProcessor::reconstructByDifferential( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint /* =false */ ) const
{
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
	mLaplacian.getSparseLaplacian(vII, vJJ, vSS);

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

void WaveletMeshProcessor::reconstructBySGW( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz, bool withConstraint /*= false*/ )
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
 		t_scales.push_back(40);
 		t_scales.push_back(20);
 		t_scales.push_back(10);
		computeSGW(t_scales);
	}
 	
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

	double weightI = 1;
	if (withConstraint)
	{
		SGW.push_back(vector<double>());
		SGW.back().resize(m_size, 0.0);
		SGW.back().at(pRef) = weightI;
		
		vxCoeff.push_back(posRef.x * weightI);
		vyCoeff.push_back(posRef.y * weightI);
		vzCoeff.push_back(posRef.z * weightI);
	}

// 	matlab_cgls(m_ep, SGW, vxCoeff, vx);
// 	matlab_cgls(m_ep, SGW, vyCoeff, vy);
// 	matlab_cgls(m_ep, SGW, vzCoeff, vz);

	matlab_scgls(m_ep, SGW, vxCoeff, vx);
	matlab_scgls(m_ep, SGW, vyCoeff, vy);
	matlab_scgls(m_ep, SGW, vzCoeff, vz);	
}

void WaveletMeshProcessor::computeSGW( const std::vector<double>& timescales )
{
	m_vTimescales = timescales;
	int nScales = m_vTimescales.size();

	m_vSGW.clear();
	//	scaling function
	if (true)	// must use scaling functions
	{
		for (int x = 0; x < m_size; ++x)
		{
			m_vSGW.push_back(vector<double>());
			m_vSGW.back().resize(m_size);

			for (int y = 0; y < m_size; ++y)
			{
				double itemSum = 0;
				for (int k = 0; k < mhb.m_nEigFunc; ++k)
				{
					itemSum += /*exp(-1.0) * */ exp(-mhb.m_func[k].m_val) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
				}
				m_vSGW.back().at(y) = itemSum;
			}
		}
	}	

	//	wavelet functions
	for (int s = 0; s < nScales; ++s)
	{
		for (int x = 0; x < m_size; ++x)
		{
			m_vSGW.push_back(vector<double>());
			m_vSGW.back().resize(m_size);

			for (int y = 0; y < m_size; ++y)
			{
				double itemSum = 0;
				for (int k = 0; k < mhb.m_nEigFunc; ++k)
				{
//					double coef = mhb.m_func[k].m_val * m_vTimescales[s];
					double coef = pow(mhb.m_func[k].m_val * m_vTimescales[s], 2.0);
					itemSum += coef * exp(-coef) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
				}
				m_vSGW.back().at(y) = itemSum;
			}
		}
	}
}

void WaveletMeshProcessor::getSGWSignature( double timescale, vector<double>& values ) const
{
	values.resize(m_size);

	for (int y = 0; y < m_size; ++y)
	{
		double itemSum = 0;
		for (int k = 0; k < mhb.m_nEigFunc; ++k)
		{
			double coef = mhb.m_func[k].m_val * timescale;
			itemSum += coef * exp(-coef) * mhb.m_func[k].m_vec[y] * mhb.m_func[k].m_vec[y];
		}
		values[y] = itemSum;
	}
}

void WaveletMeshProcessor::filterBySGW( std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz )
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

void WaveletMeshProcessor::deform( const std::vector<int>& vHandleIdx, const std::vector<Vector3D>& vHandelPos, const std::vector<int>& vFreeIdx, std::vector<Vector3D>& vDeformedPos, DeformType dfType )
{
	vDeformedPos.clear();

	if (dfType == Simple)
	{
		int hIdx = vHandleIdx[0];
		Vector3D handleTrans = vHandelPos[0] - mesh->getVertex_const(hIdx)->getPos();
		
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
			Vector3D newPos = mesh->getVertex_const(vFreeIdx[i])->getPos() + handleTrans * (1.0 - vDist2Handle[i]/distMax);
			vDeformedPos.push_back(newPos);
		}
		
	}


}
