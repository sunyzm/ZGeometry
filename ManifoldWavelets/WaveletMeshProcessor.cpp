#include "WaveletMeshProcessor.h"
#include <fstream>

using namespace std;

void WaveletMeshProcessor::computeMexicanHatWavelet( std::vector<double>& vMHW, int scale, int wtype /*= 1*/ )
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

void WaveletMeshProcessor::computeExperimentalWavelet( std::vector<double>& vExp, int scale )
{
	vExp.resize(m_size);

	for (int i = 0; i < m_size; ++i)
	{
		double sum = 0;
		for (int k = 0; k < mhb.m_nEigFunc; ++k)
		{
			double coef = mhb.m_func[k].m_val * scale;
			////mexican hat
			sum += mhb.m_func[k].m_val  * exp(-mhb.m_func[k].m_val * scale) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[pRef];
			//// haar
			//sum += pow(1.0-exp(-mhb.m_func[k].m_val * scale), 2.0) / (mhb.m_func[k].m_val * scale) *  mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[pRef];
			//// Hermitian 1,2
			//sum += coef * exp(-coef * coef) * mhb.m_func[k].m_vec[i] * mhb.m_func[k].m_vec[pRef];
		}
		vExp[i] = sum;
	}

	ofstream fout("output/frequency.txt", ios::app);
	for (int k = 0; k < mhb.m_nEigFunc; ++k)
	{
		double coef = pow(mhb.m_func[k].m_val * scale, 2.0);
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
