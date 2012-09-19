#include "ManifoldMeshProcessor.h"
#include <fstream>

using namespace std;

void ManifoldMeshProcessor::computeMexicanHatWavelet( std::vector<double>& vMHW, int scale, int wtype /*= 1*/ )
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

void ManifoldMeshProcessor::computeExperimentalWavelet( std::vector<double>& vExp, int scale )
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

void ManifoldMeshProcessor::calGeometryDWT()
{
// output: 1. x-coordinates of all vertices
//         2. wavelets coefficients of the x-coordinates
	ofstream ofs("output/coordinates.txt");
	for (int i = 0; i < m_size; ++i)
	{
		ofs << mesh->getVertex(i)->getPos().x << ' ';
	}
	ofs.close();

	ofs.open("output/wavelets.txt");

	double scale = 10;
	for (int x = 0; x < m_size; ++x)
	{
		for (int y = 0; y < m_size; ++y)
		{
			double itemSum = 0;
			for (int k = 0; k < mhb.m_nEigFunc; ++k)
			{
				itemSum += pow(mhb.m_func[k].m_val * scale, 2) * exp(-pow(mhb.m_func[k].m_val * scale,2)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
			}
			ofs << itemSum << ' ';
		}
		ofs << endl;

		for (int y = 0; y < m_size; ++y)
		{
			double itemSum = 0;
			for (int k = 0; k < mhb.m_nEigFunc; ++k)
			{
				itemSum += exp(-1.0) * exp(-pow(mhb.m_func[k].m_val * scale / 3, 4)) * mhb.m_func[k].m_vec[x] * mhb.m_func[k].m_vec[y];
			}
			ofs << itemSum << ' ';
		}
		ofs << endl;

	}
}
