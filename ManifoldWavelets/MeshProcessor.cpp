#include "MeshProcessor.h"
#include <fstream>

using namespace std;

MeshProcessor::MeshProcessor(void)
{
	ep = NULL;
	mesh = NULL;
	isMHBuilt = false;
	pRef = -1;
	m_size = 0;
}


MeshProcessor::~MeshProcessor(void)
{

}

void MeshProcessor::init(CMesh* tm, Engine* e)
{
	mesh = tm;
	ep = e;
	m_size = mesh->getVerticesNum();
	pRef = 0;
}

void MeshProcessor::decomposeLaplacian(int nEigFunc)
{
	if (mhb.decompLaplacian(this->ep, this->mesh, nEigFunc, LBO_COT))
		this->isMHBuilt = true;
}

void MeshProcessor::readMHB( std::string path )
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

}

void MeshProcessor::writeMHB(std::string path)
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

void MeshProcessor::computeMexicanHatWavelet( std::vector<double>& vMHW, int scale, int wtype /*= 1*/ )
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

void MeshProcessor::computeExperimentalWavelet( std::vector<double>& vExp, int scale )
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

void MeshProcessor::computeCurvature( std::vector<double>& vCurvature, int curvatureType /*= 0*/ )
{
	vCurvature.resize(m_size);
	if (curvatureType == 0)
	{
		for (int i = 0; i < m_size; ++i)
		{
//			mesh->calVertexCurvature(i);
			vCurvature[i] = mesh->getVertex(i)->m_vMeanCurvature;
		}
	}
	else if (curvatureType == 1)
	{
		for (int i = 0; i < m_size; ++i)
			vCurvature[i] = mesh->getVertex(i)->m_vGaussCurvature;
	}
	
}

void MeshProcessor::normalizeFrom(const std::vector<double>& vFrom)
{
	if (vFrom.empty()) return;	
	auto iResult = minmax_element(vFrom.begin(), vFrom.end());
	double sMin = *iResult.first, sMax = *iResult.second;

	this->vDisplaySignature.clear();
	for (vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
	}
}

void MeshProcessor::logNormalizeFrom( const std::vector<double>& vFrom )
{
	if (vFrom.empty()) return;	

	std::vector<double> vLog;
	vLog.reserve(vFrom.size());
	for (std::vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		vLog.push_back(std::log(*iter + 1));
	}

	auto iResult = minmax_element(vLog.begin(), vLog.end());
	double sMin = *iResult.first, sMax = *iResult.second;

	this->vDisplaySignature.clear();
	this->vDisplaySignature.reserve(vLog.size());
	for (std::vector<double>::const_iterator iter = vLog.begin(); iter != vLog.end(); ++iter)
	{
		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
	}
}

void MeshProcessor::bandCurveFrom( const std::vector<double>& vFrom, double lowend, double highend )
{
	assert(lowend < highend);
	this->vDisplaySignature.clear();
	this->vDisplaySignature.reserve(vFrom.size());
	for (std::vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		if (*iter <= lowend)
			vDisplaySignature.push_back(0.0);
		else if (*iter >=highend)
			vDisplaySignature.push_back(1.0);
		else 
			vDisplaySignature.push_back((*iter - lowend)/(highend - lowend));
	}
}

