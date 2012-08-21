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

void MeshProcessor::normalizeFrom(const std::vector<double>& vFrom)
{
	if (vFrom.empty()) return;	
	auto iResult = minmax_element(vFrom.begin(), vFrom.end());
	vector<double>::const_iterator iMin = iResult.first, iMax = iResult.second;
	double sMin = *iMin, sMax = *iMax;

	this->vDisplaySignature.clear();
	for (vector<double>::const_iterator iter = vFrom.begin(); iter != vFrom.end(); ++iter)
	{
		vDisplaySignature.push_back((*iter - sMin)/(sMax - sMin));
	}
}

