#include <ctime>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include "Laplacian.h"

using namespace std;

void ManifoldHarmonics::decompLaplacian( Engine *ep, const CMesh *tmesh, int nEigFunc, short lbo_type /*= LBO_COT*/ )
{
	clock_t start, end;
	start = clock();

	m_func.clear();

	const int nVertex = tmesh->getVerticesNum();
	m_size = nVertex;
	m_nEigFunc = std::min(m_size, nEigFunc);
	
	mxArray *II, *JJ, *SS, *AA, *evecs, *evals, *Numv;
	AA = mxCreateDoubleMatrix(nVertex, 1, mxREAL);
	double *aa = mxGetPr(AA);
	Numv = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *numv = mxGetPr(Numv);	
	numv[0] = m_nEigFunc;			// number of eigen vectors to be computed
	vector<int> IIv;
	vector<int> JJv;
	vector<double> SSv;
	vector<double> diagW;
	diagW.resize(nVertex, 0);

	if(lbo_type == LBO_COT)
	{
		for(int i = 0; i < nVertex; i++)	//for each vertex
		{
			double Av;
			tmesh->CalVertexLBO(i, IIv, JJv, SSv, Av, diagW);
			aa[i] = Av;		//mixed area
		}
	}
	for(int i = 0; i < nVertex; i++)
	{
		IIv.push_back(i+1);
		JJv.push_back(i+1);
		SSv.push_back(diagW[i]);
	}

	int ns = (int) IIv.size();
	II = mxCreateDoubleMatrix(ns, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(ns, 1, mxREAL);
	SS = mxCreateDoubleMatrix(ns, 1, mxREAL);
	double *ii = mxGetPr(II);
	double *jj = mxGetPr(JJ);
	double *ss = mxGetPr(SS);

	for(int k = 0; k < nVertex; k++)
	{
		IIv[k] = ii[k];
		JJv[k] = jj[k];
		SSv[k] = ss[k];
	}
//	std::copy(IIv.begin(), IIv.end(), ii);
//	std::copy(JJv.begin(), JJv.end(), jj);
//	std::copy(SSv.begin(), SSv.end(), ss);

	engPutVariable(ep, "II", II);
	engPutVariable(ep, "JJ", JJ);
	engPutVariable(ep, "SS", SS);
	engPutVariable(ep, "AA", AA);
	engPutVariable(ep, "Numv", Numv);

	engEvalString(ep, "[evecs,evals] = hspeigs(II,JJ,SS,AA,Numv);");

	evecs = engGetVariable(ep, "evecs");		
	double *evec = mxGetPr(evecs);				//eigenvectors
	evals = engGetVariable(ep, "evals");		
	double *eval = mxGetPr(evals);				//eigenvalues

	m_func.reserve(m_nEigFunc);
	for(int i = 0; i < m_nEigFunc; i++)
	{
		m_func.push_back(ManifoldBasis());

		m_func[i].m_vec.reserve(nVertex);
		for(int j = 0; j < nVertex; j++)
		{
			m_func[i].m_vec.push_back(evec[i*nVertex+j]);
		}
		m_func[i].m_val = abs(eval[i]);		// always non-negative
	}

	mxDestroyArray(evecs);
	mxDestroyArray(evals);
	mxDestroyArray(AA);
	mxDestroyArray(II);
	mxDestroyArray(JJ);
	mxDestroyArray(SS);
	mxDestroyArray(Numv);

	end = clock();
	double tt = double(end-start) / CLOCKS_PER_SEC;
	cout << "--Decomposition time: " << tt << endl;
}

void ManifoldHarmonics::write( const std::string& meshPath ) const
{
	ofstream ofs(meshPath.c_str(), ios::trunc);
	ofs << m_nEigFunc << endl;
	ofs << m_size << endl;
	for (int i = 0; i < m_nEigFunc; ++i)
	{
		ofs << m_func[i].m_val;
		for (int j = 0; j < m_size; ++j)
			ofs << ' ' << m_func[i].m_vec[j];
		ofs << endl;
	}
	ofs.close();
	cout << "MHB saved to " << meshPath << endl;
}

void ManifoldHarmonics::read( const std::string& meshPath )
{
	// 	FILE* file = NULL;
	// 	fopen_s(&file, meshPath.c_str(), "rb");
	// 	int MAX = 536870912;
	// 	char *buffer = new char[MAX];
	// 	setvbuf(file, NULL, _IOFBF, 1024);
	// 	while (!feof(file)) 
	// 	{
	// 		fread(buffer, 0x1, MAX, file);
	// 	}
	// 	fclose(file);

	ifstream ifs(meshPath.c_str(), ios::binary);
	ifs.seekg (0, ios::end);
	int length = (int)ifs.tellg();
	ifs.seekg (0, ios::beg);
	char *buffer = new char[length];
	ifs.read(buffer, length);
	ifs.close();
	
	istringstream iss (buffer, istringstream::in);
	iss >> m_nEigFunc;
	iss >> m_size;
	m_func.resize(m_nEigFunc);
	for (int i = 0; i < m_nEigFunc; ++i)
	{
		m_func[i].m_vec.resize(m_size);
		iss >> m_func[i].m_val;
		for (int j = 0; j < m_size; ++j)
			iss >>m_func[i].m_vec[j];
	}
	delete []buffer;
/*
	ifstream ifs(meshPath.c_str());
	ifs >> m_nEigFunc;
	ifs >> m_size;
	m_func.resize(m_nEigFunc);
	
	for (int i = 0; i < m_nEigFunc; ++i)
	{
		m_func[i].m_vec.resize(m_size);
		ifs >> m_func[i].m_val;
		for (int j = 0; j < m_size; ++j)
			ifs >>m_func[i].m_vec[j];
	}
*/
	cout << "MHB loaded from " << meshPath << endl;
}
