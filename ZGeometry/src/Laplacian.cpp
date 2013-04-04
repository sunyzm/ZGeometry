#include <ctime>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "Laplacian.h"

using namespace std;

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

MeshFunction ManifoldHarmonics::getManifoldHarmonic( int k ) const
{
	if (k > m_nEigFunc)
		throw runtime_error("Invalid request for manifold harmonic");
	MeshFunction mf(m_size);
	for (int i = 0; i < m_size; ++i)
	{
		mf[i] = m_func[k].m_vec[i];
	}
	return mf;
}

void ManifoldHarmonics::dumpEigenValues( const std::string& pathEVL ) const
{
	ofstream ofs(pathEVL.c_str(), ios::trunc);
	for (auto iter = m_func.begin(); iter != m_func.end(); ++iter)
		ofs << iter->m_val << endl;
	ofs.close();
}

//: computer inner product of two manifold functions induced by the Laplacian matrix
double SparseMeshMatrix::innerProduct( const std::vector<double>& vf, const std::vector<double>& vg ) const
{
	assert(vf.size() == vg.size() && vf.size() == vWeights.size() && (int)vWeights.size() == this->m_size);

	double sum = 0;
	for (int i = 0; i < m_size; ++i)
	{
		sum += vf[i] * vWeights[i] * vg[i];
	}
	return sum;
}

void SparseMeshMatrix::multiply( Engine *ep, const std::vector<double>& func, std::vector<double>& result ) const
{
	assert(func.size() == m_size);
	/*	
	mxArray *II, *JJ, *SS, *AA, *evecs, *evals;

	int ns = (int) vII.size();
	II = mxCreateDoubleMatrix(ns, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(ns, 1, mxREAL);
	SS = mxCreateDoubleMatrix(ns, 1, mxREAL);
	double *ii = mxGetPr(II);
	double *jj = mxGetPr(JJ);
	double *ss = mxGetPr(SS);
	std::copy(vII.begin(), vII.end(), ii);
	std::copy(vJJ.begin(), vJJ.end(), jj);
	std::copy(vSS.begin(), vSS.end(), ss);

	engPutVariable(ep, "II", II);
	engPutVariable(ep, "JJ", JJ);
	engPutVariable(ep, "SS", SS);

	//TODO: to finish the multiplication
	*/
}

void SparseMeshMatrix::getSparseMatrix( std::vector<int>& II, std::vector<int>& JJ, std::vector<double>& SS ) const
{
	II = vII;
	JJ = vJJ;
	SS = vSS;
	int nz = II.size();

	for (int k = 0; k < nz; ++k)
	{
		SS[k] /= vWeights[II[k]];
	}
}

void SparseMeshMatrix::dumpMatrix( const std::string& path ) const
{
	int nz = vII.size();
	ofstream lout(path.c_str());
	for (int i = 0; i < nz; ++i)
	{
		lout << vII[i] << ',' << vJJ[i] << ',' << vSS[i] << '\n';
	}
	lout.close();
}

void SparseMeshMatrix::decompose( ManifoldHarmonics& mhb, int nEig, Engine *ep ) const
{
	assert(isBuilt);
	assert(nEig > 0);

	mhb.m_func.clear();
	mhb.m_size = this->m_size;
	mhb.m_nEigFunc = min(mhb.m_size, nEig);

	mxArray *II, *JJ, *SS, *AA, *evecs, *evals, *NUMV;

	AA = mxCreateDoubleMatrix(m_size, 1, mxREAL);
	double *aa = mxGetPr(AA);
	assert((int)vWeights.size() == m_size);
	std::copy(vWeights.begin(), vWeights.end(), aa);

	NUMV = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *numv = mxGetPr(NUMV);	
	numv[0] = mhb.m_nEigFunc;			// number of eigen vectors to be computed

	int ns = (int) vII.size();
	II = mxCreateDoubleMatrix(ns, 1, mxREAL);
	JJ = mxCreateDoubleMatrix(ns, 1, mxREAL);
	SS = mxCreateDoubleMatrix(ns, 1, mxREAL);
	double *ii = mxGetPr(II);
	double *jj = mxGetPr(JJ);
	double *ss = mxGetPr(SS);
	std::copy(vII.begin(), vII.end(), ii);
	std::copy(vJJ.begin(), vJJ.end(), jj);
	std::copy(vSS.begin(), vSS.end(), ss);

	engPutVariable(ep, "II", II);
	engPutVariable(ep, "JJ", JJ);
	engPutVariable(ep, "SS", SS);
	engPutVariable(ep, "AA", AA);
	engPutVariable(ep, "Numv", NUMV);

	engEvalString(ep, "[evecs,evals] = hspeigs(II,JJ,SS,AA,Numv);");

	evecs = engGetVariable(ep, "evecs");		
	double *evec = mxGetPr(evecs);				//eigenvectors
	evals = engGetVariable(ep, "evals");		
	double *eval = mxGetPr(evals);				//eigenvalues

	mhb.m_func.reserve(mhb.m_nEigFunc);
	for(int i = 0; i < mhb.m_nEigFunc; i++)
	{
		mhb.m_func.push_back(ManifoldBasis());
		mhb.m_func[i].m_vec.reserve(m_size);
		for(int j = 0; j < m_size; j++)
		{
			mhb.m_func[i].m_vec.push_back(evec[i*m_size+j]);
		}
		mhb.m_func[i].m_val = std::fabs(eval[i]);		// always non-negative
	}

	mxDestroyArray(evecs);
	mxDestroyArray(evals);
	mxDestroyArray(AA);
	mxDestroyArray(II);
	mxDestroyArray(JJ);
	mxDestroyArray(SS);
	mxDestroyArray(NUMV);
}

void Laplacian::constructFromMesh( const CMesh* tmesh )
{
//	m_laplacianType = CotFormula;	// already set in contructors

	this->m_size = tmesh->getVerticesNum();
	
	vector<double> diagW;
	diagW.resize(m_size, 0);

	vII.clear();
	vJJ.clear();
	vSS.clear();
	vWeights.clear();

	if(m_laplacianType == CotFormula)
	{
// 		for(int k = 0; k < m_size; k++)	//for each vertex
// 		{
// 			double Av;
// 			tmesh->calVertexLBO2(k, vII, vJJ, vSS, Av, diagW);
// 			vWeights.push_back(Av);		//mixed area as weight of each vertex
// 		}
// 		for(int k = 0; k < m_size; k++)
// 		{
// 			vII.push_back(k+1);
// 			vJJ.push_back(k+1);
// 			vSS.push_back(diagW[k]);
// 		}

		tmesh->calLBO(vII, vJJ, vSS, vWeights);
	}
	else if (m_laplacianType == Umbrella)
	{
		for (int i = 0; i < m_size; ++i)
		{
			const CVertex* vi = tmesh->getVertex_const(i);
			vector<int> vNeighbors;
			tmesh->VertexNeighborRing(i, 1, vNeighbors);
			int valence = vNeighbors.size();

			for (int j = 0; j < valence; ++j)
			{
				vII.push_back(i+1);
				vJJ.push_back(vNeighbors[j]+1);
				vSS.push_back(-1.0);
			}
			vII.push_back(i+1);
			vJJ.push_back(i+1);
			vSS.push_back(valence);
		}

		vWeights.resize(m_size, 1.0);
	}

	isBuilt = true;
}

bool ManifoldLaplaceHarmonics::decompLaplacian( Engine *ep, const CMesh *tmesh, int nEigFunc, Laplacian::LaplacianType lbo_type /*= Laplacian::CotFormula*/ )
{
	m_func.clear();

	const int nVertex = tmesh->getVerticesNum();
	m_size = nVertex;
	m_nEigFunc = min(m_size, nEigFunc);

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

	if(lbo_type == Laplacian::CotFormula)
	{
		for(int i = 0; i < nVertex; i++)	//for each vertex
		{
			double Av;
			tmesh->calVertexLBO(i, IIv, JJv, SSv, Av, diagW);
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
	std::copy(IIv.begin(), IIv.end(), ii);
	std::copy(JJv.begin(), JJv.end(), jj);
	std::copy(SSv.begin(), SSv.end(), ss);

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

	return true;
}

void AnisotropicKernel::constructFromMesh( const CMesh* tmesh )
{
	for (int i = 0; i < m_size; ++i)
	{
		const CVertex* vi = tmesh->getVertex_const(i);
		vector<int> vNeighbors;
		tmesh->VertexNeighborRing(i, 3, vNeighbors);
		int valence = vNeighbors.size();

		vector<double> vDist;
		double avg_len = tmesh->getAvgEdgeLength();
		double dist_sum = 0.;

		for (int j = 0; j < valence; ++j)
		{
			double coeff = exp( -tmesh->getGeodesic(i, vNeighbors[j]) / avg_len);
			vDist.push_back(coeff);
			dist_sum += coeff;
		}
		
		for (int j = 0; j < valence; ++j)
		{
			vII.push_back(i+1);
			vJJ.push_back(vNeighbors[j]+1);
			vSS.push_back(-vDist[j]/dist_sum);
		}

		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(1.);
	}

	vWeights.resize(m_size, 1.0);
}
