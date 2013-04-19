#include <ctime>
#include <algorithm>
#include <cstdio>
#include <tuple>
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
	assert(m_bMatrixBuilt);
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
	this->m_size = tmesh->getVerticesNum();

	this->vII.clear();
	this->vJJ.clear();
	this->vSS.clear();
	this->vWeights.clear();

	tmesh->calLBO(vII, vJJ, vSS, vWeights);

	m_bMatrixBuilt = true;
}

bool ManifoldLaplaceHarmonics::decompLaplacian( Engine *ep, const CMesh *tmesh, int nEigFunc )
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
	
	for(int i = 0; i < nVertex; i++)	//for each vertex
	{
		double Av;
		tmesh->calVertexLBO(i, IIv, JJv, SSv, Av, diagW);
		aa[i] = Av;		//mixed area
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

void MeshLaplacian::constructFromMesh1( const CMesh* tmesh )
{
	this->m_size = tmesh->getVerticesNum();

	vII.clear();
	vJJ.clear();
	vSS.clear();
	vWeights.clear();

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

	m_bMatrixBuilt = true;
	m_laplacianType = Umbrella;
}

void MeshLaplacian::constructFromMesh2( const CMesh* tmesh )
{
	this->m_size = tmesh->getVerticesNum();

	vII.clear();
	vJJ.clear();
	vSS.clear();
	vWeights.clear();

	tmesh->calLBO(vII, vJJ, vSS, vWeights);

	m_bMatrixBuilt = true;
	m_laplacianType = CotFormula;
}

void MeshLaplacian::constructFromMesh3( const CMesh* tmesh, int ringT, double hPara1, double hPara2 )
{
	// curvature difference based

	vII.clear();
	vJJ.clear();
	vSS.clear();

	this->m_size = tmesh->getMeshSize();
	vector<std::tuple<int,int,double> > vSparseElements;

	for (int vi = 0; vi < m_size; ++vi)
	{
		const CVertex* pvi = tmesh->getVertex_const(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFacesIndex(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi)
		{
			const CFace* pfi = tmesh->getFace_const(vFaces[fi]);
			double face_area = pfi->getArea();
			for (int k = 0; k < 3; ++k)
			{
				int vki = pfi->getVertexIndex(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->getVertex_const(k);

//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				double w1 = std::exp(-(pvi->getPosition()-pvk->getPosition()).length2() / hPara1);
				double w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
//				double w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);

				double svalue = w1 * w2;
				svalue *= face_area;

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(m_size, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < m_size; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(m_size, 1.0);	

	m_bMatrixBuilt = true;
}

void MeshLaplacian::constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2)
{
	// similar to bilateral filtering
	
	vII.clear();
	vJJ.clear();
	vSS.clear();

	this->m_size = tmesh->getMeshSize();
	vector<std::tuple<int,int,double> > vSparseElements;

	ringT = 5;
	hPara1 = std::pow(tmesh->getAvgEdgeLength()*5, 2);
	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);

	for (int vi = 0; vi < m_size; ++vi)
	{
		const CVertex* pvi = tmesh->getVertex_const(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFacesIndex(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi)
		{
			const CFace* pfi = tmesh->getFace_const(vFaces[fi]);
			double face_area = pfi->getArea();
			for (int k = 0; k < 3; ++k)
			{
				int vki = pfi->getVertexIndex(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->getVertex_const(k);

				double w1 = 1., w2 = 1.;
//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
//				w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
				w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
//				w2 = std::exp((dotProduct3D(pvi->getNormal(), pfi->getNormal()) - 1) / 1.0);


				double svalue = w1 * w2;
				svalue *= face_area;

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(m_size, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < m_size; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(m_size, 1.0);	

	m_bMatrixBuilt = true;
}


void MeshLaplacian::constructFromMesh5( const CMesh* tmesh )
{
	vII.clear();
	vJJ.clear();
	vSS.clear();

	this->m_size = tmesh->getMeshSize();
	vector<std::tuple<int,int,double> > vSparseElements;

	int ringT = 5;
	double hPara1 = std::pow(tmesh->getAvgEdgeLength()*2, 2);
	
	for (int vi = 0; vi < m_size; ++vi)
	{
		const CVertex* pvi = tmesh->getVertex_const(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFacesIndex(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi)
		{
			const CFace* pfi = tmesh->getFace_const(vFaces[fi]);
			double face_area = pfi->getArea();
			for (int k = 0; k < 3; ++k)
			{
				int vki = pfi->getVertexIndex(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->getVertex_const(k);

				double svalue = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / (4 * hPara1));
				svalue = svalue * face_area / (3 * 4 * PI * hPara1 * hPara1);

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(m_size, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < m_size; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(m_size, 1.0);	

	m_bMatrixBuilt = true;

}