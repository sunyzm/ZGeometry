#include "Laplacian.h"
#include <ctime>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <ZGeom/arithmetic.h>

using namespace std;
using ZGeom::PI;

const std::string MeshLaplacian::LaplacianTypeNames[] = {"Umbrella", "CotFormula", 
                                                         "Anisotropic", "Anisotropic2", 
                                                         "IsoApproximate"};

void ManifoldHarmonics::write( const std::string& meshPath, bool binaryMode /*= true*/ ) const
{
	if (binaryMode)
	{
		ofstream ofs(meshPath.c_str(), ios::trunc|ios::binary);
		ofs.write((char*)&m_nEigFunc, sizeof(int));
		ofs.write((char*)&m_size, sizeof(int));
		double *buf = new double[m_nEigFunc*(m_size+1)];
		for (int i = 0; i < m_nEigFunc; ++i)
		{
			buf[i * (m_size+1)] = m_func[i].m_val;
			for (int j = 0; j < m_size; ++j)
			{
				buf[i*(m_size+1)+j+1] = m_func[i].m_vec[j];
			}
		}
		ofs.write((char*)buf, sizeof(double)*m_nEigFunc*(m_size+1));
		delete []buf;
		ofs.close();
	}
	else
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
	}
	
	cout << "MHB saved to " << meshPath << endl;
}

void ManifoldHarmonics::read( const std::string& meshPath, bool binaryMode /*= true*/ )
{
	if (binaryMode)
	{
		ifstream ifs(meshPath.c_str(), ios::binary);
		ifs.read((char*)&m_nEigFunc, sizeof(int));
		ifs.read((char*)&m_size, sizeof(int));
		double *buf = new double[m_nEigFunc * (m_size+1)];
		ifs.read((char*)buf,sizeof(double)*m_nEigFunc*(m_size+1));

		m_func.resize(m_nEigFunc);
		for (int i = 0; i < m_nEigFunc; ++i)
		{
			m_func[i].m_vec.resize(m_size);
			m_func[i].m_val = buf[i*(m_size+1)];
			for (int j = 0; j < m_size; ++j)
				m_func[i].m_vec[j] = buf[i*(m_size+1)+j+1];
		}

		delete []buf;
	}
	else
	{
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
	}
	
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

/*
class ManifoldLaplaceHarmonics : public ManifoldHarmonics
{
public:
    bool decompLaplacian(Engine *ep, const CMesh *tmesh, int nEigFunc);
};

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
*/

void MeshLaplacian::decompose( ManifoldHarmonics& mhb, int nEig, Engine *ep ) const
{
    assert(mLaplacianConsructed);
    assert(nEig > 0);

    std::vector<int> vII, vJJ;
    std::vector<double> vSS, vWeights;
    mLS.convertToCOO(vII, vJJ, vSS, ZGeom::MAT_FULL);
    mW.getDiagonal(vWeights);

    mhb.m_func.clear();
    mhb.m_size = this->mOrder;
    mhb.m_nEigFunc = min(mhb.m_size, nEig);

    mxArray *II, *JJ, *SS, *AA, *evecs, *evals, *NUMV;

    AA = mxCreateDoubleMatrix(mOrder, 1, mxREAL);
    double *aa = mxGetPr(AA);
    assert((int)vWeights.size() == mhb.m_size);
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
        mhb.m_func[i].m_vec.reserve(mOrder);
        for(int j = 0; j < mOrder; j++)
        {
            mhb.m_func[i].m_vec.push_back(evec[i*mOrder+j]);
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

void MeshLaplacian::constructFromMesh1( const CMesh* tmesh )
{
	mOrder = tmesh->getVerticesNum();

    std::vector<unsigned> vII, vJJ;
    std::vector<double> vSS;

	for (int i = 0; i < mOrder; ++i)
	{
		const CVertex* vi = tmesh->getVertex_const(i);
		vector<int> vNeighbors;
		tmesh->VertexNeighborRing(i, 1, vNeighbors);
		int valence = vNeighbors.size();

		for (int j = 0; j < valence; ++j)
		{
			vII.push_back(i+1);
			vJJ.push_back(vNeighbors[j]+1);
			vSS.push_back(1.0);
		}
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(-valence);
	}

	std::vector<double> vWeights(mOrder, 1.0);

    mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
    mW.convertFromDiagonal(vWeights);

	mLaplacianConsructed = true;
	m_laplacianType = Umbrella;
}

void MeshLaplacian::constructFromMesh2( const CMesh* tmesh )
{
    mOrder = tmesh->getVerticesNum();

    std::vector<int> vII, vJJ;
    std::vector<double> vSS;
	std::vector<double> vWeights;

	tmesh->calLBO(vII, vJJ, vSS, vWeights);
/*
	double scaling = (tmesh->getAvgEdgeLength() * tmesh->getAvgEdgeLength())/2;
	transform(vWeights.begin(), vWeights.end(), vWeights.begin(), [&](double v){return v/scaling;});
*/
    mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
    mW.convertFromDiagonal(vWeights);

	mLaplacianConsructed = true;
	m_laplacianType = CotFormula;
}

void MeshLaplacian::constructFromMesh3( const CMesh* tmesh, int ringT, double hPara1, double hPara2 )
{
	// curvature difference based
    mOrder = tmesh->getVerticesNum();

    std::vector<int> vII, vJJ;
    std::vector<double> vSS;
    std::vector<double> vWeights;

	vector<std::tuple<int,int,double> > vSparseElements;

	ringT = 1;
	hPara1 = std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);

	for (int vi = 0; vi < mOrder; ++vi)
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
//				double w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
				double w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);

				double svalue = w1 * w2;
				svalue *= face_area;

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(mOrder, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < mOrder; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(mOrder, 1.0);	

    mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
    mW.convertFromDiagonal(vWeights);

	mLaplacianConsructed = true;
}

void MeshLaplacian::constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2)
{
	// similar to bilateral filtering
    mOrder = tmesh->getVerticesNum();

    std::vector<int> vII, vJJ;
    std::vector<double> vSS;
    std::vector<double> vWeights;

	vector<std::tuple<int,int,double> > vSparseElements;

	ringT = 1;
	hPara1 = std::pow(tmesh->getAvgEdgeLength(), 2);
//	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = tmesh->getAvgEdgeLength();

	for (int vi = 0; vi < mOrder; ++vi)
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
//				w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
				w2 = std::exp(dotProduct3D(pvi->getNormal(), pvk->getPosition() - pvi->getPosition()) / hPara2);
				
				double svalue = w1 * w2;
				svalue *= face_area;

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(mOrder, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < mOrder; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(mOrder, 1.0);	
    mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
    mW.convertFromDiagonal(vWeights);
	mLaplacianConsructed = true;
}

void MeshLaplacian::constructFromMesh5( const CMesh* tmesh )
{
    mOrder = tmesh->getVerticesNum();

    std::vector<int> vII, vJJ;
    std::vector<double> vSS;
    std::vector<double> vWeights;
	vector<std::tuple<int,int,double> > vSparseElements;

	int ringT = 5;
	double hPara1 = std::pow(tmesh->getAvgEdgeLength()*2, 2);
	
	for (int vi = 0; vi < mOrder; ++vi)
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

	vector<double> vDiag(mOrder, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < mOrder; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(mOrder, 1.0);	

    mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
    mW.convertFromDiagonal(vWeights);
	mLaplacianConsructed = true;
}