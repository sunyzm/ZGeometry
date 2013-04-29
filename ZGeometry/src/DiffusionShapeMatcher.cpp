#include <fstream>
#include <sstream>
#include <exception>
#include <set>
#include <algorithm>
#include <cassert>
#include "DiffusionShapeMatcher.h"
#include "OutputHelper.h"
#include "SimpleConfigLoader.h"
#include <ppl.h>
#include <functional>
#include <limits>

#define LOCAL_ANCHORS_NUM 8

using namespace std;

extern OutputHelper qout;
extern SimpleConfigLoader g_configMgr;

const double DiffusionShapeMatcher::DEFAULT_C_RATIO				= 0.2;
const double DiffusionShapeMatcher::DEFAULT_RANK_EPSILON		= 1e-4;
const double DiffusionShapeMatcher::SPARSIFY_EPSILON			= 1e-6;
const double DiffusionShapeMatcher::DEFAULT_FEATURE_TIMESCALE	= 30.0;
const double DiffusionShapeMatcher::DEFAULT_T_MULTIPLIER		= 3.0;
const double DiffusionShapeMatcher::DEFAULT_MATCH_TIME_LOW		= 10.0;
const double DiffusionShapeMatcher::DEFAULT_MATCH_THRESH        = 0.52;
const double DiffusionShapeMatcher::DEFAULT_EXTREAMA_THRESH     = 0.04;
const double DiffusionShapeMatcher::DEFAULT_REGISTER_TIMESCALE	= 80;
const int	 DiffusionShapeMatcher::NUM_OF_EIGVAL_FOR_ESTIMATE	= 50;
const int	 DiffusionShapeMatcher::DEFAULT_PYRAMID_LEVELS		= 3;
const int	 DiffusionShapeMatcher::MAXIMAL_PYRAMID_LEVELS		= 5;

double findTmax( const CMesh* tmesh, int s )
{
	if (tmesh->getBoundaryVertexNum() == 0 ) return -1;
	double geo = tmesh->getGeodesicToBoundary(s) / tmesh->getAvgEdgeLength();	//normalized by average edge length; requires isHole initialized
	return geo*geo/4.0;
}

void PointParam::clear()
{
	m_size = 0; 
	m_votes = 0.0; 
	if(!empty()) 
		delete m_vec; 
	m_vec = NULL;
}

PointParam& PointParam::operator=( const PointParam& hkp )
{
	reserve(hkp.m_size);
	for (int i = 0; i < m_size; ++i)
		m_vec[i] = hkp.m_vec[i];
	m_votes = hkp.m_votes;
	return *this;
}

void ParamManager::computeHKParam( const std::vector<int>& anchors, double t /*= 30.0*/ )
{
	const int fineSize = pMP->getMesh_const()->getVerticesNum();
	const int pn = (int)anchors.size();

	vCoord.resize(fineSize);

 	for (int v = 0; v < fineSize; ++v)
 	{
 		PointParam& hkp = vCoord[v];
 		hkp.reserve(pn);
 		for (int i = 0; i < pn; ++i)
 			hkp.m_vec[i] = pMP->calHK(v, anchors[i], t);
 
 		hkp.m_votes = hkp.length();
 	}
}

void ParamManager::computeBHParam( const std::vector<int>& anchors )
{
	const int fineSize = pMP->getMesh_const()->getVerticesNum();
	const int pn = (int)anchors.size();
	vCoord.resize(fineSize);

	for (int v = 0; v < fineSize; ++v)
	{
		PointParam& hkp = vCoord[v];
		hkp.reserve(pn);
		for (int i = 0; i < pn; ++i)
			hkp.m_vec[i] = pMP->calBiharmonic(v, anchors[i]);

		hkp.m_votes = hkp.length();
	}
}

void ParamManager::para_computeHKC( const std::vector<int>& anchors, double t /*= 30.0*/ )
{
	const int fineSize = pMP->getMesh_const()->getVerticesNum();
	const int pn = (int)anchors.size();
	vCoord.resize(fineSize);

	Concurrency::parallel_for(0, fineSize, [&](int v){
	 	PointParam& hkp = vCoord[v];
	 	hkp.reserve(pn);
	 	for (int i = 0; i < pn; ++i)
	 		hkp.m_vec[i] = pMP->calHK(v, anchors[i], t);
	 	hkp.m_votes = hkp.length();
	});
}

void ParamManager::para_computeHKS( const std::vector<double>& times )
{
	const int fineSize = pMP->getMesh_const()->getVerticesNum();
	const int sn = (int)times.size();
	vSignature.resize(fineSize);

	Concurrency::parallel_for(0, fineSize, [&](int v){
		vSignature[v].reserve(sn);
		for (int i = 0; i < sn; ++i)
			vSignature[v].m_vec[i] = pMP->calHK(v, v, times[sn]);
	});
}

double distFeature(const HKSFeature& hf1, const HKSFeature& hf2, const VectorND& sig1, const VectorND& sig2, double& tl, int& tn)
{
	tl = max(hf1.m_tl, hf2.m_tl);	// now both are 10.0
	double tu = min(hf1.m_tu, hf2.m_tu);
	tn = int(std::log(tu/tl)/std::log(2.0) + 1.5);	//only consider the overlapping times

	if(tn == 0) return 1.0;

	double dist = 0.0;
	for(int i = 0; i < tn; i++)
	{
		dist += (sig1.m_vec[i] - sig2.m_vec[i]) * (sig1.m_vec[i] - sig2.m_vec[i]);
	}
	return dist/tn;
}

double distFeaturePair(const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const MatchPair& mp1, const MatchPair& mp2)
{
	double t = max(mp1.m_tl, mp2.m_tl); // take the bigger one
	int n = min(mp1.m_tn, mp2.m_tn);

	int x1 = mp1.m_idx1;
	int x2 = mp1.m_idx2;
	int y1 = mp2.m_idx1;
	int y2 = mp2.m_idx2;

	VectorND v1(n), v2(n);
	for(int i = 0; i < n; i++)
	{
		v1.m_vec[i] = pmp1->calHK(x1, y1, t);
		v2.m_vec[i] = pmp2->calHK(x2, y2, t);
		t *= 2.0;
	}

	return v1.calDistance2(v2) / (v1.length2() + v2.length2());
}

class ExtremaPoint
{
public:
	int index;
	double val;
	ExtremaPoint(int i, double v) : index(i), val(v) {}

	friend bool operator<(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val < ep2.val) return true;
		else return false;
	}

	friend bool operator<=(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val <= ep2.val) return true;
		else return false;
	}


	friend bool operator>(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val > ep2.val) return true;
		else return false;
	}

	friend bool operator>=(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val >= ep2.val) return true;
		else return false;
	}
};

DiffusionShapeMatcher::DiffusionShapeMatcher()
{
	m_ep = NULL;
	pOriginalProcessor[0] = pOriginalProcessor[1] = NULL; 
	pOriginalMesh[0] = pOriginalMesh[1] = NULL;
	m_bInitialized = false;
	m_bPyramidBuilt = false;
	m_nAlreadyRegisteredLevel = -1;
}

DiffusionShapeMatcher::~DiffusionShapeMatcher()
{
	if (!m_bInitialized) return;

	int mpSize1 = liteMP[0].size(), mpSize2 = liteMP[1].size();

	for (int k = 1; k < mpSize1; ++k)
	{
		delete liteMP[0][k];
	}
	for (int k = 1; k < mpSize2; ++k)
	{
		delete liteMP[1][k];
	}
}

void DiffusionShapeMatcher::initialize( DifferentialMeshProcessor* pMP1, DifferentialMeshProcessor* pMP2, Engine *ep )
{
	m_ep = ep;
	pOriginalProcessor[0] = pMP1;
	pOriginalProcessor[1] = pMP2;
	pOriginalMesh[0] = pMP1->getMesh();
	pOriginalMesh[1] = pMP2->getMesh();

	liteMP[0].push_back(pMP1);
	liteMP[1].push_back(pMP2);

	meshPyramids[0].setInitialMesh(pOriginalMesh[0]);
	meshPyramids[1].setInitialMesh(pOriginalMesh[1]);

	m_ParamMgr[0].initialize(pMP1);
	m_ParamMgr[1].initialize(pMP2);

	m_bInitialized = true;
}

CMesh* DiffusionShapeMatcher::getMesh( int obj, int level /*= 0*/ ) const
{
	if (level == 0)
		return this->pOriginalMesh[obj];
	else
		return meshPyramids[obj].getMesh(level);
}

void DiffusionShapeMatcher::constructPyramid( int n, double ratio, std::ostream& ostr )
{
	assert(n >= 1);

	m_nPyramidLevels = n;

	meshPyramids[0].setLevel(n, ratio);
	meshPyramids[1].setLevel(n, ratio);

	meshPyramids[0].construct(ostr);
	meshPyramids[1].construct(ostr);

	for (int k = 1; k < n; ++k)
	{
		liteMP[0].push_back(new DifferentialMeshProcessor(meshPyramids[0].getMesh(k)));
		liteMP[1].push_back(new DifferentialMeshProcessor(meshPyramids[1].getMesh(k)));
	}

	m_bPyramidBuilt = true;

	setRegistrationLevels(n);
}

void DiffusionShapeMatcher::detectFeatures( int obj, int ring /*= 2*/, int scale /*= 1*/, double baseTvalue /*= DEFAULT_FEATURE_TIMESCALE*/, double talpha /*= DEFAULT_T_MULTIPLIER*/, double thresh /*= DEFAULT_EXTREAMA_THRESH*/ )
{
	assert(obj == 0 || obj == 1);

	const CMesh* fineMesh = pOriginalMesh[obj];
	DifferentialMeshProcessor* pMP = pOriginalProcessor[obj];
	const int fineSize = fineMesh->getVerticesNum();

	vector<HKSFeature>& vF = vFeatures[obj];
	vF.clear();

	for (int s = scale-1; s >= 0; --s)
	{
		double tvalue = baseTvalue * std::pow(talpha, s);
		vector<double> hksv;
		pMP->calKernelSignature(tvalue, HEAT_KERNEL, hksv);

 		double sref = 4.0 * PI * tvalue;
 		transform( hksv.begin(), hksv.end(), hksv.begin(), [=](double v){ return std::log(v * sref);} );		

		vector<int> vFeatureIdx;
		if (s == 0) ring = 3;
		fineMesh->extractExtrema(hksv, ring, thresh, vFeatureIdx);
		cout << "obj " << obj << ", timescale: " << tvalue << endl;
		for (auto iter = vFeatureIdx.begin(); iter != vFeatureIdx.end(); ++iter)
		{
			cout << *iter << ", ";
		}
		cout << endl;

		priority_queue<ExtremaPoint, vector<ExtremaPoint>, greater<ExtremaPoint> > extremaPointsQueue;
		for (auto iter = vFeatureIdx.begin(); iter != vFeatureIdx.end(); ++iter)
			extremaPointsQueue.push(ExtremaPoint(*iter, std::abs(hksv[*iter])));
		
		while (!extremaPointsQueue.empty())
		{
			ExtremaPoint ep = extremaPointsQueue.top();
			extremaPointsQueue.pop();

			bool notCoveredFine = true;
			for (auto iter = vF.begin(); iter != vF.end(); ++iter)
			{
				if (fineMesh->isInNeighborRing(iter->m_index, ep.index, 3))
				{
					notCoveredFine = false;
					break;
				}
			}
			if (notCoveredFine)
			{
				vF.push_back(HKSFeature(ep.index, s));
			}
		}
	}

	/* Begin LOG-FEATURES */
	cout << "Candidate features " << obj << ":\n  ";
	std::sort(vF.begin(), vF.end(), [](HKSFeature& f1, HKSFeature& f2) { return f1.m_index < f2.m_index;});
	std::sort(vF.begin(), vF.end(), [](HKSFeature& f1, HKSFeature& f2) { return f1.m_scale > f2.m_scale;});
	for (auto iter = vF.begin(); iter != vF.end(); ++iter)
		cout << "(" << iter->m_index << "," << iter->m_scale << ")" << ", ";
	cout << endl;
	/* End LOG-FEATURES */

	double tl = DEFAULT_MATCH_TIME_LOW;	//10.0
	for (auto iter = vF.begin(); iter != vF.end(); )
	{
		const int vi = iter->m_index;
		double tu = findTmax(fineMesh, vi);	//?? max time to boundary? why geo*geo/4?
		int tn;
		if (tu < 0)	// no boundary
		{
			tn = 8;
			tu = std::pow(2.0, tn-1) * tl;				// no boundary
		}
		else
		{
			if(tu < tl) { vF.erase(iter); continue;}	// features too close to boundary are discarded
			tn = (int)ceil(log(tu/tl)/log(2.0));	    // tn is the number of timescales
		}
		iter->setTimes(tl, tu, tn);

		++iter;
	}

	MeshFeatureList *mfl = new MeshFeatureList;

	for (auto iter = vF.begin(); iter != vF.end(); ++iter)
	{
		mfl->addFeature(new HKSFeature(*iter));
	}

	mfl->setIDandName(FEATURE_MULTI_HKS, "Feature_multiple_hks");
	pMP->replaceProperty(mfl);
	pMP->setActiveFeaturesByID(FEATURE_MULTI_HKS);

	m_bFeatureDetected = true;
}

void DiffusionShapeMatcher::matchFeatures( std::ostream& flog, double matchThresh /*= DEFAULT_MATCH_THRESH*/ )
{
	const CMesh *mesh1 = pOriginalMesh[0], *mesh2 = pOriginalMesh[1];
	const vector<HKSFeature>& vftFine1 = vFeatures[0];
	const vector<HKSFeature>& vftFine2 = vFeatures[1];
	
	vector<HKSFeature> vftCoarse1, vftCoarse2;
	std::vector<MatchPair> matchedPairsCoarse, matchedPairsFine;

	flog << "Before defining vftCoarse" << endl;

	for_each(vftFine1.begin(), vftFine1.end(), [&](const HKSFeature& f){ 
		if(f.m_scale >= 2) 
			vftCoarse1.push_back(f); 
	});
	for_each(vftFine2.begin(), vftFine2.end(), [&](const HKSFeature& f){ 
		if(f.m_scale >= 2)
			vftCoarse2.push_back(f);
	});

	std::sort(vftCoarse1.begin(), vftCoarse1.end(), [](const HKSFeature& f1, const HKSFeature& f2) { return f1.m_index < f2.m_index; });
	std::sort(vftCoarse2.begin(), vftCoarse2.end(), [](const HKSFeature& f1, const HKSFeature& f2) { return f1.m_index < f2.m_index; });

	vector<MatchPair> vTmpMatchPairs;
	vector<double> vFeatureMatchScores;
	
////////////////    1. Match coarse features    ////////////////
//////////////////////////////////////////////////////////////////////////
	int size1 = (int) vftCoarse1.size();
	int size2 = (int) vftCoarse2.size();
	vector<VectorND> vsig1(size1), vsig2(size2);
	
	flog << "-- vftCoarse1: " << size1 << "; vftCoarse2: " << size2 << "\nBefore calVertexSignature" << endl;
	flog << "-- Coarse features 1: ";
	for (int i = 0; i < size1; ++i)
		flog << vftCoarse1[i].m_index << ' ';
	flog << "\n-- Coarse features 2: ";
	for (int i = 0; i < size2; ++i)
		flog << vftCoarse2[i].m_index << ' ';
	flog << endl;

// 	vector<HKSFeature>* vftCoarse[2] = {&vftCoarse1, &vftCoarse2};
// 	vector<VectorND>* vsig[2] = {&vsig1, &vsig2};
// 
// 	Concurrency::parallel_for(0, 1, [](int obj) {
// 		int size = vftCoarse[obj]->size();
// 		for(int i = 0; i < size; i++)
// 			calVertexSignature(pOriginalProcessor[obj], vftCoarse[obj]->at(i), vsig[obj]->at(i));
// 	});

#if (0)//(_MSC_VER >= 1600)
 	Concurrency::parallel_for (0, size1, [&](int i1) {
 		calVertexSignature(pOriginalProcessor[0], vftCoarse1[i1], vsig1[i1]);
 	});
  	Concurrency::parallel_for(0, size2, [&](int i2) {
 		calVertexSignature(pOriginalProcessor[1], vftCoarse2[i2], vsig2[i2]);
 	});
#else
  	for(int i1 = 0; i1 < size1; i1++)
  		calVertexSignature(pOriginalProcessor[0], vftCoarse1[i1], vsig1[i1]);
    for(int i2 = 0; i2 < size2; i2++)
   		calVertexSignature(pOriginalProcessor[1], vftCoarse2[i2], vsig2[i2]);
#endif

	flog << "-- candidates: ";
	double sigma1 = 4.0 * matchThresh;
	for(int i = 0; i < size1; i++)
	{
		for(int j = 0; j < size2; j++)
		{
			if (vftCoarse1[i].m_scale != vftCoarse2[j].m_scale) continue;
			double tl = 0.0;
			int tn = 0;
			double d = distFeature(vftCoarse1[i], vftCoarse2[j], vsig1[i], vsig2[j], tl, tn);	//average on each overlapped scale
			if(d < matchThresh)
			{
				flog << '(' << vftCoarse1[i].m_index <<','<< vftCoarse2[j].m_index << ") "; 

				double score = std::exp(-d/sigma1);
				vTmpMatchPairs.push_back(MatchPair(vftCoarse1[i].m_index, vftCoarse2[j].m_index, tl, tn, score));
				vFeatureMatchScores.push_back(score);
			}
		}
	}
	flog << endl;

	/* ---- create affinity matrix (compatibility of each candidate match) ---- */
	const int affinitySize = (int)vTmpMatchPairs.size();
    flog << "\nAffinity Matrix size: " << affinitySize << endl;

	mxArray *AM, *VM, *VA;
	AM = mxCreateDoubleMatrix(affinitySize, affinitySize, mxREAL);
	double *am = mxGetPr(AM);
	double sigma2 = 0.02;//0.02; 0.1

	for(int i = 0; i < affinitySize; i++)
	{
		am[i*affinitySize+i] = vFeatureMatchScores[i]; // diagonal

		for(int j = i+1; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i].m_idx1 == vTmpMatchPairs[j].m_idx1)^(vTmpMatchPairs[i].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor conflict; 1-to-multiple not desired? 
			{
				am[i*affinitySize+j] = am[j*affinitySize+i] = 0.0;
				continue;
			}
			double ds = 0.0;
			ds = distFeaturePair(pOriginalProcessor[0], pOriginalProcessor[1], vTmpMatchPairs[i], vTmpMatchPairs[j]);
			am[i*affinitySize+j] = am[j*affinitySize+i] = exp(-ds/sigma2);
		}
	}
	flog << "Affinity matrix constructed!" << endl;

	/* ---- solving the greatest eigenvector (PCA?) ---- */
	engPutVariable(m_ep, "AM", AM);
	engEvalString(m_ep, "[VM,VA] = spectral_embedding(AM);");	//computing leading eigenvector of A using svd
	VM = engGetVariable(m_ep, "VM");
	double *vm = mxGetPr(VM);
	VA = engGetVariable(m_ep, "VA");
	double *va = mxGetPr(VA);		//?? not referenced?

	const double c_thresh = 0.01;	//was 0.1
	std::vector<MatchPair> mpc1;

	while(1)
	{
		int i_max = -1;
		double v_max = 0;
		for(int i = 0; i < affinitySize; i++)
		{
			if(abs(vm[i]) > v_max)
			{
				v_max = abs(vm[i]);
				i_max = i;
			}
		}

		if (v_max <= c_thresh) 
		{
			flog << "-- Discarded v_max: " << v_max << endl;
			break;
		}

		flog << "  " << mpc1.size() << ": " << i_max << ',' << v_max << endl;

		bool hit = true;
		const int curMatchSize = (int)mpc1.size();
//*
		if (curMatchSize >= 4)
		{
			int refIndex[4];
			for (int i = 0; i < 4; ++i)
				refIndex[i] =  int( double(curMatchSize)/4.0 * (i + (double)rand()/(double)RAND_MAX)); 
			
			double errorSum(0);
			for (int i = 0; i < 4; ++i)
			{
				int selectIdx = refIndex[i];
				const MatchPair &mp1 = mpc1[selectIdx], &mp2 = vTmpMatchPairs[i_max];
				double geodist1 = mesh1->getGeodesic(mp1.m_idx1, mp2.m_idx1),
					   geodist2 = mesh2->getGeodesic(mp1.m_idx2, mp2.m_idx2);
				double distError = std::abs(geodist1 - geodist2);

				if ( distError >= (geodist1 < 10*mesh1->getAvgEdgeLength() ? 3*mesh1->getAvgEdgeLength() : 0.3*geodist1) )
				{
					hit = false;
					break;
				}
				errorSum += distError;
			}
			if (hit)
			{
				vector<ExtremaPoint> mpCandidates;
				for (int j = 0; j < affinitySize; ++j)
				{
					if (j != i_max && 
						vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1 && 
						vTmpMatchPairs[i_max].m_idx2 != vTmpMatchPairs[j].m_idx2 &&
						abs(vm[j]) > c_thresh
						)
					{
						mpCandidates.push_back(ExtremaPoint(j, abs(vm[j])));
					}
				}

				int candSize = (int)mpCandidates.size();

				for (int k = 0; k < candSize; ++k)
				{
					double errorSumK(0);
					bool pass = true;
					for (int i = 0; i < 4; ++i)
					{
						int selectIdx = refIndex[i];

						const MatchPair &mp1 = mpc1[selectIdx], &mp2 = vTmpMatchPairs[mpCandidates[k].index];
						double geodist1 = mesh1->getGeodesic(mp1.m_idx1, mp2.m_idx1),
							   geodist2 = mesh2->getGeodesic(mp1.m_idx2, mp2.m_idx2);
						double distError = std::abs(geodist1 - geodist2);

						if ( distError >= (geodist1 < 10*mesh1->getAvgEdgeLength() ? 3*mesh1->getAvgEdgeLength() : 0.3*geodist1) )
						{
							pass = false;
							break;
						}
						errorSumK += distError;
					}
					if (pass && errorSumK < errorSum)
					{
						i_max = mpCandidates[k].index;
						errorSum = errorSumK;
					}
				}
			}
		}
//*/
		vm[i_max] = 0.0;
		if (!hit) continue;
		// hit
		mpc1.push_back(vTmpMatchPairs[i_max]);
		
		// now a max candidate is found, remove conflicting candidates
		for(int j = 0; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1) ^ (vTmpMatchPairs[i_max].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor, conflict
				vm[j] = 0.0;
		}
	} // end of while()

	matchedPairsCoarse = mpc1;

	flog << "Coarse Match computed!" << endl;

	mxDestroyArray(VM);
	mxDestroyArray(VA);
	mxDestroyArray(AM);

	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); /*++iter*/ )
	{
		if (iter->m_idx1 == iter->m_idx2)
		{
			iter++;
			continue;
		}

		if (!mesh2->isInNeighborRing(iter->m_idx1, iter->m_idx2, 5))	
		{
			iter->m_note = -1;
			iter = matchedPairsCoarse.erase(iter);
			//			iter++;
		}
		else iter++;
	}

//////////////////////////////////////////////////////////////////////////
////////    2. Match rest features if possible    ////
	flog << "------------------------------------------" << endl;
	flog << "\nMatch rest features if possible" << endl;
	flog << "  Mesh 1 features: ";
	for (vector<HKSFeature>::const_iterator hiter = vftFine1.begin(); hiter != vftFine1.end(); ++hiter)
		flog << hiter->m_index  << ',' << hiter->m_scale << "  ";
	flog << endl;
	flog << "  Mesh 2 features: ";
	for (vector<HKSFeature>::const_iterator hiter = vftFine2.begin(); hiter != vftFine2.end(); ++hiter)
		flog << hiter->m_index  << ',' << hiter->m_scale << "  ";
	flog << endl;	
	
	matchedPairsFine = matchedPairsCoarse;

	const int coarseMatchSize = (int) matchedPairsCoarse.size();
	vector<int> coarseFeat1, coarseFeat2;
	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); ++iter)
	{
		coarseFeat1.push_back(iter->m_idx1);
		coarseFeat2.push_back(iter->m_idx2);
	}

	flog << "!!test ok 1" << endl;

	vector<int> restFeat1, restFeat2;
	map<int, vector<HKSFeature>::const_iterator> restIdxToHKS1, restIdxToHKS2;
	for (vector<HKSFeature>::const_iterator iter = vftFine1.begin(); iter != vftFine1.end(); ++iter)
	{
		if (find(coarseFeat1.begin(), coarseFeat1.end(), iter->m_index) == coarseFeat1.end())
		{
			restFeat1.push_back(iter->m_index);
			restIdxToHKS1[iter->m_index] = iter;
		}
	}
	flog << "!!test ok 2" << endl;

	for (vector<HKSFeature>::const_iterator iter = vftFine2.begin(); iter != vftFine2.end(); ++iter)
	{
		if (find(coarseFeat2.begin(), coarseFeat2.end(), iter->m_index) == coarseFeat2.end())
		{
			restFeat2.push_back(iter->m_index);
			restIdxToHKS2[iter->m_index] = iter;
		}
	}

	const int restSize1 = (int) restFeat1.size(), restSize2 = (int) restFeat2.size();
	flog << "Rest Size1: " << restSize1 << '\t' << "Rest Size2: " << restSize2 << endl;

	vector<VectorND> vHKC1, vHKC2;
	vHKC1.resize(restSize1);
	vHKC2.resize(restSize2);
	for (int j = 0; j < restSize1; ++j)
	{
		VectorND& vnd = vHKC1[j];
		vnd.reserve(coarseMatchSize);
		for (int i = 0; i < coarseMatchSize; ++i)
		{
			vnd.m_vec[i] = pOriginalProcessor[0]->calHK(restFeat1[j], matchedPairsCoarse[i].m_idx1, 80);
//			vnd.m_vec[i] = mesh1->CalGeodesic(restFeat1[j], matchedPairsCoarse[i].m_idx1);
		}
	}
	flog << "!!test ok 3" << endl;

	for (int j = 0; j < restSize2; ++j)
	{
		VectorND& vnd = vHKC2[j];
		vnd.reserve(coarseMatchSize);
		for (int i = 0; i < coarseMatchSize; ++i)
		{
			vnd.m_vec[i] = pOriginalProcessor[1]->calHK(restFeat2[j], matchedPairsCoarse[i].m_idx2, 80);
//			vnd.m_vec[i] = mesh2->CalGeodesic(restFeat2[j], matchedPairsCoarse[i].m_idx2);
		}
	}
	flog << "!!test ok 4" << endl;

	vector<MatchPair> restMatches;
	for (int i = 0; i < restSize1; ++i)
	{
		for (int j = 0; j < restSize2; ++j)
		{
			if (vftFine2[restFeat1[i]].m_scale != vftFine2[restFeat2[j]].m_scale)
				continue;
			double s = std::exp(-vHKC1[i].calDistance(vHKC2[j]));
			restMatches.push_back(MatchPair(restFeat1[i], restFeat2[j], s));
		}
	}

	flog << "Candidate rest matches: " << restMatches.size() << endl;
	while (!restMatches.empty())
	{
		vector<MatchPair>::iterator maxIter = restMatches.end();
		double maxVal = -1;
		for (vector<MatchPair>::iterator iter = restMatches.begin(); iter != restMatches.end(); ++iter)
		{
			if (iter->m_score > maxVal)
			{
				maxVal = iter->m_score;
				maxIter = iter;
			}
		}
		int maxIdx1 = maxIter->m_idx1, maxIdx2 = maxIter->m_idx2;

		bool pass = true;
		int matchSize = (int) matchedPairsFine.size();
		if (matchSize >= 5)	// from already registered features to judge the validity of new matches
		{
			for (int k = 0; k < 5; ++k)
			{
				int selectIdx = int( double(matchSize)/5.0 * (k + double(rand()/double(RAND_MAX)))); 
				//				int selectIdx = rand() % matchSize;
				const MatchPair &mp1 = matchedPairsFine[selectIdx];
				double geodist1 = mesh1->getGeodesic(mp1.m_idx1, maxIdx1),
					   geodist2 = mesh2->getGeodesic(mp1.m_idx2, maxIdx2);
				double distError = std::abs(geodist1 - geodist2);

				if ( distError >= (geodist1 < 10*mesh1->getAvgEdgeLength() ? 2*mesh1->getAvgEdgeLength() : 0.2*geodist1) ) 
				{
					pass = false;
					break;
				}
			}
		}
		if (pass) 
		{
			matchedPairsFine.push_back(MatchPair(maxIdx1, maxIdx2));
			flog << "  add new match: " << maxIdx1 << ", " << maxIdx2 << endl;
		}
		else 
		{
			restMatches.erase(maxIter);
			continue;
		}

		for (vector<MatchPair>::iterator iter = restMatches.begin(); iter != restMatches.end(); )
		{
			if (iter->m_idx1 == maxIdx1 || iter->m_idx2 == maxIdx2)
				restMatches.erase(iter);
			else iter++;
		}
	}

	flog << "  Fine matches computed" << endl;
	flog << " --- Total Matched Features: " << matchedPairsFine.size() << endl;
	flog << " --- Mesh 1 features num: " << vftFine1.size() << endl;
	flog << " --- Mesh 2 features num: " << vftFine2.size() << endl;

	flog << "\nResults: " << endl;
	for(auto iter = matchedPairsFine.begin(); iter != matchedPairsFine.end(); ++iter)
		flog << "(" << iter->m_idx1 << ", " << iter->m_idx2 << ")" << endl;

	//* add note to matchings with great discrepancy *//	
	for (vector<MatchPair>::iterator iter = matchedPairsFine.begin(); iter != matchedPairsFine.end(); /*++iter*/ )
	{
		for (vector<HKSFeature>::const_iterator hiter = vftFine1.begin(); hiter != vftFine1.end(); ++hiter)
		{
			if (hiter->m_scale == 3 && hiter->m_index == iter->m_idx1)
				flog << " --- Scale-3 Feature Pair: " << iter->m_idx1 << ", " << iter->m_idx2 << endl;
		}

		if (iter->m_idx1 == iter->m_idx2)
		{
			iter++;
			continue;
		}

		if (!mesh2->isInNeighborRing(iter->m_idx1, iter->m_idx2, 5))	
		{
			iter->m_note = -1;
			iter = matchedPairsFine.erase(iter);
		}
		else iter++;
	}

	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); ++iter )
	{
		if (iter->m_idx1 == iter->m_idx2) 
			continue;

		if (!mesh2->isInNeighborRing(iter->m_idx1, iter->m_idx2, 5))	
			iter->m_note = -1;
		//iter = matchedPairs.erase(iter);
	}

	m_bFeatureMatched = true;
	m_nAlreadyMatchedLevel = m_nRegistrationLevels;
	vFeatureMatchingResults[m_nAlreadyMatchedLevel] = matchedPairsFine;
} // DiffusionShapmeMatcher::matchFeatures()

void DiffusionShapeMatcher::calVertexSignature( const DifferentialMeshProcessor* pOriginalProcessor, const HKSFeature& hf, VectorND& sig)
{
	double t = hf.m_tl;
	sig.reserve(hf.m_tn);
	for(int i = 0; i < hf.m_tn; i++)
	{
		double eref = 4.0*PI*t;
		sig.m_vec[i] = std::log(pOriginalProcessor->calHK(hf.m_index, hf.m_index, t) * eref);	//normalization to balance HKS at different t
		t *= 2.0;
	}
}

void DiffusionShapeMatcher::refineRegister( std::ostream& flog )
{
	if (m_nAlreadyRegisteredLevel == 0) return;

	const int current_level = max(m_nAlreadyRegisteredLevel - 1, 0);	//the registration level we currently working on
	
	flog << "\n-------- registration level " << current_level << " --------" << endl;

	double anchorThresh;
	double regT = 20.0 * pow(2.0, current_level);
	if (current_level == 1) 
	{
		regT = 80;
		anchorThresh = 1e-4;
	}
	if (current_level == 0)
	{
		regT = 20;
		anchorThresh = 2e-4;
	}
	anchorThresh = 1e-4;
//	thresh = (double)featureMatch.size() * 1e-5;

	CMesh *tmesh1 = meshPyramids[0].getMesh(current_level), *tmesh2 = meshPyramids[1].getMesh(current_level);
	const int coarseSize1 = tmesh1->getVerticesNum(), coarseSize2 = tmesh2->getVerticesNum();
	vector<int> vMatch1(coarseSize1, -1), vMatch2(coarseSize2, -1);
	vector<double> vMatchScore1(coarseSize1, 0), vMatchScore2(coarseSize2, 0);

	const vector<MatchPair> oldAnchorMatch = getMatchedFeaturesResults(m_nAlreadyMatchedLevel);
	const vector<MatchPair> oldReg = (current_level == m_nRegistrationLevels - 1) ?
									 oldAnchorMatch : getRegistrationResults(m_nAlreadyRegisteredLevel);

	vector<MatchPair> newAnchorMatch;
	vector<MatchPair> tmpReg1, tmpReg2;

	/************************************************************************/
	/*                        1. Initialize/Compute HKC                     */
	/************************************************************************/
	vector<int> vFeatureID1, vFeatureID2;
	for (auto iter = begin(oldAnchorMatch); iter != end(oldAnchorMatch); ++iter)
	{
		vFeatureID1.push_back(iter->m_idx1);
		vFeatureID2.push_back(iter->m_idx2);
	}
	
	flog << "Anchors IDs: " << endl;
	for (auto iter = begin(oldAnchorMatch); iter != end(oldAnchorMatch); ++iter)
		flog << "(" << iter->m_idx1 << ',' << iter->m_idx2 << ") ";
	flog << endl << endl;

	CStopWatch timer;
	timer.startTimer();
	m_ParamMgr[0].para_computeHKC(vFeatureID1, regT); 
	m_ParamMgr[1].para_computeHKC(vFeatureID2, regT); 
	timer.stopTimer();
	cout << "HKParam computation time: " << timer.getElapsedTime() << 's' << endl;

	std::vector<PointParam>& vCoordinates1 = m_ParamMgr[0].vCoord;
	std::vector<PointParam>& vCoordinates2 = m_ParamMgr[1].vCoord;

	flog << "Multi-anchor coordinates initialized!" << endl;
	
	set<MatchPair> uniquePairs;
	for (vector<MatchPair>::const_iterator mpIter = oldReg.begin(); mpIter != oldReg.end(); ++mpIter)
	{
		const int vid_i = mpIter->m_idx1, vid_j = mpIter->m_idx2;
		const int vi = id2Index(0, vid_i, current_level), vj = id2Index(1, vid_j, current_level);

		// filter out repeated matches
		if (uniquePairs.find(MatchPair(vi, vj)) == uniquePairs.end())
		{
			uniquePairs.insert(MatchPair(vi, vj));
			tmpReg1.push_back(*mpIter);			
		}		
	}

	for (auto iter = tmpReg1.begin(); iter != tmpReg1.end(); ++iter)
	{
		const int vid_i = iter->m_idx1, vid_j = iter->m_idx2;
		const int vi = id2Index(0, vid_i, current_level), vj = id2Index(1, vid_j, current_level);

		if (find(oldAnchorMatch.begin(), oldAnchorMatch.end(), *iter) != oldAnchorMatch.end())
			iter->m_note = 1;

		if (iter->m_note == 1) 
			iter->m_score = 1.;	// matched feature pairs have the highest score 1.0
		else
		{
			double score = computeMatchScore(iter->m_idx1, iter->m_idx2);
			iter->m_score = score;
		}		

		vMatch1[vi] = vj;
		vMatch2[vj] = vi;
		vMatchScore1[vi] = iter->m_score;
		vMatchScore2[vj] = iter->m_score;
	}

	flog << "\n1. HKC computed!\n" 
		<< "  Old anchor size: " << oldAnchorMatch.size() 
		<< "\n  Old reg size: " << oldReg.size() 
		<< "\n  Filtered old reg size: " << tmpReg1.size() << endl;

	/************************************************************************/
	/*              2. Registration                                         */
	/************************************************************************/
	NoteQueue qscan1;
	for (auto iter = tmpReg1.begin(); iter != tmpReg1.end(); ++iter)
	{
		int vid_1 = iter->m_idx1, vid_2 = iter->m_idx2;
		int v1 = id2Index(0, vid_1, current_level), v2 = id2Index(1, vid_2, current_level);
		//MyNote mn(v1, v2, iter->m_score);	//id, score
		//MyNote mn(v1, tmesh1->m_pVertex[v1].m_vParam.m_votes);
		MyNote mn(v1, v2, vCoordinates1[vid_1].m_votes/* * iter->m_score*/);
		qscan1.push(mn);
	}

	while (!qscan1.empty())
	{
		MyNote qt = qscan1.top();
		qscan1.pop();

		int vi = qt.m_idx1, vj = qt.m_idx2;		//index on this level

		if (vMatch2[vj] != vi) continue;

		int vid_i = tmesh1->getVertex_const(vi)->getVID(),
			vid_j = tmesh2->getVertex_const(vj)->getVID();

		for (int ei = 0; ei < tmesh1->getVertex_const(vi)->getValence(); ei++)
		{
			const CHalfEdge* he = tmesh1->getVertex_const(vi)->getHalfEdge_const(ei);
			const int vt = he->getVertexIndex(1);
			const int vid_t = tmesh1->getVertex_const(vt)->getVID();

			if( vMatch1[vt] >= 0 && vMatch1[vt] < coarseSize2 ) 
				continue;  // already registered

			if (tmesh1->getVertex_const(vt)->isOnBoundary())
				continue;
			//  if( matcher1.m_vHKParamFine[vt].m_votes < thresh )
			//		continue;	// low priority points skipped

			double score;
			int vm = searchVertexMatch(vt, vj, current_level, 2, score);
			if(vm >= 0) 
			{
				int vid_m = tmesh2->getVertex_const(vm)->getVID();
				vMatch1[vt] = vm;
				vMatchScore1[vt] = score;
				tmpReg1.push_back(MatchPair(vid_t, vid_m, score));

				if (vMatch2[vm] == -1 || vMatchScore2[vm] < score)
				{
					vMatch2[vm] = vt;
					vMatchScore2[vm] = score;
				}

				//MyNote mn(vt, vm, score);
				MyNote mn(vt, vm, vCoordinates1[vid_t].m_votes);
				qscan1.push(mn);
			}
			else flog << "No match 3!" << endl;
		}
	} //while(!qscan1.empty())

	flog << "\n2. Crude registration finished!" << "\n  Tmp reg size: " << tmpReg1.size() << endl;

	/************************************************************************/
	/*              3. Adjust and update anchor points                      */
	/************************************************************************/

	newAnchorMatch = oldAnchorMatch;
	if(current_level > 0)
	{
		for (auto iter = begin(tmpReg1); iter != end(tmpReg1); ++iter)
		{
			int vid_1 = iter->m_idx1, vid_2 = iter->m_idx2;
			int v1 = id2Index(0, vid_1, current_level), v2 = id2Index(1, vid_2, current_level);
			if (vMatch1[v1] == v2 && vMatch2[v2] == v1)		// select only the pair considered best mutually
				tmpReg2.push_back(*iter);
		}

		priority_queue<MatchPair, vector<MatchPair>, greater<MatchPair> > qReg;
		for (auto iter = begin(tmpReg2); iter != end(tmpReg2); ++iter)
		{
			const int vi = iter->m_idx1, vj = iter->m_idx2;
			double score;
			const int vm1 = searchVertexMatch(vi, vj, 0, 2, score, current_level);	// vi, vj here are vertex id as well as level-0 index
			qReg.push(MatchPair(vi, vm1, vCoordinates1[vi].m_votes * score));
		}

		tmpReg1.clear();
		while(!qReg.empty())
		{
			tmpReg1.push_back(qReg.top());
			qReg.pop();
		}

		// insert additional anchor points
		for (auto riter = tmpReg1.begin(); riter != tmpReg1.end(); ++riter)
		{
			if (vCoordinates1[riter->m_idx1].m_votes < 5*anchorThresh ||
				vCoordinates2[riter->m_idx2].m_votes < 5*anchorThresh
				) 
				continue;

			bool pass = true;
			for (auto fiter = newAnchorMatch.begin(); fiter != newAnchorMatch.end(); ++fiter)
			{
				if (fiter->m_idx1 == riter->m_idx1 || fiter->m_idx2 == riter->m_idx2 ||
					pOriginalMesh[0]->isInNeighborRing(fiter->m_idx1, riter->m_idx1, 6) ||
					pOriginalMesh[1]->isInNeighborRing(fiter->m_idx2, riter->m_idx2, 6)
					) 
				{ 
					pass = false; 
					break; 
				}
			}
			if (pass)
				newAnchorMatch.push_back(*riter);
		}
	}

	const vector<MatchPair>& newReg = tmpReg1;

	flog << "\n3. Registration adjustment and anchor set augmentation finished!\n"
		<< "  New reg size: " << newReg.size() << "; new anchor size: " << newAnchorMatch.size() << endl;

	if (m_nAlreadyMatchedLevel > 0) m_nAlreadyMatchedLevel--;
	if (m_nAlreadyRegisteredLevel > 0) m_nAlreadyRegisteredLevel--;
	vFeatureMatchingResults[m_nAlreadyMatchedLevel] = newAnchorMatch;
	vRegistrationResutls[m_nAlreadyRegisteredLevel] = newReg;

} // refineRegister

void DiffusionShapeMatcher::setRegistrationLevels( int val )
{
	m_nRegistrationLevels = min(val, m_nPyramidLevels);

	vFeatureMatchingResults.resize(m_nRegistrationLevels + 1);
	vRegistrationResutls.resize(m_nRegistrationLevels);

	m_nAlreadyRegisteredLevel = m_nRegistrationLevels;
	m_nAlreadyMatchedLevel = m_nRegistrationLevels+1;
}

const std::vector<MatchPair>& DiffusionShapeMatcher::getMatchedFeaturesResults ( int level ) const
{
 	if (level < -1 || level > m_nRegistrationLevels) 
		throw logic_error("Selected level out of bound!");
 	else if (level == -1) return vFeatureMatchingResults[m_nRegistrationLevels];
 	else return vFeatureMatchingResults[level];	
}

const std::vector<MatchPair>& DiffusionShapeMatcher::getRegistrationResults( int level ) const
{
	if (level < -1 || level >= m_nRegistrationLevels) 
		throw logic_error("Selected level out of bound!");
	else if (level == -1) return vRegistrationResutls[m_nRegistrationLevels-1];
	else return vRegistrationResutls[level];	
}

double DiffusionShapeMatcher::computeMatchScore( int idx1, int idx2 ) const
{
	const PointParam &hkp1 = m_ParamMgr[0].vCoord[idx1],
		             &hkp2 = m_ParamMgr[1].vCoord[idx2];

	return std::exp(-hkp1.calDistance(hkp2, 1));
}

int DiffusionShapeMatcher::searchVertexMatch( const int vt, const int vj, const int level, const int ring, double& match_score, int upper_level )
{
	assert( level < m_nRegistrationLevels);

	// search vi's match in vj's neighborhood
	const CMesh* tmesh1 = getMesh(0, level);
	const CMesh* tmesh2 = getMesh(1, level);
	int vid_i = tmesh1->getVertex_const(vt)->getVID();
	int vid_j = tmesh2->getVertex_const(vj)->getVID();

	/* find vj's neighbors/covered vertices */
	list<int> vNeighbor;
	set<int> marked_set;
	marked_set.insert(vj);
	vNeighbor.push_back(vj);

	if (upper_level > level)
	{
		vector<int> vCover;
		vCover.push_back(vid_j);

		for (int l = upper_level-1; l >= level; --l)
		{
			vector<int> vCoverTmp = vCover;
			for (vector<int>::iterator iter = vCover.begin(); iter != vCover.end(); ++iter)
			{
				int idxL = id2Index(1, *iter, l);
				list<int> idxCovered = getMeshPyramid(1).getCoveredVertexList(l, idxL);
				for (list<int>::iterator citer = idxCovered.begin(); citer != idxCovered.end(); ++citer)
				{
					int newId = getMesh(1, l)->getVertex_const(*citer)->getVID();
					vCoverTmp.push_back(newId);
				}
			}
			vCover = vCoverTmp;
		}
		for (vector<int>::iterator iter = vCover.begin(); iter != vCover.end(); ++iter)
		{
			int idxL = id2Index(1, *iter, level);
			if (marked_set.find(idxL) != marked_set.end())
			{
				vNeighbor.push_back(idxL);
				marked_set.insert(idxL);
			}
		}
	}

	list<int> nb1, nb2;
	nb1 = vNeighbor;
	for (int r = 0; r < ring; r++)
	{
		nb2.clear();
		for (list<int>::iterator iter = nb1.begin(); iter != nb1.end(); ++iter)
		{
			int idx = *iter;
			for (int l = 0; l < tmesh2->getVertex_const(idx)->getValence(); ++l)
			{
				const CHalfEdge* he = tmesh2->getVertex_const(idx)->getHalfEdge_const(l);
				int vt = he->getVertexIndex(1);
				if (marked_set.find(vt) == marked_set.end())
				{
					marked_set.insert(vt);
					vNeighbor.push_back(vt);
					nb2.push_back(vt);
				}				
			}
		}
		nb1 = nb2;
	}
	// now vNeighbor contains the index set of indices covered by vj in current level

	/* ---- find the maximum match ---- */
	int vmatch = -1;
	double smax = -DBL_MAX;
	for (list<int>::iterator iter = vNeighbor.begin(); iter != vNeighbor.end(); ++iter)
	{
		int vt = *iter;
// 		if (tmesh2->m_pVertex[vt].m_vMatched >= 0 && tmesh2->m_pVertex[vt].m_vMatched < tmesh1->m_nVertex)
// 			continue;	//injection, not a good idea
		int vid_t = tmesh2->getVertex_const(vt)->getVID();
		if (tmesh2->getVertex_const(vt)->isOnBoundary())
			continue;
		double dt = computeMatchScore(vid_i, vid_t);
		if(dt > smax)
		{
			vmatch = vt;
			smax = dt;
		}
	}

	match_score = smax;
	return vmatch;
}

void DiffusionShapeMatcher::ComputeTensorFeature( const DifferentialMeshProcessor* pmp, int i, int j, int k, double t, double* sang)
{

	if(i==j || i==k || j==k)		//not a triangle
	{
		sang[0] = -10.0;
		sang[1] = -10.0;
		sang[2] = -10.0;
		return;
	}

	double d1 = pmp->calHK(i,j,t);
	double d2 = pmp->calHK(j,k,t);
	double d3 = pmp->calHK(k,i,t);

	if(d1<0.0) d1 = 1e-15;
	if(d2<0.0) d2 = 1e-15;
	if(d3<0.0) d3 = 1e-15;

	d1 = sqrt(-4.0*t*log(d1));
	d2 = sqrt(-4.0*t*log(d2));
	d3 = sqrt(-4.0*t*log(d3));

	// cot

	//sang[0] = (d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
	//sang[1] = (d2*d2+d3*d3-d1*d1)/(2.0*d2*d3);
	//sang[2] = (d3*d3+d1*d1-d2*d2)/(2.0*d3*d1);

	// sine theorem: area = 2*a*b*sin(C)
	double s = (d1+d2+d3)/2.0;
	double a = sqrt(s*(s-d1)*(s-d2)*(s-d3));

	sang[0] = a/(2*d1*d2);
	sang[1] = a/(2*d2*d3);
	sang[2] = a/(2*d3*d1);
}


double DiffusionShapeMatcher::TensorMatching(Engine *ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, Cluster& ct1, Cluster& ct2, vector<MatchPair>& matched, double t, double thresh)
{
	// generate triangles
	int vsize1 = (int) ct1.m_member.size();	// input feature size 1
	int vsize2 = (int) ct2.m_member.size(); // input feature size 2

	vector<int> triangs;
	int i,j,k,tsize1,tsize2;

	// ***************************************
	// improve to local triangles
	// ***************************************/
	
	if(vsize1 > 8) 
	{
		for(i=0; i<vsize1; i++)
		{
			for(j=0; j<8; j++)
			{
				if(i==j) continue;
				for(k=vsize1-1; k>vsize1-9; k--)
				{
					if(i==k || j==k) continue;
					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}
			}
		}
	}
	else
	{
		for(i=0; i<vsize1; i++)
		{
			for(j=0; j<vsize1; j++)
			{
				if(i==j) continue;
				for(k=0; k<vsize1; k++)
				{
					if(i==k || j==k) continue;
					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}
			}
		}
	}

	tsize1 = (int)triangs.size() / 3;
	tsize2 = vsize2 * vsize2 * vsize2;

	// compute feature descriptors

	mxArray *feat1, *feat2, *tris, *numbs, *mX2, *vX2, *score;
	feat1 = mxCreateDoubleMatrix(3, tsize1, mxREAL);
	double *pfeat1 = mxGetPr(feat1);
	feat2 = mxCreateDoubleMatrix(3, tsize2, mxREAL);
	double *pfeat2 = mxGetPr(feat2);
	tris = mxCreateDoubleMatrix(3, tsize1, mxREAL);
	double *ptris = mxGetPr(tris);
	numbs = mxCreateDoubleMatrix(1, 4, mxREAL);
	double *pnumbs = mxGetPr(numbs);

	pnumbs[0] = vsize1; // nP1
	pnumbs[1] = vsize2; // nP2
	pnumbs[2] = tsize1; // nT
	if(tsize1 > 60) pnumbs[3] = 60;     // nNN
	else pnumbs[3] = tsize1*0.5;

	for(i=0; i<tsize1*3; i++)
		ptris[i] = triangs[i];

	for(i=0; i<tsize1; i++)
	{
		int vi = triangs[i*3];
		int vj = triangs[i*3+1];
		int vk = triangs[i*3+2];
		ComputeTensorFeature(pmp1, ct1.m_member[vi], ct1.m_member[vj], ct1.m_member[vk], t, &pfeat1[i*3]);
	}

	for(i=0; i<vsize2; i++)
	{
		for(j=0; j<vsize2; j++)
		{
			for(k=0; k<vsize2; k++)
			{
				ComputeTensorFeature(pmp2, ct2.m_member[i], ct2.m_member[j], ct2.m_member[k], t, &pfeat2[((i*vsize2+j)*vsize2+k)*3]);
			}
		}
	}

	// invoke matlab for tensor matching

	engPutVariable(ep,"feat1",feat1);
	engPutVariable(ep,"feat2",feat2);
	engPutVariable(ep,"tris",tris);
	engPutVariable(ep,"numbs",numbs);
//	system("PAUSE");	// for manually testing tensorMat
	engEvalString(ep, "[vX2,mX2,score]=tensorMat(feat1,feat2,tris,numbs);");
	vX2 = engGetVariable(ep, "vX2");
	double *pv = mxGetPr(vX2);	// best match of each shape-1 feature in the set of shape-2 features
	mX2 = engGetVariable(ep, "mX2");
	double *px = mxGetPr(mX2);  // corresponding match scores
	score = engGetVariable(ep, "score");
	double *ps = mxGetPr(score);

	// interpret results

	double result = ps[0];	// tensor matching score
	result = 0.0;

	int count = 0;
	while(count++ < vsize1)
	{
		double *pmax = max_element(pv, pv+vsize1);	// max feature match score
		int imax = pmax - pv;	// shape-1 feature index

		cout << imax << ": " << *pmax << endl;
		if(*pmax < thresh) 
			break;

		result += pv[imax];
		MatchPair mpt;
		mpt.m_idx1 = ct1.m_member[imax];
		int ind = (int)px[imax];		// shape-2 feature index matched to imax
		//if(ind<0 || ind>vsize2) continue;
		mpt.m_idx2 = ct2.m_member[ind-1];
		matched.push_back(mpt);

		pv[imax] = 0.0;
		// clear conflicted
		for(i=0; i<vsize1; i++)
		{
			if(px[i] == px[imax])
				pv[i] = 0.0;
		}
	}

	mxDestroyArray(feat1);
	mxDestroyArray(feat2);
	mxDestroyArray(tris);
	mxDestroyArray(numbs);
// 	mxDestroyArray(mX2);
// 	mxDestroyArray(vX2);
// 	mxDestroyArray(score);

	return result;
}

void DiffusionShapeMatcher::matchFeaturesTensor( std::ostream& flog, double timescale, double thresh )
{
	const CMesh *mesh1 = pOriginalMesh[0], *mesh2 = pOriginalMesh[1];
	const vector<HKSFeature>& vftFine1 = vFeatures[0];
	const vector<HKSFeature>& vftFine2 = vFeatures[1];
	
	vector<HKSFeature> vftCoarse1, vftCoarse2;
	std::vector<MatchPair> matchedPairsCoarse, matchedPairsFine;
	vector<MatchPair> vTmpMatchPairs;
	vector<double> vFeatureMatchScores;

	for_each(vftFine1.begin(), vftFine1.end(), [&](const HKSFeature& f){ 
//		if(f.m_scale >= 2)
			vftCoarse1.push_back(f); 
	});
	for_each(vftFine2.begin(), vftFine2.end(), [&](const HKSFeature& f){ 
//		if(f.m_scale >= 2)
			vftCoarse2.push_back(f);
	});

	std::sort(vftCoarse1.begin(), vftCoarse1.end(), [](const HKSFeature& f1, const HKSFeature& f2) { return f1.m_index < f2.m_index; });
	std::sort(vftCoarse2.begin(), vftCoarse2.end(), [](const HKSFeature& f1, const HKSFeature& f2) { return f1.m_index < f2.m_index; });
		
	int size1 = (int) vftCoarse1.size();
	int size2 = (int) vftCoarse2.size();
	
	flog << "vftCoarse1(" << size1 << "): ";
	for (auto iter = vftCoarse1.begin(); iter != vftCoarse1.end(); ++iter)
		flog << iter->m_index << " ";

	flog << "\nvftCoarse2(" << size2 << "): ";
	for (auto iter = vftCoarse2.begin(); iter != vftCoarse2.end(); ++iter)
		flog << iter->m_index << " ";
	flog << endl;

	Cluster cluster1, cluster2;
	for (auto iter = vftCoarse1.begin(); iter != vftCoarse1.end(); ++iter)
		cluster1.m_member.push_back(iter->m_index);
	for (auto iter = vftCoarse2.begin(); iter != vftCoarse2.end(); ++iter)
		cluster2.m_member.push_back(iter->m_index);

	vector<MatchPair> vPairs;
	double matchScore = TensorMatching(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], cluster1, cluster2, vPairs, timescale, thresh);
	flog << "Match score: " << matchScore << endl;

	vFeatureMatchingResults[m_nRegistrationLevels] = vPairs;
	m_nAlreadyMatchedLevel = m_nRegistrationLevels;
	m_bFeatureMatched = true;
}

void DiffusionShapeMatcher::getVertexCover( int obj, int vidx, int level, int upper_level, int ring, std::vector<int>& vCoveredIdx ) const
{
	const CMesh* tmesh = getMesh(obj, level);
	int vid = tmesh->getVertex_const(vidx)->getVID();

	list<int> vNeighbor;
	set<int> marked_set;
	marked_set.insert(vidx);
	vNeighbor.push_back(vidx);

	if (upper_level > level)
	{
		vector<int> vCover;
		vCover.push_back(vid);

		for (int l = upper_level-1; l >= level; --l)
		{
			vector<int> vCoverTmp = vCover;
			for (vector<int>::iterator iter = vCover.begin(); iter != vCover.end(); ++iter)
			{
				int idxL = id2Index(obj, *iter, l);
				list<int> idxCovered = getMeshPyramid(obj).getCoveredVertexList(l, idxL);
				for (list<int>::iterator citer = idxCovered.begin(); citer != idxCovered.end(); ++citer)
				{
					int newId = getMesh(obj, l)->getVertex_const(*citer)->getVID();
					vCoverTmp.push_back(newId);
				}
			}
			vCover = vCoverTmp;
		}
		for (vector<int>::iterator iter = vCover.begin(); iter != vCover.end(); ++iter)
		{
			int idxL = id2Index(obj, *iter, level);
			if (marked_set.find(idxL) != marked_set.end())
			{
				vNeighbor.push_back(idxL);
				marked_set.insert(idxL);
			}
		}
	}

	list<int> nb1, nb2;
	nb1 = vNeighbor;
	for (int r = 0; r < ring; r++)
	{
		nb2.clear();
		for (list<int>::iterator iter = nb1.begin(); iter != nb1.end(); ++iter)
		{
			int idx = *iter;
			for (int l = 0; l < tmesh->getVertex_const(idx)->getValence(); ++l)
			{
				const CHalfEdge* he = tmesh->getVertex_const(idx)->getHalfEdge_const(l);
				int vt = he->getVertexIndex(1);
				if (marked_set.find(vt) == marked_set.end())
				{
					marked_set.insert(vt);
					vNeighbor.push_back(vt);
					nb2.push_back(vt);
				}				
			}
		}
		nb1 = nb2;
	}

	vCoveredIdx.clear();
	for (auto iter = nb1.begin(); iter != nb1.end(); ++iter)
	{
		vCoveredIdx.push_back(*iter);
	}
}

void DiffusionShapeMatcher::readInRandPair( const std::string& filename )
{
	m_randPairs.clear();
	ifstream ifs(filename.c_str());

	double x, y;
	int count(0);
	ifs >> count;
	m_randPairs.reserve(count);
	for (int i = 0; i < count; ++i)
	{
		ifs >> x >> y;
		m_randPairs.push_back(make_pair(x, y));
		
	}
	cout << "## rand pair size: " << m_randPairs.size() << endl;
}

double DiffusionShapeMatcher::evaluateDistortion( const std::vector<MatchPair>& vIdMatchPair, const CMesh* mesh1, const CMesh* mesh2, const std::vector<std::pair<double, double> >& vRandPair, int rand_start /*= 0*/ )
{
	int matchSize = vIdMatchPair.size(), rand_size = vRandPair.size();
	
	double distortSum = 0.;
	double avg_len_ratio = mesh1->getAvgEdgeLength() / mesh2->getAvgEdgeLength();

	const int total_run = 200;
	int count = 0;
	for (int k = rand_start; count < total_run || k < rand_size; ++k)
	{
		int pick1 = matchSize * vRandPair[k].first, pick2 = matchSize * vRandPair[k].second;
		int idx_11 = vIdMatchPair[pick1].m_idx1, idx_12 = vIdMatchPair[pick1].m_idx2;
		int idx_21 = vIdMatchPair[pick2].m_idx1, idx_22 = vIdMatchPair[pick2].m_idx2;

		if (idx_11 == idx_21) continue;

		double dist1 = mesh1->getGeodesic(idx_11, idx_21), dist2 = mesh2->getGeodesic(idx_12, idx_22) * avg_len_ratio;

		distortSum += abs(dist1 - dist2)/dist1;
		count++;
	}

	return distortSum / double(count);
}

double DiffusionShapeMatcher::evaluateDistance( const DifferentialMeshProcessor& mp1, const DifferentialMeshProcessor& mp2, DistanceType distType, const std::vector<double>& vParam, const std::vector<std::pair<double, double> >& vRandPair, int rand_start /*= 0*/ )
{
	int mesh_size = mp1.getMesh_const()->getVerticesNum();
	int rand_size = vRandPair.size();
	
	std::function<double(const DifferentialMeshProcessor&, int, int, const std::vector<double>&)> fDist;

	switch (distType)
	{
	case DISTANCE_GEODESIC:
		fDist = [](const DifferentialMeshProcessor& mp, int v1, int v2, const std::vector<double>& vParam) { 
			return mp.getMesh_const()->getGeodesic(v1, v2) / mp.getMesh_const()->getAvgEdgeLength(); };
		break;
	case DISTANCE_HK:
		fDist = [](const DifferentialMeshProcessor& mp, int v1, int v2, const std::vector<double>& vParam) {
			return mp.calHK(v1, v2, vParam[0]); };
		break;
	case DISTANCE_BIHARMONIC:
		fDist = [](const DifferentialMeshProcessor& mp, int v1, int v2, const std::vector<double>& vParam) {
			return mp.calBiharmonic(v1, v2); };
		break;
	}

	double distort_sum = 0.;
	int count = 0;
	const int total_run = 200;

	for (int k = rand_start; k < rand_start + total_run || k < rand_size; ++k)
	{
		int v1 = mesh_size * vRandPair[k].first, v2 = mesh_size * vRandPair[k].second;
		if (v1 == v2) continue;
		double dist1, dist2;

		dist1 = fDist(mp1, v1, v2, vParam);
		dist2 = fDist(mp2, v1, v2, vParam);

		if (dist1 < 1e-3) 
			cout << "Dist(" << v1 << ',' << v2 << "): " << dist1 << endl;

		distort_sum += abs(dist1 - dist2) / abs(dist1);
		count++;
	}

	return distort_sum / count;
}

void DiffusionShapeMatcher::refineRegister2( std::ostream& flog )
{
	const int current_level = max(m_nAlreadyRegisteredLevel - 1, 0);	//the registration level we currently working on

	flog << "\n-------- registration level " << current_level << " --------" << endl;

	double anchorThresh;
	double regT = 20.0 * pow(2.0, current_level);
	if (current_level == 1) 
	{
		regT = 80;
		anchorThresh = 1e-4;
	}
	if (current_level == 0)
	{
		regT = 20;
		anchorThresh = 2e-4;
	}
	anchorThresh = 1e-4;
	//	thresh = (double)featureMatch.size() * 1e-5;

	CMesh *tmesh1 = meshPyramids[0].getMesh(current_level), *tmesh2 = meshPyramids[1].getMesh(current_level);
	const int coarseSize1 = tmesh1->getVerticesNum(), coarseSize2 = tmesh2->getVerticesNum();
	vector<int> vMatch1(coarseSize1, -1), vMatch2(coarseSize2, -1);
	vector<double> vMatchScore1(coarseSize1, 0), vMatchScore2(coarseSize2, 0);

	const vector<MatchPair>& oldAnchorMatch = getMatchedFeaturesResults(m_nAlreadyMatchedLevel);
	const vector<MatchPair>& oldReg = (current_level == m_nRegistrationLevels - 1) ?
									 oldAnchorMatch : getRegistrationResults(m_nAlreadyRegisteredLevel);

	vector<MatchPair> newAnchorMatch;
	vector<MatchPair> tmpReg1, tmpReg2;

	/************************************************************************/
	/*                        1. Initialize/Compute HKC                     */
	/************************************************************************/
	vector<int> vFeatureID1, vFeatureID2;
	for (auto iter = begin(oldAnchorMatch); iter != end(oldAnchorMatch); ++iter)
	{
		vFeatureID1.push_back(iter->m_idx1);
		vFeatureID2.push_back(iter->m_idx2);
	}

	flog << "Anchors IDs: " << endl;
	for (auto iter = begin(oldAnchorMatch); iter != end(oldAnchorMatch); ++iter)
		flog << "(" << iter->m_idx1 << ',' << iter->m_idx2 << ") ";
	flog << endl << endl;

	CStopWatch timer;
	timer.startTimer();
	m_ParamMgr[0].para_computeHKC(vFeatureID1, regT); // parallel computing supported
	m_ParamMgr[1].para_computeHKC(vFeatureID2, regT); 
//	m_ParamMgr[0].para_computeHKS()
	timer.stopTimer();
	cout << "HKParam computation time: " << timer.getElapsedTime() << 's' << endl;

	std::vector<PointParam>& vCoordinates1 = m_ParamMgr[0].vCoord;
	std::vector<PointParam>& vCoordinates2 = m_ParamMgr[1].vCoord;

	flog << "Multi-anchor coordinates initialized!" << endl;

	if (current_level > 0)
	{
		set<MatchPair> uniquePairs;
		for (vector<MatchPair>::const_iterator mpIter = oldReg.begin(); mpIter != oldReg.end(); ++mpIter)
		{
			const int vid_i = mpIter->m_idx1, vid_j = mpIter->m_idx2;
			const int vi = id2Index(0, vid_i, current_level), vj = id2Index(1, vid_j, current_level);

			// filter out repeated matches
			if (uniquePairs.find(MatchPair(vi, vj)) == uniquePairs.end())
			{
				uniquePairs.insert(MatchPair(vi, vj));
				tmpReg1.push_back(*mpIter);			
			}		
		}
	}
	else tmpReg1 = oldReg;

	for (auto iter = tmpReg1.begin(); iter != tmpReg1.end(); ++iter)
	{
		const int vid_i = iter->m_idx1, vid_j = iter->m_idx2;
		const int vi = id2Index(0, vid_i, current_level), vj = id2Index(1, vid_j, current_level);

		if (find(oldAnchorMatch.begin(), oldAnchorMatch.end(), *iter) != oldAnchorMatch.end())
			iter->m_note = 1;

		if (iter->m_note == 1) 
			iter->m_score = 1.;	// matched feature pairs have the highest score 1.0
		else
		{
			double score = computeMatchScore(iter->m_idx1, iter->m_idx2);
			iter->m_score = score;
		}		

		vMatch1[vi] = vj;
		vMatch2[vj] = vi;
		vMatchScore1[vi] = iter->m_score;
		vMatchScore2[vj] = iter->m_score;
	}

	flog << "\n1. HKC computed!\n" 
		<< "  Old anchor size: " << oldAnchorMatch.size() 
		<< "\n  Old reg size: " << oldReg.size() 
		<< "\n  Filtered old reg size: " << tmpReg1.size() << endl;

	/************************************************************************/
	/*              2. Registration                                         */
	/************************************************************************/
	NoteQueue qscan1;
	for (auto iter = tmpReg1.begin(); iter != tmpReg1.end(); ++iter)
	{
		int vid_1 = iter->m_idx1, vid_2 = iter->m_idx2;
		int v1 = id2Index(0, vid_1, current_level), v2 = id2Index(1, vid_2, current_level);
		//MyNote mn(v1, v2, iter->m_score);	//id, score
		//MyNote mn(v1, tmesh1->m_pVertex[v1].m_vParam.m_votes);
		MyNote mn(v1, v2, vCoordinates1[vid_1].m_votes/* * iter->m_score*/);
		qscan1.push(mn);
	}

	while (!qscan1.empty())
	{
		MyNote qt = qscan1.top();
		qscan1.pop();

		int vi = qt.m_idx1, vj = qt.m_idx2;		//index on this level

		if (vMatch2[vj] != vi) continue;

		int vid_i = tmesh1->getVertex_const(vi)->getVID(),
			vid_j = tmesh2->getVertex_const(vj)->getVID();

		for (int ei = 0; ei < tmesh1->getVertex_const(vi)->getValence(); ei++)
		{
			const CHalfEdge* he = tmesh1->getVertex_const(vi)->getHalfEdge_const(ei);
			const int vt = he->getVertexIndex(1);
			const int vid_t = tmesh1->getVertex_const(vt)->getVID();

			if( vMatch1[vt] >= 0 && vMatch1[vt] < coarseSize2 ) 
				continue;  // already registered

			if (tmesh1->getVertex_const(vt)->isOnBoundary())
				continue;
			//  if( matcher1.m_vHKParamFine[vt].m_votes < thresh )
			//		continue;	// low priority points skipped

			double score;
			int vm = searchVertexMatch(vt, vj, current_level, 2, score);
			if(vm >= 0) 
			{
				int vid_m = tmesh2->getVertex_const(vm)->getVID();
				vMatch1[vt] = vm;
				vMatchScore1[vt] = score;
				tmpReg1.push_back(MatchPair(vid_t, vid_m, score));

				if (vMatch2[vm] == -1 || vMatchScore2[vm] < score)
				{
					vMatch2[vm] = vt;
					vMatchScore2[vm] = score;
				}

				//MyNote mn(vt, vm, score);
				MyNote mn(vt, vm, vCoordinates1[vid_t].m_votes);
				qscan1.push(mn);
			}
			else flog << "No match 3!" << endl;
		}
	} //while(!qscan1.empty())

	flog << "\n2. Crude registration finished!" << "\n  Tmp reg size: " << tmpReg1.size() << endl;

	/************************************************************************/
	/*              3. Adjust and update anchor points                      */
	/************************************************************************/

	newAnchorMatch = oldAnchorMatch;
//	if(current_level > 0)
	{
		for (auto iter = begin(tmpReg1); iter != end(tmpReg1); ++iter)
		{
			int vid_1 = iter->m_idx1, vid_2 = iter->m_idx2;
			int v1 = id2Index(0, vid_1, current_level), v2 = id2Index(1, vid_2, current_level);
			if (vMatch1[v1] == v2 && vMatch2[v2] == v1)		// select only the pair considered best mutually
				tmpReg2.push_back(*iter);
		}

		priority_queue<MatchPair, vector<MatchPair>, greater<MatchPair> > qReg;
		for (auto iter = begin(tmpReg2); iter != end(tmpReg2); ++iter)
		{
			const int vi = iter->m_idx1, vj = iter->m_idx2;
			double score;
			const int vm1 = searchVertexMatch(vi, vj, 0, 2, score, current_level);	// vi, vj here are vertex id as well as level-0 index
			qReg.push(MatchPair(vi, vm1, vCoordinates1[vi].m_votes * score));
		}

		tmpReg1.clear();
		while(!qReg.empty())
		{
			tmpReg1.push_back(qReg.top());
			qReg.pop();
		}

		// insert additional anchor points
		for (auto riter = tmpReg1.begin(); riter != tmpReg1.end(); ++riter)
		{
			if (vCoordinates1[riter->m_idx1].m_votes < 6 * anchorThresh ||
				vCoordinates2[riter->m_idx2].m_votes < 6 * anchorThresh
				) 
				continue;

			bool pass = true;
			for (auto fiter = newAnchorMatch.begin(); fiter != newAnchorMatch.end(); ++fiter)
			{
				if (fiter->m_idx1 == riter->m_idx1 || fiter->m_idx2 == riter->m_idx2 ||
					pOriginalMesh[0]->isInNeighborRing(fiter->m_idx1, riter->m_idx1, 7) ||
					pOriginalMesh[1]->isInNeighborRing(fiter->m_idx2, riter->m_idx2, 7)
					) 
				{ 
					pass = false; 
					break; 
				}
			}
			if (pass)
				newAnchorMatch.push_back(*riter);
			if (newAnchorMatch.size() >= 100) break;
		}
	}

	const vector<MatchPair>& newReg = tmpReg1;

	flog << "\n3. Registration adjustment and anchor set augmentation finished!\n"
		<< "  New reg size: " << newReg.size() << "; new anchor size: " << newAnchorMatch.size() << endl;

	if (m_nAlreadyMatchedLevel > 0) m_nAlreadyMatchedLevel--;
	if (m_nAlreadyRegisteredLevel > 0) m_nAlreadyRegisteredLevel--;
	vFeatureMatchingResults[m_nAlreadyMatchedLevel] = newAnchorMatch;
	vRegistrationResutls[m_nAlreadyRegisteredLevel] = newReg;

}

void DiffusionShapeMatcher::forceInitialAnchors( const std::vector<MatchPair>& mp )
{
	m_bFeatureMatched = true;
	m_nAlreadyMatchedLevel = m_nRegistrationLevels;
	vFeatureMatchingResults[m_nAlreadyMatchedLevel] = mp;
}
