#pragma once
#include <vector>
#include <algorithm>
#include <fstream>
#include <engine.h>
#include "DifferentialMeshProcessor.h"
#include "MeshPyramid.h"


class MatchPair
{
public:
	int		m_idx1;
	int		m_idx2;
	double	m_tl;		//lower time
	double	m_tu;		//upper time
	int		m_tn;
	double  m_score;
	int		m_note;		

public:
	MatchPair() { m_idx1 = -1; m_idx2 = -1; }
	MatchPair(int i1, int i2, double score = 0) { m_idx1 = i1; m_idx2 = i2; m_score = score; m_tl = 0.0; m_tn = 0;  m_note = 0;}
	MatchPair(int i1, int i2, double tl, int tn, double score = 0) { m_idx1 = i1; m_idx2 = i2; m_tl = tl; m_tn = tn; m_score = score; m_note = 0;}
//	bool operator== (const MatchPair& mc) const { return (m_idx1 == mc.m_idx1 && m_idx2 == mc.m_idx2); }
	friend bool operator== (const MatchPair& mp1, const MatchPair& mp2) { return (mp1.m_idx1 == mp2.m_idx1 && mp1.m_idx2 == mp2.m_idx2); }
 	friend bool operator< (const MatchPair& mp1, const MatchPair& mp2)
 	{
 		return mp1.m_score < mp2.m_score;
 	}
 
 	friend bool operator> (const MatchPair& mp1, const MatchPair& mp2)
 	{
 		return mp1.m_score > mp2.m_score;
 	}
};

class HKSFeature : public MeshFeature
{
public:
	double	m_tl;	    // lower time
	double	m_tu;	    // upper time
	int		m_tn;	    // num of scales; log_2(tu/tl)
public:
	HKSFeature() { m_index = 0; m_tu = 0.0; m_tl = 0.0; m_tn = 0; }
	HKSFeature(int i, int s) : MeshFeature(i, s) { HKSFeature(); }
	HKSFeature(int x, double tl, double tu) { m_index = x; m_tl = tl; m_tu = tu; m_tn = 0; }
	HKSFeature(int x, double tl, int n) { m_index = x; m_tn = n; m_tu = tl; }
	HKSFeature(int x, double tl, double tu, int n) { m_index = x; m_tn = n; m_tl = tl; m_tu = tu; }
	HKSFeature(int x, double tl, double tu, int n, int scale) { m_index = x; m_tn = n; m_tl = tl; m_tu = tu; m_scale = scale; }
	HKSFeature(const HKSFeature& f) { m_index = f.m_index; m_scale = f.m_scale; m_tl = f.m_tl; m_tu = f.m_tu; m_tn = f.m_tn; }
	void setTimes(double tl, double tu, double tn) { m_tl = tl; m_tu = tu; m_tn = tn; }
};

class HKParam : public VectorND
{
public:
	HKParam& operator =(const HKParam& hkp);
	double m_votes;
	virtual void clear();
};

class HKParamManager
{
public:
	DifferentialMeshProcessor* pMP;
	std::vector<HKParam> vHKParam;
	
	void initialize(DifferentialMeshProcessor* p) { pMP = p; }
	void computeHKParam(const std::vector<int>& anchors, double t = 30.0);
	const std::vector<HKParam>& getHKParam() const { return vHKParam; }
};

class DiffusionShapeMatcher
{
public:
	class Cluster{
	public:
		std::vector<int> m_member;
	};
public:
	DiffusionShapeMatcher();
	~DiffusionShapeMatcher();

	/* core functions */
	void    initialize(DifferentialMeshProcessor* pMP1, DifferentialMeshProcessor* pMP2, Engine *ep);
	void	constructPyramid(int n);
	void	detectFeatures(int obj, int ring = 2, int scale = 1, double tvalue = DEFAULT_FEATURE_TIMESCALE, double talpha = DEFAULT_T_MULTIPLIER, double thresh = DEFAULT_EXTREAMA_THRESH);
	void    matchFeatures(std::ofstream& flog, double matchThresh = DEFAULT_MATCH_THRESH);
	void    matchFeaturesTensor(std::ofstream& flog, double timescale, double thresh);
	void    refineRegister(std::ofstream& flog);
	void	evaluateRegistration();

	/* attributes access */
	DifferentialMeshProcessor* getMeshProcessor(int obj, int level) { return liteMP[obj].at(level); }
	int     getPyramidLevels() const { return m_nPyramidLevels; }
	bool	isPyramidBuilt() const { return m_bPyramidBuilt; }
	int		getTotalRegistrationLevels() const { return m_nRegistrationLevels; }
	void	setRegistrationLevels(int val);
	void	setEngine(Engine* ep) { m_ep = ep; }
	const MeshPyramid& getMeshPyramid(int obj) const { return meshPyramids[obj]; }
	CMesh*  getMesh(int obj, int level = 0) const;
	int     getAlreadyRegisteredLevel() const { return m_nAlreadyRegisteredLevel; }
	int		getAlreadyMatchedLevel() const { return m_nAlreadyMatchedLevel; }
	const std::vector<MatchPair>& getMatchedFeaturesResults(int level) const;
	const std::vector<MatchPair>& getRegistrationResults(int level) const;
	const std::vector<HKSFeature>& getSparseFeatures(int obj) const { return vFeatures[obj]; }

	int		id2Index(int obj, int vid, int level) const { return meshPyramids[obj].m_Id2IndexMap[vid][level]; }
	void    dumpIndexMap(const std::string& filename) const;

	// static constants
	static const double DEFAULT_C_RATIO;
	static const double DEFAULT_RANK_EPSILON;
	static const double SPARSIFY_EPSILON;
	static const double DEFAULT_FEATURE_TIMESCALE;
	static const double DEFAULT_T_MULTIPLIER;
	static const double DEFAULT_MATCH_TIME_LOW;
	static const double DEFAULT_MATCH_THRESH;
	static const double DEFAULT_EXTREAMA_THRESH;
	static const double DEFAULT_REGISTER_TIMESCALE;
	static const int	DEFAULT_PYRAMID_LEVELS;
	static const int	MAXIMAL_PYRAMID_LEVELS;
	static const int    NUM_OF_EIGVAL_FOR_ESTIMATE;

private:
	Engine* m_ep;
	CMesh* pOriginalMesh[2];
	DifferentialMeshProcessor* pOriginalProcessor[2];
	MeshPyramid meshPyramids[2];
	std::vector<DifferentialMeshProcessor*> liteMP[2];
	std::vector<HKSFeature> vFeatures[2];	// original detected fine features
	std::vector<std::vector<MatchPair> > vFeatureMatchingResults;
	std::vector<std::vector<MatchPair> > vRegistrationResutls;

	bool					m_bPyramidBuilt;
	bool					m_bFeatureDetected;
	bool					m_bFeatureMatched;
	int						m_nPyramidLevels;           // >= 1
	int						m_nRegistrationLevels;      // [1,..,m_nPyramidLevels]
	int						m_nAlreadyRegisteredLevel;	// [0,..,m_nRegistrationLevels]
	int						m_nAlreadyMatchedLevel;		// [1,...,m_nRegistrationLevels]
	int						m_nBaseEigensMatch, m_nBaseEigensRegister;
	double					m_registerTimescale;
	HKParamManager			m_HKParamMgr[2];


	/* helper functions */
	static void	calVertexSignature( const DifferentialMeshProcessor* pOriginalProcessor, const HKSFeature& hf, VectorND& sig );
	static void ComputeTensorFeature( const DifferentialMeshProcessor* pmp, int i, int j, int k, double t, double* sang);
	static double TensorMatching(Engine *ep, const DifferentialMeshProcessor* pmp1,  const DifferentialMeshProcessor* pmp2, Cluster& ct1, Cluster& ct2, std::vector<MatchPair>& matched, double t, double thresh);
	void    prepareHeatRegistration( double regTime );
	double  computeMatchScore(int idx1, int idx2) const;
	int		searchVertexMatch( const int vt, const int vj, const int level, const int ring, double& score, int uppper_level = -1 );
	void    getVertexCover(int obj, int vidx, int level, int upper_level, int ring,  std::vector<int>& vCoveredIdx) const;
};