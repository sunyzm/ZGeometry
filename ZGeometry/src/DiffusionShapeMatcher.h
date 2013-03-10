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
	MatchPair(int i1, int i2, double score = 0) { m_idx1 = i1; m_idx2 = i2; m_score = score; m_tl = 0.0; m_tn = 0;  m_note = 0;}
	MatchPair(int i1, int i2, double tl, int tn, double score = 0) { m_idx1 = i1; m_idx2 = i2; m_tl = tl; m_tn = tn; m_score = score; m_note = 0;}
	bool operator== (const MatchPair& mc) const;

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

class DiffusionShapeMatcher
{
public:
	DiffusionShapeMatcher();
	~DiffusionShapeMatcher();

	/* core functions */
	void    initialize(DifferentialMeshProcessor* pMP1, DifferentialMeshProcessor* pMP2, Engine *ep);
	void	constructPyramid(int n);
	void	detectFeatures(int obj, int ring = 2, int scale = 1, double tvalue = DEFAULT_FEATURE_TIMESCALE, double talpha = DEFAULT_T_MULTIPLIER, double thresh = DEFAULT_EXTREAMA_THRESH);
	void    matchFeatures(std::ofstream& flog, double matchThresh = DEFAULT_MATCH_THRESH);
	void    registerOneLevel(std::ofstream& flog);
	void	evaluateRegistration();

	/* attributes access */
	const std::vector<MatchPair>& getFeatureMatches(int level) const;
	int     getPyramidLevels() const { return m_nPyramidLevels; }
	DifferentialMeshProcessor* getMeshProcessor(int obj, int level) { return liteMP[obj].at(level); }

	/* basic */
	void	setEngine(Engine* ep) { m_ep = ep; }
	const MeshPyramid& getMeshPyramid(int obj) const { return meshPyramids[obj]; }
	CMesh*  getMesh(int obj, int level = 0) const;
	const std::vector<MatchPair>& getMatchedFeatures();

	bool	isPyramidBuilt() const { return m_bPyramidBuilt; }
	int		getRequiredLevels() const { return m_nRequiredLevels; }
	void	setRequiredLevels(int val) { m_nRequiredLevels = min(val, m_nPyramidLevels); }
	void    dumpIndexMap(const std::string& filename) const;
	int		id2Index(int obj, int vid, int level) const { return meshPyramids[obj].m_Id2IndexMap[vid][level]; }
	
	double	getHeatKernelValue( int vi, int vj, int level);
	void	initializeCoarseKernel(double regTimescale, const std::string& kernelName, const std::string& prolongMatName);	//calculate coarse heat kernel matrix on the given t, save the result into Matlab	
	void	setHKParam(int v, const std::vector<int>& anchors, int level) const;
	void    initializeHKParam( const std::vector<int>& anchors, double t = 30.0 );
	void	setRegisterDiffusionTime(double t) { this->m_registerTimescale = t; }

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
	static const int    DEFAULT_NBASE;
	static const int	DEFAULT_PYRAMID_LEVELS;
	static const int	MAXIMAL_PYRAMID_LEVELS;
	static const int    NUM_OF_EIGVAL_FOR_ESTIMATE;

private:
	Engine *m_ep;
	CMesh* pOriginalMesh[2];
	DifferentialMeshProcessor* pOriginalProcessor[2];
	MeshPyramid meshPyramids[2];
	std::vector<DifferentialMeshProcessor*> liteMP[2];
	std::vector<HKSFeature> vFeatures[2];	// original detected fine features

	std::vector<HKSFeature> m_hksFeatureFine, m_hksFeatureCoarse;

	bool					m_bPyramidBuilt;
	bool					m_bFeatureDetected;
	bool					m_bFeatureMatched;
	int						m_nPyramidLevels;
	int						m_nRequiredLevels;
	int						m_nCurrentMatchLevel;	// the level of mesh that have been registered. 

	int						m_nBaseEigensMatch, m_nBaseEigensRegister;
	double					m_registerTimescale;
	std::vector<std::vector<int> > m_matchCandidates;
	double *randArray;	


	std::vector<MatchPair> matchedPairsCoarse, matchedPairsFine;

	/* helper functions */
	void    prepareHeatRegistration( double regTime );
	void	calVertexSignature( const DifferentialMeshProcessor* pOriginalProcessor, const HKSFeature& vftCoarse1, VectorND& vsig1 ) const;

};