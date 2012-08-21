#ifndef _MESH_H_
#define _MESH_H_

#include <cstdio>
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <string>
#include <algorithm>
#include <util/util.h>
#include "Geometry.h"
#include "Quat.h"

#define MAX_VERTEX_PER_FACE 20
#define MAX_RING_NUMBER 15
#define MAX_HOLE_SIZE 20

class MyNote{
public:
	int m_idx1, m_idx2;			//here id is vertex index on a particular level
	double m_score;
public:
	MyNote(int mid, double s) { m_idx1 = mid; m_score = s; }
	MyNote(int mid1, int mid2, double s) { m_idx1 = mid1; m_idx2 = mid2; m_score = s; }
	MyNote& operator = (const MyNote& note) { m_idx1 = note.m_idx1; m_idx2 = note.m_idx2; m_score = note.m_score; return(*this); }
};

class NoteCompare {
public:
	bool operator()(const MyNote& Left, const MyNote& Right) const { return ( Left.m_score < Right.m_score ); }
};

typedef std::priority_queue<MyNote, std::vector<MyNote>, NoteCompare> NoteQueue;

class GeoNote{
public:
	int m_id;
	double m_geodesic;
public:
	GeoNote(int mid, double geo) {m_id = mid; m_geodesic = geo;}
	GeoNote& operator = (const GeoNote& note) { m_id=note.m_id; m_geodesic = note.m_geodesic; return(*this); }
};

class GeoCompare {
public:
	bool operator()(const GeoNote& Left, const GeoNote& Right) const { return ( Left.m_geodesic > Right.m_geodesic ); }
};

typedef std::priority_queue<GeoNote, std::vector<GeoNote>, GeoCompare> GeoQueue;

class HKParam : public VectorND
{
public:
	HKParam& operator =(const HKParam& hkp);
	double m_votes;
	virtual void clear();
};

class CFace;
class CHalfEdge;
class CMesh;

typedef int VID;
//////////////////////////////////////////////////////
//						CVertex  					//
//////////////////////////////////////////////////////
class CVertex
{
public:
	friend class CMesh;
	friend class MeshPyramid;

private:
	VID						m_vid;
	int						m_vIndex;
	std::vector<CHalfEdge*> m_HalfEdges;
	Vector3D				m_vNormal;         // vertex normal
	bool					m_bIsValid;
public:
	Vector3D	m_vPosition;		// vertex coordinates
	int			m_nValence;			// vertex valence
	int*		m_piEdge;			// half edge indices start from this vertex

	bool		 m_bIsBoundary;     // if boundary vertex
	bool		 m_bIsHole;
	RGBf		 m_vColor;			// vertex color
	double		 m_vMeanCurvature;	// mean curvature
	double		 m_vGaussCurvature;	// Gauss curvature
	
	int			 m_vMatched;
	double       m_matchScore;

	HKParam		 m_vParam;			// this is VectorND + votes, likely to be HKC norm
	double		 m_LocalGeodesic;	// geodesic from local vertex
	bool		 m_inheap;			// in heap or not

	int m_mark;

public:
	//constructors
	CVertex();
	CVertex(double x, double y, double z);
	CVertex(const Vector3D& v);
	CVertex(const CVertex& v);
	CVertex(double x, double y, double z, float r, float g, float b);
	CVertex& operator = (CVertex& v);	
	//destructor
	virtual ~CVertex();
	
	//operations
	int					getVID() const { return m_vid; }
	int					getIndex() const { return m_vIndex; }
	Vector3D			getNormal() const { return m_vNormal; }
	void				calcNormal();
	std::vector<CFace*> getAdjacentFaces();
	const Vector3D&		getPos() const {return m_vPosition;} 
	bool				judgeOnBoundary();
	bool				isValid() const {return m_bIsValid;}
	void				invalidate(bool flag) {m_bIsValid = flag;}
	void				translateAndScale(Vector3D translation, double s);

public:
	std::list<int>	m_lEdgeList;      // temporary list of constructing m_piEdge
};

//////////////////////////////////////////////////////
//						CHalfEdge					//
//////////////////////////////////////////////////////
class CHalfEdge
{
public:
	friend class CMesh;
	friend class MeshPyramid;
private:
	CVertex*	m_Vertices[2];	//starting and ending vertices
	CHalfEdge*	m_eTwin;		//reverse half-edge; null if boundary half edge
	CHalfEdge*	m_eNext;		//next half-edge (counterclockwise)
	CHalfEdge*	m_ePrev;
	CFace*		m_Face;			//attached face
	int			m_eIndex;		//half-edge id
	bool		m_bIsValid;
public:
	int	m_iVertex[2];		// starting and ending vertex index Vertex0 ­> Vertex1
	int	m_iTwinEdge;		// reverse half-edge index, -1 if boundary half edge
	int	m_iNextEdge;		// next half-edge index ( counter-clock wise )
	int m_iPrevEdge;
	int	m_iFace;            // attaching face index ( on the left side )
public:
	//constructions
	CHalfEdge();
	CHalfEdge(const CHalfEdge& oldE);
	CHalfEdge(int iV0, int iV1);
	virtual ~CHalfEdge();

	//operations
	CHalfEdge& operator = (const CHalfEdge& e);
	CFace* getAttachedFace(){return m_Face;}
	CHalfEdge* getTwinHalfEdge(){return m_eTwin;}
	double getLength();
	bool	isValid() { return m_bIsValid; }
};

//////////////////////////////////////////////////////
//						CFace   					//
//////////////////////////////////////////////////////
class CFace
{
public:
	friend class CMesh;
	friend class MeshPyramid;
private:
	Vector3D	m_vNormal;			// normalized face normal
	double		m_faceArea;			// face area
	int			m_fIndex;
	std::vector<CVertex*> m_Vertices;	//all vertices
	std::vector<CHalfEdge*> m_HalfEdges;//all half-edges
	bool		m_bIsValid;
public:
	short		m_nType;			// number of polygon face edges
	int*		m_piVertex;			// all vertex index
	int*		m_piEdge;			// all half-edge index

public:
	//constructions
	CFace();
	CFace(short s);
	CFace(const CFace& oldF);
	virtual ~CFace();
	CFace& operator = (const CFace& f);

	//operations
	void Create(short s);
	std::vector<double> getPlaneFunction();	
	CVertex* getVertex(int i) { return m_Vertices[i]; }
	double getArea();
	void calcNormalAndArea();
	Vector3D getNormal() { return m_vNormal; }
	bool hasVertex(int vid);
	double distanceToVertex(const CVertex* vq, std::vector<double>& baryCoord);
};

//////////////////////////////////////////////////////
//						CMesh   					//
//////////////////////////////////////////////////////
class CMesh 
{
	friend class MeshPyramid;
//attributes
private:
	std::vector<CVertex*>	m_Vertices;
	std::vector<CHalfEdge*> m_HalfEdges;
	std::vector<CFace*>		m_Faces;

	bool m_bIsPointerVectorExist;		// whether vertex, HalfEdge, Face data is local
	bool m_bIsArrayRepresentationExist;
	int  m_nBoundaryEdgeNum;
public:
	int		    m_nVertex;				// number of vertices
	CVertex*	m_pVertex;				// array pointer of vertices
	int		    m_nHalfEdge;				// number of half-edges
	CHalfEdge*	m_pHalfEdge; 				// array pointer of half-edges
	int			m_nFace;	 			// number of faces
	CFace*		m_pFace;				// array pointer of faces

	Vector3D    m_Center;
	Vector3D    m_bBox;
	double		m_edge;					// average edge length
	std::string m_meshName;				// name of the mesh

public:
//constructors
	CMesh();
	CMesh(const CMesh* pMesh);
	CMesh(const CMesh& oldMesh);
	virtual ~CMesh();
	void        clone(const CMesh& oldMesh);
	void        clone(const CMesh* oldMesh);

//operations
	bool	    Load(std::string sFileName);			// load from file
	bool	    Save(std::string sFileName);			// save to file
	
	void	    move(Vector3D translation);
	void	    scaleAreaToVertexNum();					// move to origin point and scale the mesh so that the surface area equals vertex num
	void        scaleEdgeLenToUnit();
	void	    gatherStatistics();

	CVertex*	getVertex(int i) { return m_Vertices[i]; }
	CFace*		getFace(int i) { return &m_pFace[i]; }
	double		getHalfEdgeLen(int iEdge) const;				// get the Euclidean length of the iEdge-th half-edge
	int			getVerticesNum() const { return m_nVertex; }
	int			getMatchedVerticesNum() const;
	int			getFaceNum() {return m_nFace;}
	double		getAvgEdgeLength() const { return m_edge; }
	int			GetEdgeNum(  );		// get number of edges ( not half-edge! )
	int			getBoundaryNum( ) const;    // get number of boundaries
	bool        hasBounary() const;
	int			getBoundaryVertexNum(  ) const; // get number of boundary vertices
	int			EulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int			Genus( );			// get mesh genus
	double		CalGaussianCurvatureIntegration(  );	// compute the integration of Gaussian curvature over all vertices
	bool		CalVertexLBO(int i, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, double& Av, std::vector<double>& tw) const;
//  bool	    CalVertexLBOs(int i, double h, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, std::vector<double>& Av, std::vector<double>& tw);
	bool		CalVertexArea(std::vector<double>& Av);
	bool		VertexNeighborG(int i, double ring, std::vector<GeoNote>& nbg); // geodesic vertex neighbor
	bool		VertexNeighborR(int i, int ring, std::vector<int>& nbr) const;
	bool		isInNeighborR(int ref, int query, int ring) const;
	double		CalGeodesic(int s, int t) const;
	double		CalGeodesicToBoundary(int s) const;	// return 0.0 if in a manifold
	double		CalGeodesicToBoundary(int s, std::vector<GeoNote>& nbg);
	double		CalVolume();	// calculate volume (area) of a surface
	bool		CalAreaRatio(CMesh* tmesh, std::vector<int>& ar); // for registration
	bool		CalLengthDifference(const CMesh* tmesh, std::vector<double>& ld) const;
	void		ClearVertexMark();
	void		buildArrayRepresentation();	//already have the pointer-vector representation; fill in the array represenatation
	void		extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures ) const;
	std::vector<int> getOriginalVertexIndex() const;
	std::vector<int> getNeighboringVertex(int v, int ring) const;
	std::vector<int> getVertexAdjacentFacesIndex(int vIdx);

private:
	void	clear();
	bool	construct();	// construct connectivity
	bool	loadFromOBJ(std::string sFileName);	// load mesh from .obj file
	bool	loadFromM(std::string sFileName);	// load mesh from .m file
	bool	loadFromVERT(std::string sFileName); // load mesh from .vert + .tri files
	bool	loadFromPLY(std::string sFileName);	// load mesh from .ply files
	bool	loadFromOFF(std::string sFileName);
	bool	saveToOBJ(std::string sFileName);	// save mesh to .obj file
	bool    saveToM( std::string sFileName );    // save mesh to .m file

	void	calFaceNormalAndArea(int i);			// compute i-th face's normal
	void	calVertexNormal(int i);					// compute i-th vertex's normal
	bool	calVertexCurvature( int i );			// calculate number i-th vertex's Gaussian and mean curvature
	double	calAreaMixed(double a, double b, double c, double& cotan_a, double& cotan_c);
	double  calHalfAreaMixed(double a, double b, double c, double& cotan_a) const;
	double  calLocalGeodesic(int ia, int ib, int ic) const;
	void	findHoles();

	void	buildConnectivity();		//construct vectors of pointers based on array representations
	void	assignElementsIndex();
	bool	isHalfEdgeMergeable(const CHalfEdge* halfEdge);
};	//CMesh


#endif