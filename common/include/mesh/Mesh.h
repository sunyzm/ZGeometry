#pragma once
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
#include <util/misc.h>
#include "Geometry.h"
#include "Quat.h"

#define MAX_VERTEX_PER_FACE 20
#define MAX_RING_NUMBER 15
#define MAX_HOLE_SIZE 20

typedef std::vector<int> VectorInt;

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
	int			m_nValence;			// out valence
	int*		m_piEdge;			// half edge indices start from this vertex

	bool		 m_bIsBoundary;     // if boundary vertex
	bool		 m_bIsHole;
	RGBf		 m_vColor;			// vertex color
	double		 m_vMeanCurvature;	// mean curvature
	double		 m_vGaussCurvature;	// Gauss curvature
	
	int			 m_vMatched;

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
	const Vector3D&		getPos() const { return m_vPosition; } 
	bool				judgeOnBoundary();
	bool				isValid() const { return m_bIsValid; }
	void				invalidate(bool flag) { m_bIsValid = flag; }
	void				translateAndScale(Vector3D translation, double s);
	void                setPos( double x, double y, double z );
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
	int	m_iVertex[2];		// starting and ending vertex index Vertex0 ? Vertex1
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

////////////////    attributes    ////////////////
private:
	std::vector<CVertex*>	m_vVertices;
	std::vector<CHalfEdge*> m_vHalfEdges;
	std::vector<CFace*>		m_vFaces;

	bool m_bIsPointerVectorExist;		// pointer vectors representation
	bool m_bIsIndexArrayExist;			// index array representation
	bool m_bSeperateStorage;		
	int  m_nBoundaryEdgeNum;

	int		    m_nVertex;				// number of vertices
	int		    m_nHalfEdge;			// number of half-edges
	int			m_nFace;	 			// number of faces


public:
	CVertex*	m_pVertex;				// array pointer of vertices
	CHalfEdge*	m_pHalfEdge; 			// array pointer of half-edges
	CFace*		m_pFace;				// array pointer of faces

	Vector3D    m_Center;
	Vector3D    m_bBox;
	double		m_edge;					// average edge length
	std::string m_meshName;				// name of the mesh


////////////////    methods    ////////////////
public:
//constructors
	CMesh();
	CMesh(const CMesh* pMesh);
	CMesh(const CMesh& oldMesh);
	virtual ~CMesh();
	void        cloneFrom(const CMesh& oldMesh);
//	void        clone(const CMesh* oldMesh);

//operations
	bool	    Load(std::string sFileName);			// load from file
	bool	    Save(std::string sFileName);			// save to file
	
	void	    move(Vector3D translation);
	void	    scaleAreaToVertexNum();					// move to origin point and scale the mesh so that the surface area equals vertex num
	void        scaleEdgeLenToUnit();
	void	    gatherStatistics();

	CVertex*	getVertex(int i) { return &m_pVertex[i]; }
	const CVertex* getVertex_const(int i) const { return &m_pVertex[i]; }
	CFace*		getFace(int i) { return &m_pFace[i]; }
	double		getHalfEdgeLen(int iEdge) const;				// get the Euclidean length of the iEdge-th half-edge
	int			getVerticesNum() const { return m_nVertex; }
	int			getFaceNum() const { return m_nFace; }
	int         getHalfEdgeNum() const { return m_nHalfEdge; }
	int			getMatchedVerticesNum() const;
	int			getFaceNum() {return m_nFace;}
	double		getAvgEdgeLength() const { return m_edge; }
	int			getEdgeNum();		// get number of edges ( not half-edge! )
	int			getBoundaryNum() const;    // get number of boundaries
	bool        hasBounary() const;
	int			getBoundaryVertexNum() const; // get number of boundary vertices
	Vector3D	getBoundingBox() const { return m_bBox; }
	int			getEulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int			getMeshGenus();			// get mesh genus
	double		calGaussianCurvatureIntegration();	// compute the integration of Gaussian curvature over all vertices
	bool		calVertexLBO(int i, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, double& Av, std::vector<double>& tw) const;
	
//  bool	    CalVertexLBOs(int i, double h, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, std::vector<double>& Av, std::vector<double>& tw);
	bool		calVertexArea(std::vector<double>& Av);
	bool		VertexNeighborGeo(int i, double ring, std::vector<GeoNote>& nbg); // geodesic vertex neighbor
	void		VertexNeighborRing(int i, int ring, std::vector<int>& nbr) const;
	bool		isInNeighborRing(int ref, int query, int ring) const;
	double		getGeodesic(int s, int t) const;
	double		getGeodesicToBoundary(int s) const;	// return 0.0 if in a manifold
	double		getGeodesicToBoundary(int s, std::vector<GeoNote>& nbg);
	double		getVolume();	// calculate volume (area) of a surface
	bool		calAreaRatio(CMesh* tmesh, std::vector<int>& ar);	// for registration
	bool		calLengthDifference(const CMesh* tmesh, std::vector<double>& ld) const;
	bool		calVertexCurvature( int i );			// calculate number i-th vertex's Gaussian and mean curvature
	void		clearVertexMark();
	void		extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures ) const;
	VectorInt	getOriginalVertexIndex() const;
	VectorInt	getNeighboringVertex(int v, int ring) const;
	VectorInt   getRingVertex(int v, int ring) const;
	VectorInt	getVertexAdjacentFacesIndex(int vIdx);
	void        getCoordinateFunction(int dim, std::vector<double>& vCoord) const;
	void        setVertexCoordinates(const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord);
	void		setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<Vector3D>& vNewPos);
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
	double	calAreaMixed(double a, double b, double c, double& cotan_a, double& cotan_c);
	double  calHalfAreaMixed(double a, double b, double c, double& cotan_a) const;
	double  calLocalGeodesic(int ia, int ib, int ic) const;
	void	findHoles();

	void	buildPointerVectors();		//construct vectors of pointers based on array representations
	void	buildIndexArrays();	// already have the pointer-vector representation; fill in the array represenatation
	void	assignElementsIndex();
	bool	isHalfEdgeMergeable(const CHalfEdge* halfEdge);
};	//CMesh
