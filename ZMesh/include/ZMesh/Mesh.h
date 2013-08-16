#ifndef ZMESH_MESH_H
#define ZMESH_MESH_H

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
#include <functional>
#include <unordered_map>
#include <ZUtil/color.h>
#include "Geometry.h"
#include "Quat.h"
#include "MeshAttr.h"

const int MAX_VERTEX_PER_FACE = 20;
const int MAX_RING_NUMBER = 15;
const int MAX_HOLE_SIZE = 20;

typedef std::vector<int> VectorInt;

class GeoNote
{
public:
	int m_id;
	double m_geodesic;
public:
	GeoNote(int mid, double geo) {m_id = mid; m_geodesic = geo;}
	GeoNote& operator = (const GeoNote& note) { m_id=note.m_id; m_geodesic = note.m_geodesic; return(*this); }
    friend bool operator > (const GeoNote& note1, const GeoNote& note2) { return note1.m_geodesic > note2.m_geodesic; }
};

class GeoCompare 
{
public:
	bool operator()(const GeoNote& Left, const GeoNote& Right) const { return ( Left.m_geodesic > Right.m_geodesic ); }
};

typedef std::priority_queue<GeoNote, std::vector<GeoNote>, /*GeoCompare*/std::greater<GeoNote> > GeoQueue;

class CFace;
class CHalfEdge;
class CMesh;

//////////////////////////////////////////////////////
//						CVertex  					//
//////////////////////////////////////////////////////
class CVertex
{
public:
	friend class CHalfEdge;
	friend class CFace;
	friend class CMesh;
	friend class MeshPyramid;
	
	// ---- constructors ---- //
	CVertex();
	CVertex(double x, double y, double z);
	CVertex(const Vector3D& v);
	CVertex(const CVertex& v);
	CVertex(double x, double y, double z, float r, float g, float b);
	CVertex& operator = (const CVertex& v);	
	virtual ~CVertex();

	// ---- operations ---- //
	void				calcNormal();
	int					getIndex() const { return m_vIndex; }
	int					getVID() const { return m_vid; }
	Vector3D			getNormal() const { return m_vNormal; }
	double				getMeanCurvature() const { return m_vMeanCurvature; }
	double				getGaussCurvature() const { return m_vGaussCurvature; }
	std::vector<const CFace*> getAdjacentFaces() const;
	CHalfEdge*			getHalfEdge(int ei) { return m_HalfEdges[ei]; }
	const CHalfEdge*    getHalfEdge_const(int ei) const { return m_HalfEdges[ei]; }
	const Vector3D&		getPosition() const { return m_vPosition; } 
	int					getValence() const { return m_nValence; }
	bool				judgeOnBoundary();
	bool				isOnBoundary() const { return m_bIsBoundary; }
	bool				isHole() const { return m_bIsHole; }
	bool				isValid() const { return m_bIsValid; }
	void				invalidate(bool flag) { m_bIsValid = flag; }
	void				translateAndScale(Vector3D translation, double s);
	void                setPosition( double x, double y, double z );

private:
	void clone(const CVertex& v);
	// ---- fields ---- //
	int						m_vIndex;
	int						m_vid;				// ID of the vertex from original mesh
	std::vector<CHalfEdge*> m_HalfEdges;		// all half-edges from the vertex
	int*					m_piEdge;			// half edge indices start from this vertex
	int						m_nValence;		    // out valence
	Vector3D				m_vPosition;		// vertex coordinates
	Vector3D				m_vNormal;          // vertex normal
	RGBf					m_vColor;			// vertex color
	bool					m_bIsValid;
	bool					m_bIsBoundary;      // if boundary vertex
	bool					m_bIsHole;
	double					m_vMeanCurvature;	// mean curvature
	double					m_vGaussCurvature;	// Gauss curvature

	double					m_LocalGeodesic;	// geodesic from local vertex
	bool					m_inheap;			// in heap or not
	int						m_mark;
};

//////////////////////////////////////////////////////
//						CHalfEdge					//
//////////////////////////////////////////////////////
class CHalfEdge
{
public:
	friend class CVertex;
	friend class CFace;
	friend class CMesh;
	friend class MeshPyramid;

	// -- constructors -- //
	CHalfEdge();
	CHalfEdge(const CHalfEdge& oldE);
	CHalfEdge(int iV0, int iV1);
	virtual ~CHalfEdge();

	// -- operations -- //
	CHalfEdge& operator = (const CHalfEdge& e);
	CFace*		getAttachedFace() { return m_Face; }
	const CFace* getAttachedFace_const() const { return m_Face; }
	CHalfEdge*  getTwinHalfEdge() { return m_eTwin; }
	bool		isBoundaryEdge() const { return (m_iTwinEdge == -1); }
	int			getVertexIndex(int i) const { return m_iVertex[i]; }
	double		getLength();
	bool		isValid() const { return m_bIsValid; }
	int         getIndex() const { return m_eIndex; }
private:
	// -- fields -- //
	CVertex*	m_Vertices[2];	//starting and ending vertices
	CHalfEdge*	m_eTwin;		//reverse half-edge; null if boundary half edge
	CHalfEdge*	m_eNext;		//next half-edge (counterclockwise)
	CHalfEdge*	m_ePrev;
	CFace*		m_Face;			//attached face
	int			m_eIndex;		//half-edge id
	bool		m_bIsValid;

	int			m_iVertex[2];	// starting and ending vertex index
	int			m_iTwinEdge;	// reverse half-edge index, -1 if boundary half edge
	int			m_iNextEdge;	// next half-edge index ( counter-clock wise )
	int			m_iPrevEdge;	// previous half-edge index
	int			m_iFace;        // attaching face index ( on the left side )
};

//////////////////////////////////////////////////////
//						CFace   					//
//////////////////////////////////////////////////////
class CFace
{
public:
	friend class CVertex;
	friend class CHalfEdge;
	friend class CMesh;
	friend class MeshPyramid;

	// -- constructors -- //
	CFace();
	CFace(short s);
	CFace(const CFace& oldF);
	virtual ~CFace();
	CFace& operator= (const CFace& f);

	// -- operations -- //
	void					Create(short s);
	std::vector<double>		getPlaneFunction();	
	CVertex*				getVertex(int i) { return m_Vertices[i]; }
	const CVertex*			getVertex_const(int i) const { return m_Vertices[i]; }
	int						getVertexIndex(int i) const { return m_Vertices[i]->getIndex(); }
	double					getArea() const { return m_faceArea; }
	double					computeArea() const { return TriArea(m_Vertices[0]->getPosition(), m_Vertices[1]->getPosition(), m_Vertices[2]->getPosition()); }
	void					calcNormalAndArea();
	Vector3D				getNormal() const { return m_vNormal; }
	bool					hasVertex(int vidx) const;
	bool					hasHalfEdge() const { return (m_piEdge != NULL); }
	double					distanceToVertex(const CVertex* vq, std::vector<double>& baryCoord);
	int						getFaceIndex() const { return m_fIndex; }
private:
	// ---- fields ---- // 
	Vector3D				m_vNormal;		// normalized face normal
	double					m_faceArea;		// face area
	int						m_fIndex;
	std::vector<CVertex*>	m_Vertices;		//all vertices
	std::vector<CHalfEdge*> m_HalfEdges;	//all half-edges
	bool					m_bIsValid;

	short					m_nType;		// number of polygon face edges
	int*					m_piVertex;		// all vertex index
	int*					m_piEdge;		// all half-edge index
};

//////////////////////////////////////////////////////
//						CMesh   					//
//////////////////////////////////////////////////////
class CMesh 
{
	friend class MeshPyramid;

////////////////   fields    ////////////////
private:
	std::vector<CVertex*>	m_vVertices;
	std::vector<CHalfEdge*> m_vHalfEdges;
	std::vector<CFace*>		m_vFaces;

	bool		m_bIsPointerVectorExist;		// pointer vectors representation
	bool		m_bIsIndexArrayExist;			// index array representation
	bool		m_bSeparateStorage;		
	int			m_nBoundaryEdgeNum;

	int		    m_nVertex;				// number of vertices
	int		    m_nHalfEdge;			// number of half-edges
	int			m_nFace;	 			// number of faces

	Vector3D    m_Center;
	Vector3D    m_bBox;
	double		m_avgEdgeLen;			// average edge length

	CVertex*	m_pVertex;				// array pointer of vertices
	CHalfEdge*	m_pHalfEdge; 			// array pointer of half-edges
	CFace*		m_pFace;				// array pointer of faces

	std::string m_meshName;				// name of the mesh
    std::unordered_map<std::string, MeshAttrBase*> mAttributes;

////////////////    methods    ////////////////
public:
	/* ---- constructors ---- */
	CMesh();
	CMesh(const CMesh* pMesh);
	CMesh(const CMesh& oldMesh);
	virtual ~CMesh();

	/* ---- Mesh IO and processing ---- */
	void        cloneFrom(const CMesh& oldMesh);
	void		cloneFrom(const CMesh* oldMesh);
	bool	    Load(std::string sFileName);			// load from file
	bool	    Save(std::string sFileName);			// save to file
	void	    move(Vector3D translation);				// translate mesh
	void	    scaleAreaToVertexNum();					// move to origin and scale the mesh so that the surface area equals number of vertices
	void        scaleEdgeLenToUnit();					// move to origin and scale the mesh so that the average edge length is 1

	/* ---- attributes access ---- */
	const std::string&	getMeshName() const { return m_meshName; }
	CVertex*			getVertex(int i) { return m_vVertices[i]; }
	const CVertex*		getVertex_const(int i) const { return m_vVertices[i]; }
	CFace*				getFace(int i) { return &m_pFace[i]; }
	const CFace*		getFace_const(int i) const { return m_vFaces[i]; }
	const CHalfEdge*	getHalfEdge_const(int i) const { return m_vHalfEdges[i]; }
	double				getHalfEdgeLen(int iEdge) const;				// get the Euclidean length of the iEdge-th half-edge
	int					getMeshSize() const { return m_nVertex; }
	int					getVerticesNum() const { return m_nVertex; }
	int					getFaceNum() const { return m_nFace; }
	int					getHalfEdgeNum() const { return m_nHalfEdge; }
	int					getFaceNum() { return m_nFace; }
	double				getAvgEdgeLength() const { return m_avgEdgeLen; }
	int					getEdgeNum();		// get number of edges ( not half-edge! )
	int					getBoundaryNum() const;    // get number of boundary edges
	int					getBoundaryVertexNum() const; // get number of boundary vertices
	const Vector3D&		getBoundingBox() const { return m_bBox; }
	const Vector3D&		getCenter() const { return m_Center; }
	const Vector3D&		getVertexPosition(int idx) const { return m_vVertices[idx]->m_vPosition; }

    VectorInt           getOriginalVertexIndex() const;
    VectorInt	        getNeighborVertexIndex(int v, int ring) const;
    VectorInt           getRingVertexIndex(int v, int ring) const;
    VectorInt	        getVertexAdjacentFacesIndex(int vIdx, int ring = 1) const;
    void                getCoordinateFunction(int dim, std::vector<double>& vCoord) const;
    void                setVertexCoordinates(const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord);
    void		        setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<Vector3D>& vNewPos);

    template<typename T> 
    void addAttr(AttrRate rate, const std::string& name) {
        auto iter = mAttributes.find(name);
        if (iter != mAttributes.end()) {
            delete iter->second;
            mAttributes.erase(iter);
        }
        mAttributes.insert(std::make_pair(name, new MeshAttr<T>(rate, name)));        
    }

    template<typename T> 
    void addAttr(const T& data, AttrRate rate, const std::string& name) {
        auto iter = mAttributes.find(name);
        if (iter != mAttributes.end()) {
            delete iter->second;
            mAttributes.erase(iter);
        }
        mAttributes.insert(std::make_pair(name, new MeshAttr<T>(data, rate, name)));        
    }

    template<typename T>
    MeshAttr<T>* getAttr(const std::string& name) {
        auto iter = mAttributes.find(name);
        if (iter != mAttributes.end()) return dynamic_cast<MeshAttr<T>*>(iter->second);
        else return nullptr;
    }

    template<typename T>
    T& getAttrValue(const std::string& name) {
        auto iter = mAttributes.find(name);
        if (iter != mAttributes.end()) 
            return dynamic_cast<MeshAttr<T>*>(iter->second)->getValue();
        else throw std::runtime_error("Requested mesh attribute " + name + " does not exist!");
    }

    /* geometry query and processing */
	void	    gatherStatistics();
	void        gatherStatistics2();
	bool        hasBounary() const;
	int			getEulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int			getMeshGenus();			// get mesh genus
	double		calGaussianCurvatureIntegration();	// compute the integration of Gaussian curvature over all vertices
	bool		calVertexLBO(int i, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, double& Av, std::vector<double>& tw) const;
	bool		calVertexLBO2(int i, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, double& Av, std::vector<double>& tw) const;
	void		calLBO(std::vector<int>& vII, std::vector<int>& vJJ, std::vector<double>& vSS, std::vector<double>& vArea) const;
	bool		calVertexArea(std::vector<double>& Av);
	bool		VertexNeighborGeo(int i, double ring, std::vector<GeoNote>& nbg); // geodesic vertex neighbor
	void		VertexNeighborRing(int i, int ring, std::vector<int>& nbr) const;
	bool		isInNeighborRing(int ref, int query, int ring) const;
	double		calGeodesic(int s, int t) const;
	double		getGeodesicToBoundary(int s) const;	// return 0.0 if in a manifold
	double		getGeodesicToBoundary(int s, std::vector<GeoNote>& nbg);
	double		getVolume();	// calculate volume (area) of a surface
	void		calAreaRatio(CMesh* tmesh, std::vector<int>& ar);	// for registration
	void		calLengthDifference(const CMesh* tmesh, std::vector<double>& ld) const;
	bool		calVertexCurvature( int i );			// calculate number i-th vertex's Gaussian and mean curvature
	void		clearVertexMark();
	void		extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures ) const;
	void        extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary = 1) const;

	static double calAreaMixed(double a, double b, double c, double& cotan_a, double& cotan_c);
	static double calHalfAreaMixed(double a, double b, double c, double& cotan_a);

private:
	void	clearMesh();
	bool	construct();	// construct connectivity
	bool	loadFromOBJ(std::string sFileName);	// load mesh from .obj file
	bool	loadFromM(std::string sFileName);	// load mesh from .m file
	bool	loadFromVERT(std::string sFileName); // load mesh from .vert + .tri files
	bool	loadFromPLY(std::string sFileName);	// load mesh from .ply files
	bool	loadFromOFF(std::string sFileName);
	bool	saveToOBJ(std::string sFileName);	// save mesh to .obj file
	bool    saveToM(const std::string& sFileName );    // save mesh to .m file

	void	calFaceNormalAndArea(int i);			// compute i-th face's normal
	void	calVertexNormal(int i);					// compute i-th vertex's normal
	double  calLocalGeodesic(int ia, int ib, int ic) const;
	void	findHoles();

	void	buildPointerVectors();		//construct vectors of pointers based on array representations
	void	buildIndexArrays();			// already have the pointer-vector representation; fill in the array represenatation
	void	assignElementsIndex();
	bool	isHalfEdgeMergeable(const CHalfEdge* halfEdge);

};	//CMesh

#endif