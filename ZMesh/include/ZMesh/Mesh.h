#ifndef ZMESH_MESH_H
#define ZMESH_MESH_H
#include <cstdio>
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <queue>
#include <string>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <ZGeom/Color.h>
#include "ZGeom/VecN.h"
#include "ZGeom/Vec3.h"
#include "Geometry.h"
#include "Quat.h"
#include "MeshAttr.h"

const int MAX_VERTEX_PER_FACE = 20;
const int MAX_RING_NUMBER = 15;
const int MAX_HOLE_SIZE = 20;

class MeshCoordinates
{
public:
	MeshCoordinates() : mSize(0) {}
	MeshCoordinates(int meshSize) { resize(meshSize); }
	MeshCoordinates(int meshSize, double *cx, double *cy, double *cz) {
		resize(meshSize);
		std::copy_n(cx, meshSize, mCoordX.c_ptr());
		std::copy_n(cy, meshSize, mCoordY.c_ptr());
		std::copy_n(cz, meshSize, mCoordZ.c_ptr());
	}
	int size() const { return mSize; }
	void resize(int n) {
		mSize = n; 
		mCoordX.resize(mSize, 0); 
		mCoordY.resize(mSize, 0); 
		mCoordZ.resize(mSize, 0); 
	}

	void add(double *cx, double *cy, double *cz) {
		mCoordX.add(cx);
		mCoordY.add(cy);
		mCoordZ.add(cz);
	}

	const ZGeom::VecNd& getCoordFunc(int i) const {
		switch (i)
		{
		case 0: return mCoordX;
		case 1: return mCoordY;
		case 2: return mCoordZ;
		default: throw std::logic_error("Invalid mesh coordinate");
		}
	}

	ZGeom::VecNd& getCoordFunc(int i) {
		switch (i) 
		{
		case 0: return mCoordX;
		case 1: return mCoordY;
		case 2: return mCoordZ;
		default: throw std::logic_error("Invalid mesh coordinate");
		}
	}

	ZGeom::VecNd& getXCoord() { return mCoordX; }
	ZGeom::VecNd& getYCoord() { return mCoordY; }
	ZGeom::VecNd& getZCoord() { return mCoordZ; }

	ZGeom::Vec3d getCoordinate(int k) const {
		if (k < 0 || k >= mSize) throw std::logic_error("Vertex index out of bound!");
		return ZGeom::Vec3d(mCoordX[k], mCoordY[k], mCoordZ[k]);
	}
	
	ZGeom::Vec3d operator [] (int k) const { return getCoordinate(k); }

private:
	int mSize;
	ZGeom::VecNd mCoordX, mCoordY, mCoordZ;
};

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

typedef std::priority_queue<GeoNote, std::vector<GeoNote>, std::greater<GeoNote> > GeoQueue;

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
	int					getIndex() const { return m_vIndex; }
	int					getVID() const { return m_vid; }
	std::vector<CFace*> getAdjacentFaces() const;
	CHalfEdge*			getHalfEdge(int ei) { return m_HalfEdges[ei]; }
	const CHalfEdge*    getHalfEdge_const(int ei) const { return m_HalfEdges[ei]; }
	const Vector3D&		getPosition() const { return m_vPosition; } 
	int					getOutValence() const { return (int)m_HalfEdges.size(); }
	bool				judgeOnBoundary();
	bool				isOnBoundary() const { return m_bIsBoundary; }
	bool				isValid() const { return m_bIsValid; }
	void				invalidate(bool flag) { m_bIsValid = flag; }
	void                translateAndScale(const Vector3D& translation, double s);
	void                setPosition( double x, double y, double z );

private:
	void clone(const CVertex& v);

	// ---- fields ---- //
	int						m_vIndex;           // index of the vertex 0-based
	int						m_vid;				// ID of the vertex from original mesh 0-based
	std::vector<CHalfEdge*> m_HalfEdges;		// all half-edges from the vertex
	int*					m_piEdge;			// half edge indices start from this vertex
	int						mOutValence;		// out valence
	Vector3D				m_vPosition;		// vertex coordinates
	bool					m_bIsValid;
	bool					m_bIsBoundary;      // if boundary vertex

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
	CFace*		getAttachedFace() const { return m_Face; }
	CHalfEdge*  getTwinHalfEdge() const { return m_eTwin; }
	bool		isBoundaryEdge() const { return m_eTwin == NULL; }
	int			getVertexIndex(int i) const { return m_iVertex[i]; }
	double		getLength() const;
	bool		isValid() const { return m_bIsValid; }
	int         getIndex() const { return m_eIndex; }
	CVertex*	vert(int i) const { return m_Vertices[i]; }

private:
	// -- fields -- //
	int			m_eIndex;		//half-edge id
	bool		m_bIsValid;

	CVertex*	m_Vertices[2];	//starting and ending vertices
	CHalfEdge*	m_eTwin;		//reverse half-edge; null if boundary half edge
	CHalfEdge*	m_eNext;		//next half-edge (counterclockwise)
	CHalfEdge*	m_ePrev;
	CFace*		m_Face;			//attached face

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
	CFace(int s);
	CFace(const CFace& oldF);
	virtual ~CFace();
	CFace& operator= (const CFace& f);

	// -- operations -- //
	void					Create(int s);
	std::vector<double>		getPlaneFunction();	
	CVertex*				getVertex(int i) const { return m_Vertices[i]; }
	int						getVertexIndex(int i) const { return m_Vertices[i]->getIndex(); }
	double					computeArea() const;
	Vector3D				calcNormal();
	bool					hasVertex(int vidx) const;
	bool					hasHalfEdge() const { return (m_piEdge != NULL); }
	double					distanceToVertex(const CVertex* vq, std::vector<double>& baryCoord);
	int						getFaceIndex() const { return m_fIndex; }

private:
	// ---- fields ---- // 
	int						m_fIndex;
	bool					m_bIsValid;
	int						m_nType;		// number of polygon face edges

	std::vector<CVertex*>	m_Vertices;		//all vertices
	std::vector<CHalfEdge*> m_HalfEdges;	//all half-edges
		
	int*					m_piVertex;		// all vertex index
	int*					m_piEdge;		// all half-edge index
};

typedef MeshAttr< std::vector<ZGeom::Colorf> > AttrVertColor;

//////////////////////////////////////////////////////
//						CMesh   					//
//////////////////////////////////////////////////////
class CMesh 
{
public:
	friend class MeshPyramid;

	/* attribute strings */
	static const std::string StrAttrBoundaryVertCount;
	static const std::string StrAttrBoundaryCount;
	static const std::string StrAttrAvgEdgeLength;
	static const std::string StrAttrMeshCenter;
	static const std::string StrAttrMeshBBox;
	static const std::string StrAttrVertColors;
	static const std::string StrAttrVertGaussCurvatures;
	static const std::string StrAttrVertMeanCurvatures;
	static const std::string StrAttrVertNormal;
	static const std::string StrAttrVertOnHole;
	static const std::string StrAttrVertOnBoundary;
	static const std::string StrAttrFaceNormal;

////////////////   fields    ////////////////
private:
	std::vector<CVertex*>	m_vVertices;
	std::vector<CHalfEdge*> m_vHalfEdges;
	std::vector<CFace*>		m_vFaces;

	bool		m_bIsPointerVectorExist;		// pointer vectors representation
	bool		m_bIsIndexArrayExist;			// index array representation
	bool		m_bSeparateStorage;				// indicate whether the point vectors and index arrays are stored separately 

	int		    m_nVertex;				// number of vertices
	int		    m_nHalfEdge;			// number of half-edges
	int			m_nFace;	 			// number of faces

	CVertex*	m_pVertex;				// array pointer of vertices
	CHalfEdge*	m_pHalfEdge; 			// array pointer of half-edges
	CFace*		m_pFace;				// array pointer of faces

	std::string m_meshName;				// name of the mesh
	std::unordered_map<std::string, MeshAttrBase*> mAttributes;

////////////////    methods    ////////////////
public:
	/* ---- constructors ---- */
	CMesh();
	CMesh(const CMesh& oldMesh);
	virtual ~CMesh();

	/* ---- Mesh IO and processing ---- */
	void        cloneFrom(const CMesh& oldMesh);
	//void		cloneFrom(const CMesh* oldMesh);
	bool		Load(const std::string& sFileName);		// load from file
	bool	    Save(std::string sFileName);			// save to file
	void        move(const Vector3D& translation);		// translate mesh
	void	    scaleAreaToVertexNum();					// move to origin and scale the mesh so that the surface area equals number of vertices
	void        scaleEdgeLenToUnit();					// move to origin and scale the mesh so that the average edge length is 1

	/* ---- attributes access ---- */
	const std::string&	getMeshName() const { return m_meshName; }
	int					vertCount() const { return (int)m_vVertices.size(); }
	int					faceCount() const { return (int)m_vFaces.size(); }
	int					halfEdgeCount() const { return (int)m_vHalfEdges.size(); }
	CVertex*			getVertex(int i) { return m_vVertices[i]; }
	const CVertex*		getVertex(int i) const { return m_vVertices[i]; }
	CFace*				getFace(int i) { return m_vFaces[i]; }
	const CFace*		getFace(int i) const { return m_vFaces[i]; }
	CHalfEdge*			getHalfEdge(int i) { return m_vHalfEdges[i]; }
	const CHalfEdge*	getHalfEdge(int i) const { return m_vHalfEdges[i]; }
	int					getMeshSize() const { return (int)m_vVertices.size(); }
	/*************************************************************************/
	
	/* geometry query, analysis and processing */
	const Vector3D&		         getVertexPosition(int idx) const { return m_vVertices[idx]->m_vPosition; }
	double		     			 getHalfEdgeLen(int iEdge) const;				// get the Euclidean length of the iEdge-th half-edge
	double						 calFaceArea(int i) const;
	const Vector3D&			     getBoundingBox() const { return getAttrValue<Vector3D>(StrAttrMeshBBox); }
	const Vector3D&		         getCenter() const { return getAttrValue<Vector3D>(StrAttrMeshCenter); }
	double						 getAvgEdgeLength() const;
	const std::vector<double>&   getMeanCurvature();
	const std::vector<double>&   getGaussCurvature();
	const std::vector<Vector3D>& getFaceNormals();
	const std::vector<Vector3D>& getVertNormals();
	const std::vector<Vector3D>& getVertNormals_const() const;
	const std::vector<bool>&	 getVertOnHole();
	const std::vector<bool>&     getVertOnHole_const() const;
	const std::vector<bool>&	 getVertOnBoundary();

	void				vertRingNeighborVerts(int vIndex, int ring, std::set<int>& nbr, bool inclusive = false) const;
	void				vertRingNeighborVerts(int i, int ring, std::vector<int>& nbr, bool inclusive = false) const;
	bool				isInNeighborRing(int ref, int query, int ring) const;
	std::vector<int>	getVertNeighborVerts(int v, int ring, bool inclusive = false) const;
	std::vector<int>    getVertIsoNeighborVerts(int v, int ring) const;
	std::vector<int>	getVertexAdjacentFaces(int vIdx, int ring = 1) const;
	bool				vertGeoNeighborVerts(int i, double ring, std::vector<GeoNote>& nbg); // geodesic vertex neighbor

	void                getCoordinateFunction(int dim, std::vector<double>& vCoord) const;
	void				getVertCoordinates(MeshCoordinates& coords) const;
	void				setVertCoordinates(const MeshCoordinates& coords);
	void                setVertexCoordinates(const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord);
	void		        setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<Vector3D>& vNewPos);
	void				gatherStatistics();
	bool				hasBoundary() const;
	int					calBoundaryNum();    // compute number of (connective) boundaries
	int					calBoundaryVert();	 // get number of boundary vertices; set BoundaryVertCount and VertIsOnBoundary attributes
	int					calEulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int					calEdgeCount();		    // get number of edges ( not half-edge! )
	int					calMeshGenus();			// get mesh genus
	double				calGaussianCurvatureIntegration();	// compute the integration of Gaussian curvature over all vertices
	bool				calVertexLBO(int i, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, double& Av, std::vector<double>& tw) const;
	bool				calVertexLBO2(int i, std::vector<int>& Iv, std::vector<int>& Jv, std::vector<double>& Sv, double& Av, std::vector<double>& tw) const;
	void				calLBO(std::vector<int>& vII, std::vector<int>& vJJ, std::vector<double>& vSS, std::vector<double>& vArea) const;
	bool				calVertexArea(std::vector<double>& Av);
	double				calGeodesic(int s, int t) const;
	double				getGeodesicToBoundary(int s) const;	// return 0.0 if in a manifold
	double				getGeodesicToBoundary(int s, std::vector<GeoNote>& nbg);
	double				getVolume();	// calculate volume (area) of a surface
	void				calAreaRatio(CMesh* tmesh, std::vector<int>& ar);	// for registration
	void				calLengthDifference(const CMesh* tmesh, std::vector<double>& ld) const;
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures ) const;
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary = 1) const;
	std::vector<int>    getOriginalVertexIndex() const;
	/*************************************************************************/

	/* MeshAttr functions */
	template<typename T> 
	MeshAttr<T>& addAttr(AttrRate rate, const std::string& name, AttrType attrType = AttrType::UNKNOWN_TYPE) {
		removeAttr(name);
		mAttributes.insert(std::make_pair(name, new MeshAttr<T>(rate, name, attrType)));
		auto iter = mAttributes.find(name);
		return *dynamic_cast<MeshAttr<T>*>(iter->second);        
	}

	template<typename T> 
	void addAttr(const T& data, AttrRate rate, const std::string& name, AttrType attrType = AttrType::UNKNOWN_TYPE) {
		removeAttr(name);
		mAttributes.insert(std::make_pair(name, new MeshAttr<T>(data, rate, name, attrType)));        
	}

	void removeAttr(const std::string& name) {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) {
			delete iter->second;
			mAttributes.erase(iter);
		}
	}
	
	bool hasAttr(const std::string& name) const {
		auto iter = mAttributes.find(name);
		return iter != mAttributes.end();
	}

	template<typename T>
	MeshAttr<T>* getAttr(const std::string& name) {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) return dynamic_cast<MeshAttr<T>*>(iter->second);
		else return nullptr;
	}

	template<typename T>
	const MeshAttr<T>* getAttr(const std::string& name) const {
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

	template<typename T>
	const T& getAttrValue(const std::string& name) const {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) 
			return dynamic_cast<MeshAttr<T>*>(iter->second)->getValue();
		else throw std::runtime_error("Requested mesh attribute " + name + " does not exist!");
	}

	void copyAttributes(const std::unordered_map<std::string, MeshAttrBase*>& attributeMaps) {
		if (&mAttributes == &attributeMaps) return;

		for (auto ma : attributeMaps) {
			MeshAttrBase* a = ma.second->clone();
			mAttributes.insert(std::make_pair(a->getAttrName(), a));
		}
	}

	MeshAttr< std::vector<ZGeom::Colorf> >& addColorAttr(const std::string& colorAttrName)
	{
		return addAttr< std::vector<ZGeom::Colorf> >(AttrRate::VERTEX, colorAttrName, AttrType::CPP_VECTOR_COLOR);
	}

	std::vector<ZGeom::Colorf>& getVertColors(const std::string& colorAttrName)
	{
		return getAttrValue<std::vector<ZGeom::Colorf> >(colorAttrName);
	}

	std::vector<MeshAttr< std::vector<ZGeom::Colorf> > *> getColorAttrLists() {
		std::vector<MeshAttr< std::vector<ZGeom::Colorf> > *> vColorAttr;
		for (auto ap : mAttributes) {
			if (ap.second->attrType() == AttrType::CPP_VECTOR_COLOR) {
				vColorAttr.push_back(dynamic_cast<MeshAttr<std::vector<ZGeom::Colorf> > *>(ap.second));
			}
		}
		return vColorAttr;
	}

	/*************************************************************************/

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

	void	calFaceNormals();			// compute face normals
	void    calVertNormals();			// compute vertex normals
	void	calCurvatures();			// calculate Gaussian and mean curvatures
	void	findHoles();
	double  calLocalGeodesic(int ia, int ib, int ic) const;	

	void	buildPointerVectors();		//construct vectors of pointers based on array representations
	void	buildIndexArrays();			// already have the pointer-vector representation; fill in the array represenatation
	void	assignElementsIndex();
	bool	isHalfEdgeMergeable(const CHalfEdge* halfEdge);
	void	clearVertexMark();
};	//CMesh

#endif