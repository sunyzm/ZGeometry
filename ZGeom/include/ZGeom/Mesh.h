﻿#ifndef ZMESH_MESH_H
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
#include "VecN.h"
#include "Vec3.h"
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
	MeshCoordinates(int meshSize, double *cx, double *cy, double *cz) 
	{
		resize(meshSize);
		std::copy_n(cx, meshSize, mCoordX.c_ptr());
		std::copy_n(cy, meshSize, mCoordY.c_ptr());
		std::copy_n(cz, meshSize, mCoordZ.c_ptr());
	}

	MeshCoordinates(const MeshCoordinates& mc)
	{
		this->mSize = mc.mSize;
		this->mCoordX = mc.mCoordX;
		this->mCoordY = mc.mCoordY;
		this->mCoordZ = mc.mCoordZ;
	}

	MeshCoordinates& operator = (const MeshCoordinates& mc)
	{
		this->mSize = mc.mSize;
		this->mCoordX = mc.mCoordX;
		this->mCoordY = mc.mCoordY;
		this->mCoordZ = mc.mCoordZ;
		return *this;
	}

	bool empty() const { return mSize == 0; }
	int size() const { return mSize; }

	//resize coordinates size and initialize to 0
	void resize(int n)	
	{
		mSize = n; 
		mCoordX.resize(mSize, 0); 
		mCoordY.resize(mSize, 0); 
		mCoordZ.resize(mSize, 0); 
	}

	void add(double *cx, double *cy, double *cz) 
	{
		mCoordX.add(cx);
		mCoordY.add(cy);
		mCoordZ.add(cz);
	}

	const ZGeom::VecNd& getCoordFunc(int i) const 
	{
		switch (i)
		{
		case 0: return mCoordX;
		case 1: return mCoordY;
		case 2: return mCoordZ;
		default: throw std::logic_error("Invalid mesh coordinate");
		}
	}

	ZGeom::VecNd& getCoordFunc(int i) 
	{
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

	ZGeom::Vec3d getVertCoordinate(int v) const 
	{
		if (v < 0 || v >= mSize) throw std::logic_error("Vertex index out of bound!");
		return ZGeom::Vec3d(mCoordX[v], mCoordY[v], mCoordZ[v]);
	}
	
	ZGeom::Vec3d operator [] (int v) const { return getVertCoordinate(v); }

	double difference(const MeshCoordinates& mc2) const
	{
		assert(this->mSize == mc2.mSize);

		double errorSum(0);
		for (int i = 0; i < mSize; ++i) {
			errorSum += std::pow((this->getVertCoordinate(i) - mc2.getVertCoordinate(i)).length(), 2);
		}
		errorSum = std::sqrt(errorSum);
		return errorSum;
	}

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
	CVertex& operator = (const CVertex& v);	
	virtual ~CVertex();

	// ---- operations ---- //
	int					getIndex() const { return m_vIndex; }
	int					getVID() const { return m_vid; }
	std::vector<const CFace*> getAdjacentFaces() const;
	CHalfEdge*			getHalfEdge(int ei) { return m_HalfEdges[ei]; }
	const CHalfEdge*    getHalfEdge(int ei) const { return m_HalfEdges[ei]; }
	const Vector3D&		getPosition() const { return m_vPosition; } 
	int					outValence() const { return (int)m_HalfEdges.size(); }
	bool				judgeOnBoundary();
	bool				isOnBoundary() const { return m_bIsBoundary; }
	bool				isValid() const { return m_bIsValid; }
	void                translateAndScale(const Vector3D& translation, double s);
	void                setPosition( double x, double y, double z );

private:
	void init();
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
	CHalfEdge& operator = (const CHalfEdge& e);

	// -- operations -- //	
	const CFace*	 getAttachedFace() const { return m_Face; }
	const CHalfEdge* twinHalfEdge() const { return m_eTwin; }
	const CHalfEdge* nextHalfEdge() const { return m_eNext; }
	const CHalfEdge* prevHalfEdge() const { return m_ePrev; }
	const CVertex*	vert(int i) const { return m_Vertices[i]; }

	bool		isBoundaryEdge() const { return m_eTwin == NULL; }
	int			getVertIndex(int i) const { return m_iVertex[i]; }
	double		getLength() const;
	bool		isValid() const { return m_bIsValid; }
	int         getIndex() const { return m_eIndex; }

private:
	void clone(const CHalfEdge& oldEdge);

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
	void					create(int s);
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
	void clone(const CFace& f);

	// ---- fields ---- // 
	int						m_fIndex;
	bool					m_bIsValid;
	int						m_nType;		// number of polygon face edges

	std::vector<CVertex*>	m_Vertices;		//all vertices
	std::vector<CHalfEdge*> m_HalfEdges;	//all half-edges
		
	int*					m_piVertex;		// all vertex index
	int*					m_piEdge;		// all half-edge index
};

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
	~CMesh();

	/* ---- Mesh IO and processing ---- */
	void        cloneFrom(const CMesh& oldMesh, const std::string nameSuffix = ".clone");
	void		load(const std::string& sFileName);		// load from file
	void	    save(std::string sFileName);			// save to file
	void        move(const Vector3D& translation);		// translate mesh
	void	    scaleAreaToVertexNum();					// move to origin and scale the mesh so that the surface area equals number of vertices
	void        scaleEdgeLenToUnit();					// move to origin and scale the mesh so that the average edge length is 1
	void		scaleAndTranslate(const Vector3D& center, double scale);
	void		saveToMetis(const std::string& sFileName) const; // save mesh to .mtm Metis-compatible mesh file
	void		getGraphCSR(std::vector<int>& xadj, std::vector<int>& adjncy) const;
	/* ---- geometry primitives access ---- */
	void				setMeshName(const std::string& meshName) { m_meshName = meshName; }
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
	/*************************************************************************/
	
	/* geometry query, analysis and processing */
	const Vector3D&		         getVertexPosition(int iVert) const { return m_vVertices[iVert]->m_vPosition; }
	double		     			 getHalfEdgeLen(int iEdge) const;				// get the Euclidean length of the iEdge-th half-edge
	double						 calFaceArea(int i) const;
	const Vector3D&			     getBoundingBox() const { return getAttrValue<Vector3D,AR_UNIFORM>(StrAttrMeshBBox); }
	const Vector3D&		         getCenter() const { return getAttrValue<Vector3D,AR_UNIFORM>(StrAttrMeshCenter); }
	double						 getAvgEdgeLength() const;
	const std::vector<double>&   getMeanCurvature();
	const std::vector<double>&   getMeanCurvature() const;
	const std::vector<double>&   getGaussCurvature();
	const std::vector<Vector3D>& getFaceNormals();
	const std::vector<Vector3D>& getVertNormals();
	const std::vector<Vector3D>& getVertNormals() const;
	const std::vector<bool>&	 getVertsOnHole();
	const std::vector<bool>&     getVertsOnHole_const() const;
	const std::vector<bool>&	 getVertsOnBoundary();
	Vector3D			         calMeshCenter() const;
	Vector3D			         calBoundingBox(const Vector3D& center) const;
	double				         calSurfaceArea() const;
	double				         calVolume() const;

	void				vertRingNeighborVerts(int vIndex, int ring, std::set<int>& nbr, bool inclusive = false) const;
	void				vertRingNeighborVerts(int i, int ring, std::vector<int>& nbr, bool inclusive = false) const;
	bool				isInNeighborRing(int ref, int query, int ring) const;
	std::vector<int>	getVertNeighborVerts(int v, int ring, bool inclusive = false) const;
	std::vector<int>    getVertIsoNeighborVerts(int v, int ring) const;	// get vertices at the distances w.r.t. the given vertex
	std::vector<int>	getVertexAdjacentFaceIdx(int vIdx, int ring = 1) const;
	bool				vertGeoNeighborVerts(int i, double ring, std::vector<GeoNote>& nbg); // geodesic vertex neighbor

	std::vector<int>    getOriginalVertexIndex() const;
	void                getVertCoordinateFunction(int dim, std::vector<double>& vCoord) const;
	void				getVertCoordinates(MeshCoordinates& coords) const;
	void				setVertCoordinates(const MeshCoordinates& coords);
	void                setVertexCoordinates(const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord);
	void		        setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<Vector3D>& vNewPos);
	void				diffCoordinates(const MeshCoordinates& coordsToCompare, std::vector<double>& vDiff) const;

	void				gatherStatistics();
	bool				hasBoundary() const;
	int					calBoundaryNum();    // compute number of (connective) boundaries
	int					calBoundaryVert();	 // get number of boundary vertices; set BoundaryVertCount and VertIsOnBoundary attributes
	int					calEulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int					calEdgeCount();		    // get number of edges ( not half-edge! )
	int					calMeshGenus();			// get mesh genus
	double				calGaussianCurvatureIntegration();	// compute the integration of Gaussian curvature over all vertices
	bool				calVertexArea(std::vector<double>& Av);
	double				calGeodesic(int s, int t) const;
	double				getGeodesicToBoundary(int s) const;	// return 0.0 if in a manifold
	double				getGeodesicToBoundary(int s, std::vector<GeoNote>& nbg);
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures ) const;
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary = 1) const;

	void                partitionToSubMeshes(const std::vector<std::vector<int>*>& vSubMappedIdx, std::vector<CMesh*>& vSubMeshes) const;
	/*************************************************************************/

	/************************************************************************/
	/* MeshAttr methods                                                     */
	/************************************************************************/
	bool hasAttr(const std::string& name) const {
		auto iter = mAttributes.find(name);
		return iter != mAttributes.end();
	}

	template<typename T, AttrRate R> 
	MeshAttr<T,R>& addAttr(const std::string& name, AttrType attrType = AttrType::AT_UNKNOWN) {
		removeAttr(name);
		mAttributes.insert(std::make_pair(name, new MeshAttr<T,R>(name, attrType)));
		auto iter = mAttributes.find(name);
		return *dynamic_cast<MeshAttr<T,R>*>(iter->second);        
	}

	template<typename T, AttrRate R> 
	void addAttr(const T& data, const std::string& name, AttrType attrType = AttrType::AT_UNKNOWN) {
		removeAttr(name);
		mAttributes.insert(std::make_pair(name, new MeshAttr<T,R>(data, name, attrType)));        
	}

	void removeAttr(const std::string& name) {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) {
			delete iter->second;
			mAttributes.erase(iter);
		}
	}
		
	template<typename T, AttrRate R>
	MeshAttr<T,R>* getAttr(const std::string& name) {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) return dynamic_cast<MeshAttr<T,R>*>(iter->second);
		else return nullptr;
	}
	
	template<typename T, AttrRate R>
	const MeshAttr<T,R>* getAttr(const std::string& name) const {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) return dynamic_cast<MeshAttr<T,R>*>(iter->second);
		else return nullptr;
	}

	template<typename T, AttrRate R>
	T& getAttrValue(const std::string& name) {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) 
			return dynamic_cast<MeshAttr<T,R>*>(iter->second)->attrValue();
		else throw std::runtime_error("Requested mesh attribute " + name + " does not exist!");
	}

	template<typename T, AttrRate R>
	const T& getAttrValue(const std::string& name) const {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) 
			return dynamic_cast<MeshAttr<T,R>*>(iter->second)->attrValue();
		else throw std::runtime_error("Requested mesh attribute " + name + " does not exist!");
	}

	void copyAttributes(const std::unordered_map<std::string, MeshAttrBase*>& attributeMaps) {
		if (&mAttributes == &attributeMaps) return;

		for (auto ma : attributeMaps) {
			MeshAttrBase* a = ma.second->clone();
			mAttributes.insert(std::make_pair(a->attrName(), a));
		}
	}

	std::vector<std::string> getAttrNamesList() {
		std::vector<std::string> vAttrNames;
		for (auto ap : mAttributes) vAttrNames.push_back(ap.second->attrName());
		std::sort(vAttrNames.begin(), vAttrNames.end());
		return vAttrNames;
	}

	/************************************************************************/
	/* Color attributes methods                                           */
	/************************************************************************/
	AttrVertColors& getColorAttr(const std::string& colorAttrName) {
		return *getAttr<std::vector<ZGeom::Colorf>, AR_VERTEX>(colorAttrName);
	}

	AttrVertColors& addColorAttr(const std::string& colorAttrName) {
		if (hasAttr(colorAttrName)) return getColorAttr(colorAttrName);
		else return addAttr<std::vector<ZGeom::Colorf>, AR_VERTEX>(colorAttrName, AttrType::AT_VEC_COLOR);
	}

	std::vector<ZGeom::Colorf>& getVertColors(const std::string& colorAttrName) {
		return getAttrValue<std::vector<ZGeom::Colorf>, AR_VERTEX>(colorAttrName);
	}

	std::vector<AttrVertColors*> getColorAttrList() {
		std::vector<AttrVertColors*> vColorAttr;
		for (auto ap : mAttributes) {
			if (ap.second->attrType() == AttrType::AT_VEC_COLOR) {
				vColorAttr.push_back(dynamic_cast<AttrVertColors*>(ap.second));
			}
		}
		return vColorAttr;
	}

	/************************************************************************/
	/* Mesh feature attributes methods                                      */
	/************************************************************************/
	AttrMeshFeatures& addAttrMeshFeatures(const std::string& name) {
		return addAttr<MeshFeatureList, AR_UNIFORM>(name, AT_FEATURES);
	}

	void addAttrMeshFeatures(const MeshFeatureList& mfl, const std::string& name) {
		addAttr<MeshFeatureList, AR_UNIFORM>(mfl, name, AT_FEATURES);
	}

	const MeshFeatureList& getMeshFeatures(const std::string& name) const {
		return getAttrValue<MeshFeatureList, AR_UNIFORM>(name);
	}

	std::vector<AttrMeshFeatures*> getMeshFeatureList() {
		std::vector<AttrMeshFeatures*> vMeshFeatures;
		for (auto ap : mAttributes) {
			if (ap.second->attrType() == AttrType::AT_FEATURES) {
				vMeshFeatures.push_back(dynamic_cast<AttrMeshFeatures*>(ap.second));
			}
		}
		return vMeshFeatures;
	}

	/************************************************************************/
	/* Vertex scalar attributes methods                                     */
	/************************************************************************/
	AttrVertScalars& addAttrVertScalars(const std::string& name) {
		return addAttr<std::vector<double>, AR_VERTEX>(name, AT_VEC_DBL);
	}

	void addAttrVertScalars(const std::vector<double>& vScalars, const std::string& name) {
		assert(vScalars.size() == vertCount());
		addAttr<std::vector<double>, AR_VERTEX>(vScalars, name, AT_VEC_DBL);
	}

	std::vector<double>& getVertScalars(const std::string& name) {
		return getAttrValue<std::vector<double>, AR_VERTEX>(name);
	}
		
	//////////////////////////////////////////////////////////////////////////	

private:
	void	clearMesh();
	void	construct();	// construct connectivity
	void	loadFromOBJ(std::string sFileName);	// load mesh from .obj file
	void	loadFromM(std::string sFileName);	// load mesh from .m file
	void	loadFromVERT(std::string sFileName); // load mesh from .vert + .tri files
	void	loadFromPLY(std::string sFileName);	// load mesh from .ply files
	void	loadFromOFF(std::string sFileName);
	void	saveToOBJ(std::string sFileName);	// save mesh to .obj file
	void	saveToM(const std::string& sFileName );    // save mesh to .m file

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
	
	/* helper functions */
	static double calAreaMixed(double a, double b, double c, double& cotan_a, double& cotan_c);
	static double calHalfAreaMixed(double a, double b, double c, double& cotan_a);
}; //CMesh

#endif
