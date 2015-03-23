#ifndef ZGEOM_MESH_H
#define ZGEOM_MESH_H
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <set>
#include <string>
#include <queue>
#include <unordered_map>
#include <vector>
#include "Vec3.h"
#include "VecN.h"
#include "Plane.h"
#include "PointCloud.h"
#include "MeshAttr.h"
#include "MeshCoordinates.h"

class CFace;
class CHalfEdge;
class CMesh;

//////////////////////////////////////////////////////
//						CVertex  					//
//////////////////////////////////////////////////////
class CVertex
{
public:
	CVertex();
	CVertex(double x, double y, double z);
	CVertex(const ZGeom::Vec3d& v);
	CVertex(const CVertex& v);
	CVertex& operator = (const CVertex& v);	

	int					    getIndex() const            { return m_vIndex; }
	CHalfEdge*			    getHalfEdge(int ei)         { return m_HalfEdges[ei]; }
	const CHalfEdge*        getHalfEdge(int ei) const   { return m_HalfEdges[ei]; }    
    const std::vector<CHalfEdge*>& getHalfEdges() const  { return m_HalfEdges;  }
	const ZGeom::Vec3d&	    pos() const                 { return m_vPosition; } 
	int					    outValence() const          { return (int)m_HalfEdges.size(); }
	bool				    isValid() const             { return m_bIsValid; }
    void                    addHalfEdge(CHalfEdge* he)  { m_HalfEdges.push_back(he); }
    void                    removeHalfEdge(CHalfEdge *he);
	void                    init();
    bool				    judgeOnBoundary() const;
    void                    setPosition(const ZGeom::Vec3d& inVec) { m_vPosition = inVec; }
    void                    setPosition(double x, double y, double z);
    void                    clone(const CVertex& v); 
    void                    translateAndScale(const ZGeom::Vec3d& translation, double s);
    std::vector<CFace*>     getAdjacentFaces() const;
    CHalfEdge*              adjacentTo(CVertex* v2) const;

public:
	std::vector<CHalfEdge*> m_HalfEdges;		// all half-edges from the vertex

    int						m_vIndex;           // index of the vertex 0-based
    ZGeom::Vec3d			m_vPosition;		// vertex coordinates
    bool					m_bIsValid;
};

//////////////////////////////////////////////////////
//						CHalfEdge					//
//////////////////////////////////////////////////////
class CHalfEdge
{
public:
	CHalfEdge();
	CHalfEdge(const CHalfEdge& oldE);
	virtual ~CHalfEdge();
	CHalfEdge& operator = (const CHalfEdge& e);

    void                clone(const CHalfEdge& oldEdge);
	CFace*	            getAttachedFace() const     { return m_Face; }
	CHalfEdge*          twinHalfEdge() const        { return m_eTwin; }
	CHalfEdge*          nextHalfEdge() const        { return m_eNext; }
	CHalfEdge*          prevHalfEdge() const;
    CVertex*    	    vert(int i) const           { return (i == 0) ? vert0() : vert1(); }
    CVertex*            vert0() const               { return m_vert; }
    CVertex*            vert1() const               { return nextHalfEdge()->vert0(); }
	bool		        isBoundaryEdge() const      { return m_eTwin == NULL; }
    int			        getVertIndex(int i) const   { return vert(i)->getIndex(); }
	double		        length() const;
	int                 getIndex() const            { return m_eIndex; }	
    void                setVertOrigin(CVertex *v1) { m_vert = v1; }
    ZGeom::Vec3d        midEdge() const { return 0.5 * (vert0()->pos() + vert1()->pos()); }

    static void         makeTwins(CHalfEdge* e1, CHalfEdge* e2);
    static void         makeLoop(CHalfEdge* e1, CHalfEdge* e2, CHalfEdge* e3);

public:
    CVertex*            m_vert;         // starting vert
    CHalfEdge*	        m_eNext;		// next half-edge (counterclockwise)
    CHalfEdge*	        m_eTwin;		// reverse half-edge; null if boundary half edge
	CFace*		        m_Face;			// attached face

    int			        m_eIndex;		// half-edge id
    bool		        m_bIsValid;
};

//////////////////////////////////////////////////////
//						CFace   					//
//////////////////////////////////////////////////////
class CFace
{
public:
	CFace();
	CFace(int s);
	CFace(const CFace& oldF);
	virtual ~CFace();
	CFace& operator= (const CFace& f);
    
    void                    clone(const CFace& f);
    void					create(int s);
    CHalfEdge*              getHalfEdge(int i) { return m_HalfEdges[i]; }
    CHalfEdge*              getHalfEdge(int i) const { return m_HalfEdges[i]; }
    const std::vector<CHalfEdge*>& getAllHalfEdges() const { return m_HalfEdges; }
    std::vector<double>     getAllEdgeLengths() const;
    int                     edgeCount() const { return (int)m_HalfEdges.size(); }
    CVertex*				vert(int i) const { return getHalfEdge(i)->vert(0); }
    std::vector<CVertex*>   getAllVerts() const;
    int						vertIdx(int i) const { return vert(i)->getIndex(); }
    std::vector<int>        getAllVertIdx() const;
    ZGeom::Vec3d            vertPos(int i) const { return vert(i)->pos(); }
    std::vector<ZGeom::Vec3d> getAllVertPos() const;
    bool                    hasVertex(CVertex* pv) const;
	bool					hasVertex(int vidx) const;
	int						getFaceIndex() const { return m_fIndex; }	
	double					calArea() const;
    ZGeom::Vec3d            calBarycenter() const;
    ZGeom::Vec3d            calNormal() const;
    std::vector<double>     getPlaneFunction() const;
    ZGeom::Plane3           getPlane3() const { return ZGeom::Plane3(calNormal(), vertPos(0)); }

public:	
    std::vector<CHalfEdge*> m_HalfEdges;	//all half-edges

    int						m_fIndex;
	bool					m_bIsValid;
	int						m_nType;		// number of polygon face edges
};


//////////////////////////////////////////////////////
//						CMesh   					//
//////////////////////////////////////////////////////
class CMesh 
{
public:
	/* attribute strings */
    static const std::string StrAttrMeshName;    
    static const std::string StrAttrColorSigDefault;
    static const std::string StrAttrNamedCoordinates;
    static const std::string StrAttrCurrentCoordIdx;

    static const std::string StrAttrMeshDescription;
    static const std::string StrAttrVertNormals;
    static const std::string StrAttrFaceNormals;
    static const std::string StrAttrMeshCenter;
    static const std::string StrAttrMeshBBox;
    static const std::string StrAttrAvgEdgeLength;
	static const std::string StrAttrVertOnBoundary;
    static const std::string StrAttrBoundaryVertCount;

////////////////   fields    ////////////////
public:
	std::vector<CVertex*>	m_vVertices;
	std::vector<CHalfEdge*> m_vHalfEdges;
	std::vector<CFace*>		m_vFaces;
    std::unordered_map<std::string, MeshAttrBase*> mAttributes;

////////////////    methods    ////////////////
public:
	/* constructors */
	CMesh();
	CMesh(const CMesh & oldMesh);
    CMesh& operator = (CMesh && oldMesh);
    CMesh(CMesh && oldMesh) { *this = std::move(oldMesh); }
	~CMesh();	

    /* mesh basics */
    void	            clearMesh();
    void                clearAttributes();
    void                clearNonEssentialAttributes();
    void                initAttributes(std::string mesh_name, ZGeom::Colorf default_color);
    void                cloneFrom(const CMesh& oldMesh, const std::string nameSuffix = ".clone");

	/* Mesh IO and processing */
	void		        load(const std::string& sFileName);		// load from file
	void	            save(std::string sFileName);			// save to file
    void                construct(const std::vector<ZGeom::Vec3d>& pVertex, const std::vector<std::vector<int>>& faceVertIdx, int nType = 3);	// construct connectivity
    void	            loadFromOBJ(std::string sFileName);	// load mesh from .obj file
    void	            loadFromM(std::string sFileName);	// load mesh from .m file
    void	            loadFromVERT(std::string sFileName); // load mesh from .vert + .tri files
    void	            loadFromPLY(std::string sFileName);	// load mesh from .ply files
    void	            loadFromOFF(std::string sFileName);
    void	            saveToOBJ(std::string sFileName);	// save mesh to .obj file
    void	            saveToM(const std::string& sFileName);    // save mesh to .m file
    void                getSubMesh(const std::vector<int>& vSubIdx, std::string subMeshName, CMesh& submesh);
    void                getSubMeshFromFaces(const std::vector<int>& vSubFaces, std::string subMeshName, CMesh& submesh);

    /* primitives access */
	int					vertCount() const { return (int)m_vVertices.size(); }
	int					faceCount() const { return (int)m_vFaces.size(); }
	int					halfEdgeCount() const { return (int)m_vHalfEdges.size(); }
	CVertex*			vert(int i) { return m_vVertices[i]; }
	const CVertex*		vert(int i) const { return m_vVertices[i]; }
	CFace*				getFace(int i) { return m_vFaces[i]; }
	const CFace*		getFace(int i) const { return m_vFaces[i]; }
	CHalfEdge*			getHalfEdge(int i) { return m_vHalfEdges[i]; }
	const CHalfEdge*	getHalfEdge(int i) const { return m_vHalfEdges[i]; }
    void	            assignElementsIndex();

    /* basic geometry query, analysis and processing */
	const ZGeom::Vec3d&		        vertPos(int iVert) const { return m_vVertices[iVert]->pos(); }
    std::vector<ZGeom::Vec3d>       allVertPos() const;
	const ZGeom::Vec3d&			    getBoundingBox() const { return getAttrValue<ZGeom::Vec3d>(StrAttrMeshBBox); }
	const ZGeom::Vec3d&		        getCenter() const { return getAttrValue<ZGeom::Vec3d>(StrAttrMeshCenter); }
    double                          getAvgEdgeLength();
	const std::vector<ZGeom::Vec3d>& getFaceNormals();
    bool				            hasBoundary();
    int					            calAttrBoundaryVert();	    // get number of boundary vertices; set BoundaryVertCount and VertIsOnBoundary attributes
	const std::vector<bool>&	    getVertsOnBoundary();
    bool                            isVertOnBoundary(int vi);
	ZGeom::Vec3d			        calMeshCenter() const;
	ZGeom::Vec3d			        calBoundingBox(const ZGeom::Vec3d& center) const;
    double                          calAvgEdgeLength();
	double				            calSurfaceArea() const;
	double				            calVolume() const;
//    std::vector<double>          calPrincipalCurvature(int k);  // k = 0 or 1 or 2

    void                addVertex(CVertex *v);
    void                addHalfEdge(CHalfEdge *e);
    void                addFace(CFace *f);
    CVertex*            faceSplit3(CFace *f);
    CVertex*            faceSplit3(int fIdx);
    CVertex*            edgeSplit(int heIdx);
    CVertex*            edgeSplit(CHalfEdge *he);
    void                edgeSwap(int v1, int v2);
    void                edgeSwap(CHalfEdge* he);
    bool                relaxEdge(CHalfEdge* he);
    static void         makeFace(CHalfEdge* e1, CHalfEdge* e2, CHalfEdge* e3, CFace *f);
    static void         assoicateVertEdges(CVertex *v1, CVertex *v2, CHalfEdge *e12, CHalfEdge *e21);

    void		        scaleAndTranslate(ZGeom::Vec3d translation, double scale);
    void                move(ZGeom::Vec3d translation) { scaleAndTranslate(translation, 1.0); }
    void	            scaleAreaToVertexNum();					// move to origin and scale the mesh so that the surface area equals number of vertices
    void                scaleToUnitBox();                       // move to origin and scale the mesh to inside the unit box
    void                scaleEdgeLenToUnit();					// move to origin and scale the mesh so that the average edge length is 1
    void		        saveToMetis(const std::string& sFileName) const; // save mesh to .mtm Metis-compatible mesh file
    
	bool				isInNeighborRing(int ref, int query, int ring) const;
    std::set<int>       getVertNeighborVertSet(int vIndex, int ring, bool inclusive = false) const;
	std::vector<int>	getVertNeighborVerts(int v, int ring, bool inclusive = false) const;
	std::vector<int>    getVertIsoNeighborVerts(int v, int ring) const;	// get vertices at the distances w.r.t. the given vertex
	std::vector<int>	getVertAdjacentFaceIdx(int vIdx, int ring = 1) const;

	MeshCoordinates     getVertCoordinates() const;
	void				setVertCoordinates(const MeshCoordinates& coords);
	void                setVertCoordinates(const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord);
	void		        setPartialVertCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<ZGeom::Vec3d>& vNewPos);

	int					calEulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int					calEdgeCount();		    // get number of edges ( not half-edge! )
    void	            calAttrFaceNormals();	// compute face normals
    void                calAttrVertNormals();   // compute vert normals
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures );
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary = 1);
    bool	            isHalfEdgeMergeable(const CHalfEdge* halfEdge);
    ZGeom::PointCloud3d toPointCloud() const;
	void                partitionToSubMeshes(const std::vector<std::vector<int>*>& vSubMappedIdx, std::vector<CMesh*>& vSubMeshes) const;


	/************************************************************************/
	/* MeshAttr methods                                                     */
	/************************************************************************/
	bool hasAttr(const std::string& name) const;
    void removeAttr(const std::string& name);
    void copyAttributes(const std::unordered_map<std::string, MeshAttrBase*>& attributeMaps);
    std::vector<std::string> getAttrNamesList() const;

	template<typename T> 
	MeshAttr<T>& addAttr(const std::string& name, AttrRate attrRate, AttrType attrType = AttrType::AT_UNKNOWN) {
		removeAttr(name);
		mAttributes.insert(std::make_pair(name, new MeshAttr<T>(name, attrRate, attrType)));
		auto iter = mAttributes.find(name);
		return *dynamic_cast<MeshAttr<T>*>(iter->second);        
	}

	template<typename T> 
    MeshAttr<T>& addAttr(const T& data, const std::string& name, AttrRate attrRate, AttrType attrType = AttrType::AT_UNKNOWN) {
		removeAttr(name);
		mAttributes.insert(std::make_pair(name, new MeshAttr<T>(data, name, attrRate, attrType)));    
        auto iter = mAttributes.find(name);
        return *dynamic_cast<MeshAttr<T>*>(iter->second);
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
			return dynamic_cast<MeshAttr<T>*>(iter->second)->attrValue();
		else throw std::runtime_error("Requested mesh attribute " + name + " does not exist!");
	}

	template<typename T>
	const T& getAttrValue(const std::string& name) const {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) 
			return dynamic_cast<MeshAttr<T>*>(iter->second)->attrValue();
		else throw std::runtime_error("Requested mesh attribute " + name + " does not exist!");
	}
    	
    AttrType getAttrType(const std::string& name) const {
        auto iter = mAttributes.find(name);
        if (iter == mAttributes.end()) return AT_UNKNOWN;
        else return iter->second->attrType();
    }

	/************************************************************************/
	/* Mesh color attributes methods                                        */
	/************************************************************************/
	AttrVertColors& getColorAttr(const std::string& colorAttrName);
	AttrVertColors& addColorSigAttr(const std::string& colorAttrName);
	void addColorSigAttr(const std::string& colorAttrName, const ZGeom::ColorSignature& vColors);
    void setDefaultColor(ZGeom::Colorf color);
    ZGeom::ColorSignature& getColorSignature(const std::string& colorAttrName);
	std::vector<ZGeom::Colorf>& getVertColors(const std::string& colorAttrName);
	std::vector<AttrVertColors*> getColorAttrList();
        

	/************************************************************************/
	/* Mesh feature attributes methods                                      */
	/************************************************************************/
	AttrMeshFeatures& addAttrMeshFeatures(const std::string& name);
    void addAttrMeshFeatures(const std::vector<int>& featureIdx, const std::string& name);
	void addAttrMeshFeatures(const MeshFeatureList& mfl, const std::string& name);
	const MeshFeatureList& getMeshFeatures(const std::string& name) const;
	std::vector<AttrMeshFeatures*> getMeshFeatureList();
    void addAttrLines(const MeshLineList& vVecs, const std::string& name);
    std::vector<AttrMeshLines*> getMeshLineList();


	/* Vertex scalar attributes methods                                     */
	AttrVertScalars& addAttrVertScalars(const std::string& name);
	void addAttrVertScalars(const std::vector<double>& vScalars, const std::string& name);
	std::vector<double>& getVertScalars(const std::string& name);		

    /************************************************************************/
    /* mesh coordinates attributes methods                                  */
    /************************************************************************/
    void initNamedCoordinates();
    bool hasNamedCoordinates();
    void addNamedCoordinate(const MeshCoordinates& newCoord, const std::string& coordinate_name = "unnamed");
    const std::string& switchNextCoordinate();
    const std::string& switchPrevCoordinate();
    void revertCoordinate();

    /************************************************************************/
    /* mesh string attributes methods */
    /************************************************************************/
    void setMeshName(std::string mesh_name);
    std::string	getMeshName() const;
    void setMeshDescription(std::string descript);
    std::string getMeshDescription() const;
};  // CMesh

#endif
