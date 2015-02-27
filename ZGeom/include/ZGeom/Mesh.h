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

const int MAX_VERTEX_PER_FACE = 20;
const int MAX_RING_NUMBER = 15;
const int MAX_HOLE_SIZE = 200;

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
	int					    getVID() const              { return m_vid; }
	CHalfEdge*			    getHalfEdge(int ei)         { return m_HalfEdges[ei]; }
	const CHalfEdge*        getHalfEdge(int ei) const   { return m_HalfEdges[ei]; }    
    const std::vector<CHalfEdge*> getHalfEdges() const  { return m_HalfEdges;  }
	ZGeom::Vec3d		    pos() const                 { return (ZGeom::Vec3d)m_vPosition; } 
	int					    outValence() const          { return (int)m_HalfEdges.size(); }
	bool				    isValid() const             { return m_bIsValid; }
    void                    addHalfEdge(CHalfEdge* he)  { m_HalfEdges.push_back(he); }
    void                    removeHalfEdge(CHalfEdge *he);
	void                    init();
    bool				    judgeOnBoundary() const;
    void                    setPosition(double x, double y, double z);
    void                    clone(const CVertex& v); 
    void                    translateAndScale(const ZGeom::Vec3d& translation, double s);
    std::vector<const CFace*>   getAdjacentFaces() const;
    CHalfEdge*              adjacentTo(CVertex* v2) const;

public:
	int						m_vIndex;           // index of the vertex 0-based
	int						m_vid;				// ID of the vertex from original mesh 0-based
	std::vector<CHalfEdge*> m_HalfEdges;		// all half-edges from the vertex
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
	CHalfEdge*          prevHalfEdge() const        { return m_ePrev; }
	CVertex*    	    vert(int i) const           { return m_Vertices[i]; }
	bool		        isBoundaryEdge() const      { return m_eTwin == NULL; }
    int			        getVertIndex(int i) const   { return vert(i)->getIndex(); }
	double		        length() const;
	int                 getIndex() const            { return m_eIndex; }	
    void                setVerts(CVertex* v1, CVertex* v2) { m_Vertices[0] = v1; m_Vertices[1] = v2; }

    static void         makeTwins(CHalfEdge* e1, CHalfEdge* e2);
    static void         makeLoop(CHalfEdge* e1, CHalfEdge* e2, CHalfEdge* e3);

public:
	int			        m_eIndex;		// half-edge id
	CVertex*	        m_Vertices[2];	// starting and ending vertices
	CHalfEdge*	        m_eTwin;		// reverse half-edge; null if boundary half edge
	CHalfEdge*	        m_eNext;		// next half-edge (counterclockwise)
	CHalfEdge*	        m_ePrev;
	CFace*		        m_Face;			// attached face
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
	CVertex*				getVertex(int i) const { return m_Vertices[i]; }
	int						getVertexIndex(int i) const { return m_Vertices[i]->getIndex(); }
    int                     edgeCount() const { return m_nType; }
    bool                    hasVertex(CVertex* pv) const;
	bool					hasVertex(int vidx) const;
	double					distanceToVertex(const CVertex* vq, std::vector<double>& baryCoord);
	int						getFaceIndex() const { return m_fIndex; }	
	double					calArea() const;
    ZGeom::Vec3d            calBarycenter() const;
    ZGeom::Vec3d            calcNormal() const;
    std::vector<double>     getPlaneFunction() const;
    ZGeom::Plane3           getPlane3() const { return ZGeom::Plane3(calcNormal(), m_Vertices[0]->pos()); }
    std::vector<CVertex*>   getAllVerts() const { return m_Vertices; }
    std::vector<int>        getAllVertIdx() const;
    std::vector<ZGeom::Vec3d> getAllVertCoords() const;

public:	
	int						m_fIndex;
	bool					m_bIsValid;
	int						m_nType;		// number of polygon face edges
	std::vector<CVertex*>	m_Vertices;		//all vertices
	std::vector<CHalfEdge*> m_HalfEdges;	//all half-edges
};


//////////////////////////////////////////////////////
//						CMesh   					//
//////////////////////////////////////////////////////
class CMesh 
{
public:
	/* attribute strings */
    static const std::string StrAttrAvgEdgeLength;
    static const std::string StrAttrMeshCenter;
    static const std::string StrAttrMeshBBox;
	static const std::string StrAttrBoundaryVertCount;
    static const std::string StrAttrBoundaryLoops;
	static const std::string StrAttrVertColors;
    static const std::string StrAttrColorDefault;
	static const std::string StrAttrVertNormal;
	static const std::string StrAttrVertOnHole;
	static const std::string StrAttrVertOnBoundary;
	static const std::string StrAttrFaceNormal;

////////////////   fields    ////////////////
public:
	std::vector<CVertex*>	m_vVertices;
	std::vector<CHalfEdge*> m_vHalfEdges;
	std::vector<CFace*>		m_vFaces;
    std::string             m_meshName;				// name of the mesh
    std::unordered_map<std::string, MeshAttrBase*> mAttributes;
    bool                    m_verbose;
    ZGeom::Colorf           m_defaultColor;

////////////////    methods    ////////////////
public:
	/* constructors */
	CMesh();
	CMesh(const CMesh& oldMesh);
	~CMesh();	

    /* mesh basics */
    void	            clearMesh();
    void                clearAttributes();
    void                cloneFrom(const CMesh& oldMesh, const std::string nameSuffix = ".clone");
    void                setVerbose(bool verbose) { m_verbose = verbose; }
    void				setMeshName(const std::string& meshName) { m_meshName = meshName; }
    std::string	        getMeshName() const { return m_meshName; }

	/* Mesh IO and processing */
	void		        load(const std::string& sFileName);		// load from file
	void	            save(std::string sFileName);			// save to file
    void                construct(const std::vector<CVertex>& pVertex, const std::vector<std::vector<int>>& faceVerts, int nType = 3);	// construct connectivity
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
	const ZGeom::Vec3d&		     vertPos(int iVert) const { return m_vVertices[iVert]->pos(); }
    std::vector<ZGeom::Vec3d>    allVertPos() const;
	const ZGeom::Vec3d&			 getBoundingBox() const { return getAttrValue<ZGeom::Vec3d>(StrAttrMeshBBox); }
	const ZGeom::Vec3d&		     getCenter() const { return getAttrValue<ZGeom::Vec3d>(StrAttrMeshCenter); }
	double						 getAvgEdgeLength() const;
	const std::vector<ZGeom::Vec3d>& getFaceNormals();
	const std::vector<ZGeom::Vec3d>& getVertNormals() const;
	const std::vector<bool>&	 getVertsOnBoundary();
    std::vector<bool>            getVertsOnHoles();
    bool                         isVertOnBoundary(int vi);
	ZGeom::Vec3d			     calMeshCenter() const;
	ZGeom::Vec3d			     calBoundingBox(const ZGeom::Vec3d& center) const;
	double				         calSurfaceArea() const;
	double				         calVolume() const;
//    std::vector<double>          calPrincipalCurvature(int k);  // k = 0 or 1 or 2

    void                addVertex(CVertex *v);
    void                addHalfEdge(CHalfEdge *e);
    void                addFace(CFace *f);
    void                faceSplit(int fIdx);
    CVertex*            faceSplit(CFace* face);
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
	std::vector<int>	getVertNeighborVerts(int v, int ring, bool inclusive = false) const;
	std::vector<int>    getVertIsoNeighborVerts(int v, int ring) const;	// get vertices at the distances w.r.t. the given vertex
	std::vector<int>	getVertexAdjacentFaceIdx(int vIdx, int ring = 1) const;

	MeshCoordinates     getVertCoordinates() const;
	void				setVertCoordinates(const MeshCoordinates& coords);
	void                setVertexCoordinates(const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord);
	void		        setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<ZGeom::Vec3d>& vNewPos);
    void				vertRingNeighborVerts(int vIndex, int ring, std::set<int>& nbr, bool inclusive = false) const;
    void				vertRingNeighborVerts(int i, int ring, std::vector<int>& nbr, bool inclusive = false) const;

	bool				hasBoundary();
	int					calAttrBoundaryLoops();     // compute number of (connective) boundaries
    std::vector<std::vector<int>> getBoundaryLoopVerts(); // get the vertex indices of each boundary loop
    std::vector<std::vector<int>> getBoundaryLoopEdges();   // get the halfedge indices of each boundary loop
	int					calAttrBoundaryVert();	    // get number of boundary vertices; set BoundaryVertCount and VertIsOnBoundary attributes
	int					calEulerNum();			// get Euler number of mesh: Euler# = v - e + f
	int					calEdgeCount();		    // get number of edges ( not half-edge! )
	int					calMeshGenus();			// get mesh genus
    void	            calAttrFaceNormals();			// compute face normals
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, std::vector<int>& vFeatures );
	void				extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary = 1);
    bool	            isHalfEdgeMergeable(const CHalfEdge* halfEdge);
    ZGeom::PointCloud3d toPointCloud() const;
	void                partitionToSubMeshes(const std::vector<std::vector<int>*>& vSubMappedIdx, std::vector<CMesh*>& vSubMeshes) const;

    /*************************************************************************/

	/************************************************************************/
	/* MeshAttr methods                                                     */
	/************************************************************************/
	bool hasAttr(const std::string& name) const {
		auto iter = mAttributes.find(name);
		return iter != mAttributes.end();
	}

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

	void removeAttr(const std::string& name) {
		auto iter = mAttributes.find(name);
		if (iter != mAttributes.end()) {
			delete iter->second;
			mAttributes.erase(iter);
		}
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
	/* Mesh color attributes methods                                        */
	/************************************************************************/
	AttrVertColors& getColorAttr(const std::string& colorAttrName) {
		return *getAttr<ZGeom::ColorSignature>(colorAttrName);
	}

	AttrVertColors& addColorAttr(const std::string& colorAttrName) {
		if (hasAttr(colorAttrName)) return getColorAttr(colorAttrName);
		else return addAttr<ZGeom::ColorSignature>(colorAttrName, AttrRate::AR_VERTEX, AttrType::AT_VEC_COLOR);
	}

	void addColorAttr(const std::string& colorAttrName, const ZGeom::ColorSignature& vColors) {
		if (hasAttr(colorAttrName)) getColorAttr(colorAttrName).attrValue() = vColors;
		else {
            ZGeom::ColorSignature& vNewColor = addColorAttr(colorAttrName).attrValue();
			vNewColor = vColors;
		}
	}

    void setDefaultColor(ZGeom::Colorf color);
    void addDefaultColorAttr();

    ZGeom::ColorSignature& getColorSignature(const std::string& colorAttrName) {
        return getAttrValue<ZGeom::ColorSignature>(colorAttrName);
    }

	std::vector<ZGeom::Colorf>& getVertColors(const std::string& colorAttrName) {
        return getAttrValue<ZGeom::ColorSignature>(colorAttrName).getColors();
	}

	std::vector<AttrVertColors*> getColorAttrList() {
		std::vector<AttrVertColors*> vColorAttr;
		for (auto ap : mAttributes) {
			if (ap.second->attrType() == AttrType::AT_VEC_COLOR && ap.second->attrRate() == AttrRate::AR_VERTEX) {
				vColorAttr.push_back(dynamic_cast<AttrVertColors*>(ap.second));
			}
		}
		return vColorAttr;
	}
        

	/************************************************************************/
	/* Mesh feature attributes methods                                      */
	/************************************************************************/
	AttrMeshFeatures& addAttrMeshFeatures(const std::string& name) {
		return addAttr<MeshFeatureList>(name, AR_UNIFORM, AT_FEATURES);
	}

    void addAttrMeshFeatures(const std::vector<int>& featureIdx, const std::string& name);

	void addAttrMeshFeatures(const MeshFeatureList& mfl, const std::string& name) {
		addAttr<MeshFeatureList>(mfl, name, AR_UNIFORM, AT_FEATURES);
	}

	const MeshFeatureList& getMeshFeatures(const std::string& name) const {
		return getAttrValue<MeshFeatureList>(name);
	}

	std::vector<AttrMeshFeatures*> getMeshFeatureList() {
		std::vector<AttrMeshFeatures*> vMeshFeatures;
		for (auto ap : mAttributes) {
            if (ap.second->attrType() == AttrType::AT_FEATURES && ap.second->attrRate() == AttrRate::AR_UNIFORM) {
				vMeshFeatures.push_back(dynamic_cast<AttrMeshFeatures*>(ap.second));
			}
		}
		return vMeshFeatures;
	}

    void addAttrLines(const MeshLineList& vVecs, const std::string& name) {
        addAttr<MeshLineList>(vVecs, name, AR_UNIFORM, AT_VEC_LINE);
    }

    std::vector<AttrMeshLines*> getMeshLineList() {
        std::vector<AttrMeshLines*> vMeshLines;
        for (auto ap : mAttributes) {
            if (ap.second->attrType() == AttrType::AT_VEC_LINE)
                vMeshLines.push_back(dynamic_cast<AttrMeshLines*>(ap.second));
        }
        return vMeshLines;
    }

	/************************************************************************/
	/* Vertex scalar attributes methods                                     */
	/************************************************************************/
	AttrVertScalars& addAttrVertScalars(const std::string& name) {
		return addAttr<std::vector<double>>(name, AR_VERTEX, AT_VEC_DBL);
	}

	void addAttrVertScalars(const std::vector<double>& vScalars, const std::string& name) {
		assert(vScalars.size() == vertCount());
		addAttr<std::vector<double>>(vScalars, name, AR_VERTEX, AT_VEC_DBL);
	}

	std::vector<double>& getVertScalars(const std::string& name) {
		return getAttrValue<std::vector<double>>(name);
	}		
	//////////////////////////////////////////////////////////////////////////	

};  // CMesh

#endif
