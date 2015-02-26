#include "Mesh.h"
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <set>
#include <map>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "zassert.h"
#include "arithmetic.h"

using namespace std;
using ZGeom::Colorf;
using ZGeom::Vec3d;

/* basic Mesh attribute names */
const std::string CMesh::StrAttrAvgEdgeLength			= "mesh_average_edge_length";
const std::string CMesh::StrAttrMeshBBox				= "mesh_bounding_box";
const std::string CMesh::StrAttrMeshCenter				= "mesh_center";
const std::string CMesh::StrAttrVertColors				= "vert_color";
const std::string CMesh::StrAttrColorDefault            = "vert_color_default";
const std::string CMesh::StrAttrVertPrincipalCurvatures1 = "vert_principal_curvature_1";
const std::string CMesh::StrAttrVertPrincipalCurvatures2 = "vert_principal_curvature_2";
const std::string CMesh::StrAttrVertNormal				= "vert_normal";
const std::string CMesh::StrAttrFaceNormal				= "face_normal";
const std::string CMesh::StrAttrBoundaryVertCount       = "mesh_boundary_vert_count";
const std::string CMesh::StrAttrBoundaryLoops           = "mesh_boundary_loops";
const std::string CMesh::StrAttrVertOnHole              = "vert_on_hole";
const std::string CMesh::StrAttrVertOnBoundary          = "vert_on_boundary";

//////////////////////////////////////////////////////
//						CVertex						//
//////////////////////////////////////////////////////

CVertex::CVertex()
{
	init();
}

CVertex::CVertex( double x, double y, double z )
{
	init();
	m_vPosition = ZGeom::Vec3d(x,y,z);
}

CVertex::CVertex( const ZGeom::Vec3d& v )
{
	init();
	m_vPosition = v;
}

CVertex::CVertex( const CVertex& v )
{	
	clone(v);
}

CVertex& CVertex::operator = ( const CVertex& v )
{	
	clone(v);
	return *this;
}

void CVertex::init()
{
	m_bIsValid = true;
}

void CVertex::clone(const CVertex& v)
{
	if (this == &v) return;
	m_vIndex			= v.m_vIndex;
	m_vid				= v.m_vid;
	m_vPosition			= v.m_vPosition;
	m_bIsValid			= v.m_bIsValid;

    m_HalfEdges.clear();
}

std::vector<const CFace*> CVertex::getAdjacentFaces() const
{
	vector<const CFace*> pFaces;
	for (CHalfEdge* he : m_HalfEdges) pFaces.push_back(he->getAttachedFace());
	
	return pFaces;
}

bool CVertex::judgeOnBoundary() const
{
	for (auto he : m_HalfEdges) {
		if (he->twinHalfEdge() == NULL) return true;
	}	
	return false;
}

void CVertex::translateAndScale( const ZGeom::Vec3d& translation, double s )
{
	m_vPosition += translation;
	m_vPosition *= s;
}

void CVertex::setPosition( double x, double y, double z )
{
	m_vPosition.x = x;
	m_vPosition.y = y;
	m_vPosition.z = z;
}

CHalfEdge * CVertex::adjacentTo(CVertex* v2) const
{
    for (CHalfEdge* he : m_HalfEdges) {
        if (he->vert(1) == v2) return he;
    }
    return NULL;
}

void CVertex::removeHalfEdge(CHalfEdge *he)
{
    std::remove(m_HalfEdges.begin(), m_HalfEdges.end(), he);
}

//////////////////////////////////////////////////////
//						CHalfEdge					//
//////////////////////////////////////////////////////

CHalfEdge::CHalfEdge()
{
	m_Vertices[0] = m_Vertices[1] = NULL;
	m_eTwin= m_eNext = m_ePrev = NULL;
	m_Face = NULL;
	m_bIsValid = true;
}

CHalfEdge::CHalfEdge( const CHalfEdge& e )
{
	clone(e);
}

CHalfEdge& CHalfEdge::operator = (const CHalfEdge& e)
{
	clone(e);
	return *this;
}

CHalfEdge::~CHalfEdge()
{
}

void CHalfEdge::clone(const CHalfEdge& e)
{
	if (this == &e) return;
	m_eIndex = e.m_eIndex;
    m_bIsValid = e.m_bIsValid;

	m_Vertices[0] = m_Vertices[1] = NULL;
	m_eTwin = m_eNext = m_ePrev = NULL;
	m_Face = NULL;	
}

double CHalfEdge::length() const
{
	return (m_Vertices[0]->pos() - m_Vertices[1]->pos()).length();
}

void CHalfEdge::makeTwins(CHalfEdge* e1, CHalfEdge* e2)
{
    if (e1 == NULL || e2 == NULL) return;
    e1->m_eTwin = e2;
    e2->m_eTwin = e1;
}

void CHalfEdge::makeLoop(CHalfEdge* e1, CHalfEdge* e2, CHalfEdge* e3)
{
    e1->m_eNext = e3->m_ePrev = e2;
    e2->m_eNext = e1->m_ePrev = e3;
    e3->m_eNext = e2->m_ePrev = e1;
}



//////////////////////////////////////////////////////
//						CFace						//
//////////////////////////////////////////////////////

CFace::CFace() : m_nType(0), m_bIsValid(true)
{	
}

CFace::CFace(int s) : m_nType(s), m_bIsValid(true)
{	
}

CFace::CFace( const CFace& f )
{
	clone(f);
}

CFace& CFace::operator = (const CFace& f)
{
	clone(f);
	return *this;
}

CFace::~CFace()
{
}

void CFace::clone( const CFace& f )
{
	if (this == &f) return;
	m_fIndex	= f.m_fIndex;
	m_nType		= f.m_nType;
    m_bIsValid = f.m_bIsValid;
}

void CFace::create(int s)
{
	m_nType = s;
}

bool CFace::hasVertex( int vidx ) const
{
	for (CVertex* pv : m_Vertices) {
		if (pv->getIndex() == vidx)	return true;
	}

	return false;
}

bool CFace::hasVertex(CVertex* pv) const
{
    return std::find(m_Vertices.begin(), m_Vertices.end(), pv) != m_Vertices.end();
}

double CFace::distanceToVertex( const CVertex* vp, std::vector<double>& baryCoord )
{	
	//baryCoord.resize(3, 0);
	
	/**** adapted from WildMagic ****/
	ZGeom::Vec3d V[3] = {m_Vertices[0]->pos(), m_Vertices[1]->pos(), m_Vertices[2]->pos()};
	ZGeom::Vec3d p = vp->pos();

	ZGeom::Vec3d diff = V[0] - p;
	ZGeom::Vec3d edge0 = V[1] - V[0];
	ZGeom::Vec3d edge1 = V[2] - V[0];
	double   a = edge0.length2();
	double   b = edge0.dot(edge1);
	double   c = edge1.length2();
	double   d = diff.dot(edge0);
	double   e = diff.dot(edge1);
	double   f = diff.length2();
	double det = std::abs(a*c - b*b);
	double   s = b*e - c*d;
	double   t = b*d - a*e;
	//double s_bar = s / det, t_bar = t / det;
	double sqrDistance;
	
	if (s + t <= det) 
	{
		if (s < 0.) {
			if (t < 0.) { // region 4
				if (d < 0.) {
					if (-d >= a) {
						sqrDistance = a + 2.*d + f;	// on V1
						s = 1; t = 0;
					} else {
						sqrDistance = f - d*d/a; // on E0
						s = -d/a; t = 0;
					}
				} 
				else {
					if (e >= 0.) {
						sqrDistance = f;   // on V0
						s = 0; t = 0;
					} else if (-e >= c) {
						sqrDistance = c + 2.*e + f;	// on V2
						s = 0; t = 1;
					} else {
						sqrDistance = f - e*e/c;	//on E1
						s = 0; t = -e/c;
					}
				}
			}
			else {  // region 3
				if (e >= 0.) {
					sqrDistance = f;	// on V0
					s = 0; t = 0;
				} else if (-e >= c) {
					sqrDistance = c + 2.*e + f;	// on V2
					s = 0; t = 1;
				} else {
					sqrDistance = f - e*e/c;	//on E1
					s = 0; t = -e/c;
				}
			}
		} 
		else if (t < 0.)  { // region 5
			if (d >= 0.) {
				sqrDistance = f;	// on V0
				s = 0; t = 0;
			} else if (-d >= a) {
				sqrDistance = a + 2.*d + f;	// on V1
				s = 1; t = 0;
			} else {
				sqrDistance = d*s + f - d*d/a;	// on E0
				s = -d/a; t = 0;
			}
		}
		else  { // region 0
			// The minimum is at an interior point of the triangle.
			double invDet = 1./det;
			s *= invDet;
			t *= invDet;
			sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
		}
	} // if (s + t <= det)
	else
	{
		double tmp0, tmp1, numer, denom;

		if (s < 0.)  {// region 2
			tmp0 = b + d;
			tmp1 = c + e;
			if (tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a - 2.*b + c;
				if (numer >= denom) {
					sqrDistance = a + 2.*d + f;	// on V1?
					s = 1; t = 0;
				} else {
					s = numer/denom;
					t = 1. - s;
					sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
				}
			}
			else {
				if (tmp1 <= 0.) {
					sqrDistance = c + 2.*e + f;	//on v2
					s = 0; t = 1;
				} else if (e >= 0.) {
					sqrDistance = f;	// on v0
					s = 0; t = 0;
				} else {
					sqrDistance = f - e*e/c;	// on E1?
					s = 0; t = -e/c;
				}
			}
		}
		else if (t < 0.) { // region 6
			tmp0 = b + e;
			tmp1 = a + d;
			if (tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a - 2.*b + c;
				if (numer >= denom) {
					sqrDistance = c + 2.*e + f;	// on V2
					s = 0.; t = 1.; 
				} else {
					t = numer/denom;
					s = 1. - t;
					sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
				}
			}
			else {
				if (tmp1 <= 0.) {
					sqrDistance = a + 2.*d + f;	// on V1
					s = 1.; t = 0.;
				} else if (d >= 0.) {
					sqrDistance = f;	// on V0
					s = 0.; t = 0.;
				} else {
					sqrDistance = f - d*d/a;	// on E0
					s = -d/a; t = 0;
				}
			}
		}
		else { // region 1
			numer = c + e - b - d;
			if (numer <= 0.) {
				sqrDistance = c + 2.*e + f;		// on V2
				s = 0.; t = 1.;
			} else {
				denom = a - 2.*b + c;
				if (numer >= denom) {
					sqrDistance = a + 2.*d + f;	// on V1
					s = 1.; t = 0.;
				} else {
					s = numer/denom;
					t = 1. - s;
					sqrDistance = s*(a*s + b*t + 2.*d) + t*(b*s + c*t + 2.*e) + f;
				}
			}
		}
	} // if (s+t > det)

	baryCoord[0] = 1.0 - s - t;
	baryCoord[1] = s;
	baryCoord[2] = t;

	if (baryCoord[0] < -1e-6 || baryCoord[1] < -1e-6 || baryCoord[2] < -1e-6) {
		cout << "Illegal barycentric coordinate: (" << baryCoord[0] << ", " << baryCoord[1] << ", " << baryCoord[2] << ")" << endl;
	}

	return std::sqrt(std::abs(sqrDistance));
	
}

double CFace::calArea() const
{
	return ZGeom::triArea(m_Vertices[0]->pos(), m_Vertices[1]->pos(), m_Vertices[2]->pos());
}

ZGeom::Vec3d CFace::calcNormal() const
{
    using ZGeom::Vec3d;
    Vec3d v0 = m_Vertices[2]->pos() - m_Vertices[0]->pos();
    Vec3d v1 = m_Vertices[2]->pos() - m_Vertices[1]->pos();
    return (v0 ^ v1).normalize();
}

ZGeom::Vec3d CFace::calBarycenter() const
{
    Vec3d center(0,0,0);
    for (auto p : m_Vertices) center += p->pos();
    return center / (double)m_Vertices.size();
}

// output: parameters of plane function A*x + B*y + C*z + D = 0
std::vector<double> CFace::getPlaneFunction() const
{
	vector<double> para(4);
	ZGeom::Vec3d vNormal = calcNormal();
	para[0] = vNormal[0];
	para[1] = vNormal[1];
	para[2] = vNormal[2];
	double d = vNormal.dot(m_Vertices[0]->pos());
	para[3] = -d;
	return para;
}

std::vector<ZGeom::Vec3d> CFace::getAllVertCoords() const
{
    vector<Vec3d> vCoord;
    for (CVertex *v : m_Vertices) vCoord.push_back(v->pos());
    return vCoord;
}

std::vector<int> CFace::getAllVertIdx() const
{
    vector<int> result;
    for (CVertex* v : m_Vertices) result.push_back(v->getIndex());
    return result;
}

//////////////////////////////////////////////////////
//						CMesh						//
//////////////////////////////////////////////////////

CMesh::CMesh() : m_meshName(""), 
				 m_verbose(true), m_defaultColor(Colorf(0.53f, 0.70f, 0.93f, 1.0f))
{
}

CMesh::CMesh( const CMesh& oldMesh )
{
	cloneFrom(oldMesh);
	m_meshName = oldMesh.m_meshName;
	m_verbose = oldMesh.m_verbose;
}

CMesh::~CMesh()
{
	if (m_verbose) std::cout << "Destroying Mesh '" + m_meshName << "'... ";
	clearMesh();	
	if (m_verbose) std::cout << "Finished!" << std::endl;
}

void CMesh::cloneFrom( const CMesh& oldMesh, const std::string nameSuffix /*=".clone"*/)
{
	if (this == &oldMesh) return;
	clearMesh();

	m_meshName      = oldMesh.m_meshName + nameSuffix;
    m_defaultColor  = oldMesh.m_defaultColor;
    
	for (int i = 0; i < oldMesh.vertCount(); ++i) {
		this->m_vVertices.push_back(new CVertex(*oldMesh.m_vVertices[i]));
	}
	for (int i = 0; i < oldMesh.faceCount(); ++i) {
		this->m_vFaces.push_back(new CFace(*oldMesh.m_vFaces[i]));
	}
	for (int i = 0; i < oldMesh.halfEdgeCount(); ++i) {
		this->m_vHalfEdges.push_back(new CHalfEdge(*oldMesh.m_vHalfEdges[i]));
	}

	for (int i = 0; i < oldMesh.vertCount(); ++i) {
		CVertex* curV = this->m_vVertices[i];
		const CVertex* oldV = oldMesh.m_vVertices[i];
		for (int j = 0; j < oldV->outValence(); ++j) {
			int eidx = oldV->m_HalfEdges[j]->m_eIndex;
			curV->m_HalfEdges.push_back(this->m_vHalfEdges[eidx]);
		}
	}
	for (int i = 0; i < oldMesh.faceCount(); ++i) {
		CFace* curF = this->m_vFaces[i];
		const CFace* oldF = oldMesh.m_vFaces[i];
		for (int j = 0; j < oldF->m_nType; ++j) {
			int vidx = oldF->m_Vertices[j]->m_vIndex;
			int eidx = oldF->m_HalfEdges[j]->m_eIndex;
			curF->m_Vertices.push_back(this->m_vVertices[vidx]);
			curF->m_HalfEdges.push_back(this->m_vHalfEdges[eidx]);
		}
	}
	for (int i = 0; i < oldMesh.halfEdgeCount(); ++i) {
		CHalfEdge* curE = this->m_vHalfEdges[i];
		const CHalfEdge* oldE = oldMesh.m_vHalfEdges[i];
		int vidx0 = oldE->m_Vertices[0]->m_vIndex,
			vidx1 = oldE->m_Vertices[1]->m_vIndex,
			neidx = oldE->m_eNext->m_eIndex,
			peidx = oldE->m_ePrev->m_eIndex,
			fidx  = oldE->m_Face->m_fIndex;

		curE->m_Vertices[0] = this->m_vVertices[vidx0];
		curE->m_Vertices[1] = this->m_vVertices[vidx1];
		curE->m_eNext = this->m_vHalfEdges[neidx];
		curE->m_ePrev = this->m_vHalfEdges[peidx];
		curE->m_Face = this->m_vFaces[fidx];
        if (oldE->m_eTwin == NULL) curE->m_eTwin = NULL;
        else {
			int teidx = oldE->m_eTwin->m_eIndex;
			curE->m_eTwin = this->m_vHalfEdges[teidx];
		}
	}

    copyAttributes(oldMesh.mAttributes);
    setDefaultColor(m_defaultColor);
}

void CMesh::clearMesh()
{
    clearAttributes();

	for (CVertex* v : m_vVertices)	delete v;
	for (CHalfEdge* e : m_vHalfEdges) delete e;
	for (CFace* f : m_vFaces) delete f;
	m_vVertices.clear();
	m_vHalfEdges.clear();
	m_vFaces.clear();

	m_meshName = "";
}

void CMesh::clearAttributes()
{
    for (auto iter = mAttributes.begin(); iter != mAttributes.end(); ++iter)
        delete iter->second;
    mAttributes.clear();

    vector<Colorf> vDefaultColors(vertCount(), m_defaultColor);
    addColorAttr(StrAttrColorDefault, vDefaultColors);
}

void CMesh::load( const std::string& sFileName )
{
	clearMesh();	
	size_t dotPos = sFileName.rfind('.'), slashPos = sFileName.rfind('/');
	m_meshName = sFileName.substr(slashPos+1, dotPos-slashPos-1);
	std::string ext = sFileName.substr(dotPos, sFileName.size() - dotPos);

	if (ext == ".obj" || ext == ".OBJ" || ext == ".Obj") loadFromOBJ(sFileName);
// 	else if (ext == ".m" || ext == ".M") 
//         loadFromM(sFileName);
// 	else if (ext == ".ply" || ext == ".PLY" || ext == ".Ply") 
//         loadFromPLY(sFileName);
//     else if (ext == ".vert" || ext == ".VERT")
//         loadFromVERT(sFileName);
//     else if (ext == ".off" || ext == ".OFF" || ext == ".Off")
//         loadFromOFF(sFileName);
	else throw runtime_error("Unrecognizable file extension!");
}

void CMesh::save(string sFileName)
{
	string sExt = sFileName.substr(sFileName.length()-4);
    if (sExt == ".obj" || sExt == ".OBJ" || sExt == ".Obj") saveToOBJ(sFileName);
// 	sExt = sFileName.substr(sFileName.length()-2);
// 	if (sExt == ".m" || sExt==".M")
// 		saveToM(sFileName);
}

void CMesh::construct(const std::vector<CVertex>& m_pVertex, const std::vector<std::vector<int>>& faceVerts, int nType /*= 3*/)
{
    assert(nType == 3);
///////////////////////////////////////////////////////////////////////////////
///////////////	code to construct primitives pointer vectors	///////////////
	m_vVertices.reserve(m_pVertex.size());
	m_vFaces.reserve(faceVerts.size());
	m_vHalfEdges.reserve(faceVerts.size()*3);

    int nVertex = (int)m_pVertex.size();
    int nFace = (int)faceVerts.size();

    for (int i = 0; i < nVertex; ++i) {
		CVertex* newVertex = new CVertex(m_pVertex[i]);
		newVertex->m_bIsValid = true;
		newVertex->m_HalfEdges.clear();
		m_vVertices.push_back(newVertex);
	}

	for (int i = 0; i < nFace; ++i)
	{
        CFace* newFace = new CFace();
        newFace->create(nType);		
		newFace->m_bIsValid = true;
        this->m_vFaces.push_back(newFace);

		for (int j = 0; j < nType; ++j)
			newFace->m_Vertices.push_back(m_vVertices[faceVerts[i][j]]);

		CHalfEdge* he1 = new CHalfEdge();
		CHalfEdge* he2 = new CHalfEdge();
		CHalfEdge* he3 = new CHalfEdge();
		he1->m_eNext = he2; he1->m_ePrev = he3;
		he2->m_eNext = he3; he2->m_ePrev = he1;
		he3->m_eNext = he1; he3->m_ePrev = he2;
		
		CHalfEdge* heInFace[3] = {he1, he2, he3};
		for (int j = 0; j < 3; ++j) {
			heInFace[j]->m_eTwin = NULL;
			heInFace[j]->m_Face = newFace;
			heInFace[j]->m_bIsValid = true;

			newFace->m_HalfEdges.push_back(heInFace[j]);
			newFace->m_Vertices[j]->m_HalfEdges.push_back(heInFace[j]);
			this->m_vHalfEdges.push_back(heInFace[j]); 
		}
		
		he1->m_Vertices[0] = newFace->m_Vertices[0];
		he1->m_Vertices[1] = newFace->m_Vertices[1];
		he2->m_Vertices[0] = newFace->m_Vertices[1];
		he2->m_Vertices[1] = newFace->m_Vertices[2];
		he3->m_Vertices[0] = newFace->m_Vertices[2];
		he3->m_Vertices[1] = newFace->m_Vertices[0];

		for (int j = 0; j < newFace->m_Vertices[1]->outValence(); ++j)
		{
			if (newFace->m_Vertices[1]->m_HalfEdges[j]->m_Vertices[1] == he1->m_Vertices[0])
			{
				newFace->m_Vertices[1]->m_HalfEdges[j]->m_eTwin = he1;
				he1->m_eTwin = newFace->m_Vertices[1]->m_HalfEdges[j];
				break;
			}
		}
		for (int j = 0; j < newFace->m_Vertices[2]->outValence(); ++j)
		{
			if (newFace->m_Vertices[2]->m_HalfEdges[j]->m_Vertices[1] == he2->m_Vertices[0])
			{
				newFace->m_Vertices[2]->m_HalfEdges[j]->m_eTwin = he2;
				he2->m_eTwin = newFace->m_Vertices[2]->m_HalfEdges[j];
				break;
			}
		}
		for (int j = 0; j < newFace->m_Vertices[0]->outValence(); ++j)
		{
			if (newFace->m_Vertices[0]->m_HalfEdges[j]->m_Vertices[1] == he3->m_Vertices[0])
			{
				newFace->m_Vertices[0]->m_HalfEdges[j]->m_eTwin = he3;
				he3->m_eTwin = newFace->m_Vertices[0]->m_HalfEdges[j];
				break;
			}
		}

	} // for each face

    // re-arrange each vertex's half-edges clockwise
	for (vector<CVertex*>::iterator iter = m_vVertices.begin(); iter != m_vVertices.end();)		
	{
		CVertex* pV = *iter;
		if (pV == NULL) throw std::logic_error("Error: null CVertex pointer encountered!");
		if (pV->outValence() == 0) {
			delete pV;
			iter = m_vVertices.erase(iter);
			continue;
		}		
		++iter;
	} // for each vertex
	
	assignElementsIndex();
}

void CMesh::calAttrFaceNormals()
{
	int faceNum = faceCount();
	std::vector<ZGeom::Vec3d> vFaceNormals(faceNum);
	ZGeom::Vec3d v[2];
	for (int fIndex = 0; fIndex < faceNum; ++fIndex) {
		CFace* face = m_vFaces[fIndex];
		v[0] = face->getVertex(2)->pos() - face->getVertex(0)->pos();
		if (face->m_nType == 3) {
			v[1] = face->getVertex(2)->pos() - face->getVertex(1)->pos();
		} else {
			v[1] = face->getVertex(3)->pos() - face->getVertex(1)->pos();
		}

		vFaceNormals[fIndex] = v[0] ^ v[1];
		vFaceNormals[fIndex].normalize();
	}
	
    addAttr<std::vector<ZGeom::Vec3d>>(vFaceNormals, StrAttrFaceNormal, AR_FACE, AT_VEC_VEC3);
}
 
int CMesh::calEdgeCount()
{
    int halfedgeCount = this->halfEdgeCount();
	int twinEdgeCount = 0;
	for (CHalfEdge* he : m_vHalfEdges) {
		if (he->m_eTwin && he->m_eTwin->m_bIsValid)
            twinEdgeCount++;
	}

	if(twinEdgeCount % 2 != 0)
		throw logic_error("Error: CMesh::getEdgeNum; twinEdgeNum should be even number");

    return halfedgeCount - twinEdgeCount / 2;
}

int CMesh::calAttrBoundaryLoops()
{
    typedef vector<vector<int>> VecVecInt;
    int boundaryNum = 0;
    int vertCount = this->vertCount();
    const vector<bool>& vVertOnBoundary = getVertsOnBoundary();
    
    VecVecInt boundaryEdges;
    vector<bool> vertVisited(vertCount, false);

	for( int i = 0; i < vertCount; i++ ) {
		// find boundary loop from boundary vertex i if it is not in any loop 
        if (!vVertOnBoundary[i] || vertVisited[i]) continue;
        vector<int> edgeLoop;
        int currentIndex = i;
        do {
            vertVisited[currentIndex] = true;
			int edgeIndex = -1;
            for (CHalfEdge* he : m_vVertices[currentIndex]->m_HalfEdges) {					
				if(he->isBoundaryEdge()) {
					edgeIndex = he->getIndex(); break;
				}
			}
			currentIndex = m_vHalfEdges[edgeIndex]->getVertIndex(1);
            edgeLoop.push_back(edgeIndex);
        } while (currentIndex != i);

        boundaryEdges.push_back(edgeLoop);		
	}

    // sort by number of edges of boundaries
    std::sort(boundaryEdges.begin(), boundaryEdges.end(),
        [](const vector<int>& v1, const vector<int>& v2) { return v1.size() < v2.size(); });

    boundaryNum = (int)boundaryEdges.size();    
    addAttr<VecVecInt>(boundaryEdges, StrAttrBoundaryLoops, AR_UNIFORM);
    return boundaryNum;
}

int CMesh::calAttrBoundaryVert()
{
	vector<bool> vVertOnBoundary(vertCount(), false);
	int bNum = 0;
	for(int i = 0; i < vertCount(); ++i) {
		if(m_vVertices[i]->judgeOnBoundary()) {
			vVertOnBoundary[i] = true;
			bNum++;
		} 
	}

    addAttr<std::vector<bool>>(vVertOnBoundary, StrAttrVertOnBoundary, AR_VERTEX, AT_VEC_BOOL);
	addAttr<int>(bNum, StrAttrBoundaryVertCount, AR_UNIFORM, AT_INT);
	return bNum;
}

int CMesh::calEulerNum(  )
{
	int edgeCount = calEdgeCount();
	return vertCount() - edgeCount + faceCount();
}

int CMesh::calMeshGenus(  )
{
	int b = calAttrBoundaryLoops();
	int euler_number = calEulerNum();
	return ( 2 - euler_number - b ) / 2;
}

double CMesh::calVolume() const
{
	double vol = 0.0;
	for(int fi = 0; fi < faceCount(); fi++) {
		CFace* face = m_vFaces[fi];
		const ZGeom::Vec3d& v1 = face->getVertex(0)->pos();
		const ZGeom::Vec3d& v2 = face->getVertex(1)->pos();
		const ZGeom::Vec3d& v3 = face->getVertex(2)->pos();

        ZGeom::Vec3d vn = v1 ^ v2;
		vol += dot(vn, v3);
	}

	vol /= 6;
	return std::fabs(vol);
}

std::vector<int> CMesh::getVertexAdjacentFaceIdx( int vIdx, int ring /*= 1*/ ) const
{
	assert(ring >= 1);
	vector<int> vNeighbors = getVertNeighborVerts(vIdx, ring-1, true);
	
	set<int> markedFaces;
	for (auto iter = begin(vNeighbors); iter != end(vNeighbors); ++iter) {
		const CVertex* pv = vert(*iter);
		for (CHalfEdge* he : pv->m_HalfEdges) {
			const CFace* pf = he->getAttachedFace();
			markedFaces.insert(pf->getFaceIndex());
		}		
	}

	vector<int> vFaces;
	for(int f : markedFaces)
		vFaces.push_back(f);

	return vFaces;
}

void CMesh::assignElementsIndex()
{
	for (int i = 0; i < (int)m_vVertices.size(); ++i)
		m_vVertices[i]->m_vIndex = i;
	
	for (int i = 0; i < (int)m_vHalfEdges.size(); ++i)
		m_vHalfEdges[i]->m_eIndex = i;
	
	for (int i = 0; i < (int)m_vFaces.size(); ++i)
		m_vFaces[i]->m_fIndex = i;
}

bool CMesh::isHalfEdgeMergeable( const CHalfEdge* halfEdge )
{
	const CVertex* v1 = halfEdge->m_Vertices[0], *v2 = halfEdge->m_Vertices[1];
	list<CHalfEdge*> v1HeList, v2HeList;
	list<CVertex*> v1VList, v2VList;

	for (vector<CHalfEdge*>::const_iterator iter = v1->m_HalfEdges.begin(); iter != v1->m_HalfEdges.end(); ++iter) {
		v1HeList.push_back((*iter)->m_eNext);
		v1VList.push_back((*iter)->m_Vertices[1]);
	}
	for (vector<CHalfEdge*>::const_iterator iter = v2->m_HalfEdges.begin(); iter != v2->m_HalfEdges.end(); ++iter) {
		v2HeList.push_back((*iter)->m_eNext);
		v2VList.push_back((*iter)->m_Vertices[1]);
	}
	for (list<CHalfEdge*>::iterator iter1 = v1HeList.begin(); iter1 != v1HeList.end(); ++iter1) {
		for (list<CHalfEdge*>::iterator iter2 = v2HeList.begin(); iter2 != v2HeList.end(); ++iter2) {
			if (*iter1 == *iter2)
				return false;
		}
	}

	CVertex* vOppo1 = halfEdge->m_eNext->m_Vertices[1];
	CVertex* vOppo2(NULL);
	if (halfEdge->m_eTwin && halfEdge->m_eTwin->m_bIsValid)	{
		vOppo2 = halfEdge->m_eTwin->m_eNext->m_Vertices[1];
	}

	for (list<CVertex*>::iterator iter1 = v1VList.begin(); iter1 != v1VList.end(); ++iter1) {
		for (list<CVertex*>::iterator iter2 = v2VList.begin(); iter2 != v2VList.end(); ++iter2) {
			if ( (*iter1 == *iter2) && (*iter1 != vOppo1) && (*iter2 != vOppo2))
				return false;
		}
	}

	return true;
}

void CMesh::scaleAreaToVertexNum()
{
	ZGeom::Vec3d center = calMeshCenter();
	double totalSufaceArea = calSurfaceArea();
	double scale = sqrt( double(vertCount()) / totalSufaceArea );
    scaleAndTranslate(-center, scale);
}

void CMesh::scaleToUnitBox()
{
    move(-calMeshCenter());
    double maxCoord = 0;
    for (int i = 0; i < vertCount(); ++i) {
        for (int j = 0; j < 3; ++j)
            maxCoord = max(maxCoord, fabs(vertPos(i)[j]));
    }
    scaleAndTranslate(ZGeom::Vec3d(0, 0, 0), 1. / maxCoord);
}

void CMesh::scaleEdgeLenToUnit()
{
    ZGeom::Vec3d center = calMeshCenter();
	double length = 0.;
	int edgeNum = 0;
    vector<bool> heVisited(halfEdgeCount(), false);
	for (int i = 0; i < halfEdgeCount(); ++i) {
		if (heVisited[i]) continue;		
		edgeNum++;
		length += m_vHalfEdges[i]->length();
		const CHalfEdge* ptwin = m_vHalfEdges[i]->twinHalfEdge();
		if (ptwin != NULL) heVisited[ptwin->getIndex()] = true;
	}
	
	length /= edgeNum;
	double scale = 1.0 / length;	
    scaleAndTranslate(-center, scale);
}

void CMesh::scaleAndTranslate( ZGeom::Vec3d center, double scale )
{
	for (CVertex* v : m_vVertices)	v->translateAndScale(center, scale);	
}

void CMesh::vertRingNeighborVerts( int vIndex, int ring, std::vector<int>& nbr, bool inclusive /*= false*/ ) const
{	
	std::set<int> snb;
	vertRingNeighborVerts(vIndex, ring, snb, inclusive);	
    nbr = std::vector<int>(snb.begin(), snb.end());
}

void CMesh::vertRingNeighborVerts( int vIndex, int ring, std::set<int>& nbr, bool inclusive /*= false*/ ) const
{
	const CVertex* notei = m_vVertices[vIndex];
	nbr.clear();
	nbr.insert(vIndex);

	std::set<int> nbp = nbr;
	for (int r = 1; r <= ring; ++r) {
		std::set<int> nbn;
		for (int vn : nbp) {
			const CVertex* vStart = m_vVertices[vn];
			for (auto he : vStart->m_HalfEdges) {
				int endv = he->m_Vertices[1]->m_vIndex;
				if (nbr.find(endv) == nbr.end()) {
					nbr.insert(endv);
					nbn.insert(endv);
				}				
				// to avoid boundary vertex being ignored
				if (he->m_eNext) { 
					int endv2 = he->m_eNext->m_Vertices[1]->m_vIndex;
					if (nbr.find(endv2) == nbr.end()) {
						nbr.insert(endv2);
						nbn.insert(endv2);
					}					
				}
			}
		}
		nbp = nbn;
	}
	
	if (!inclusive) nbr.erase(vIndex);
}

bool CMesh::isInNeighborRing( int ref, int query, int ring ) const
{
	assert(ring >= 0);
	if (ref == query) return true;
	std::set<int> iNeighbor;
	vertRingNeighborVerts(ref, ring, iNeighbor, true);
	return (iNeighbor.find(ref) != iNeighbor.end());
}

std::vector<int> CMesh::getVertNeighborVerts( int v, int ring, bool inclusive ) const
{
	if (ring < 0) throw runtime_error("Error: getNeighboringVertex with ring < 0");

	std::vector<int> vn;
	vertRingNeighborVerts(v, ring, vn, inclusive);
	return vn;
}

std::vector<int> CMesh::getVertIsoNeighborVerts( int v, int ring ) const
{
	if (ring < 1) 
		throw std::logic_error("Error: CMesh::getRingVertex with ring < 1");

	std::set<int> v1, v2;
	vertRingNeighborVerts(v, ring-1, v1, true);
	vertRingNeighborVerts(v, ring, v2, true);

	std::vector<int> v3;
	for (int vIdx : v2) {
		if (v1.find(vIdx) == v1.end())
			v3.push_back(vIdx);
	}

	return v3;
}

void CMesh::extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, vector<int>& vFeatures )
{
    const vector<bool>& vertOnBoundary = getVertsOnBoundary();

	const int STATE_IDLE = 0;
	const int STATE_MIN	= -1;
	const int STATE_MAX	=  1;

	assert(vSigVal.size() == vertCount());
	vFeatures.clear();
	
	int state = STATE_IDLE;
	int state_c = STATE_IDLE;
	vector<int> nb;

	double pz = 1e-5;		//1e-5
	double nz = -1e-5;

	for(int j = 0; j < vertCount(); j++)		//m_size: size of the mesh
	{
		state = STATE_IDLE;
		if (vertOnBoundary[j]) 
			continue;  // ignore boundary vertex

		if (vSigVal[j] < lowThresh)				//too small signature discarded
			continue;

		vertRingNeighborVerts(j, ring, nb, false);	//ring == 2
		for (size_t k = 0; k < nb.size(); k++) {	//for each neighbor 
		
			int ev = nb[k];
			state_c = STATE_IDLE;
			if( vSigVal[j] - vSigVal[ev] < 0)		// low bound
				state_c = STATE_MIN;
			else if( vSigVal[j] - vSigVal[ev] > 0)	// high bound
				state_c = STATE_MAX;

			if(state == STATE_IDLE)				    // two-step change
				state = state_c;
			else if( state * state_c <= 0 ) {
				state = STATE_IDLE;
				break;
			}
		}
		if(state == STATE_IDLE) continue;

		vFeatures.push_back(j);
	}
}

void CMesh::extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary/* = 1*/ )
{
    const vector<bool>& vertOnBoundary = getVertsOnBoundary();
    const int STATE_IDLE = 0;
	const int STATE_MIN	= -1;
	const int STATE_MAX	=  1;

	assert(vSigVal.size() == vertCount());
	vFeatures.clear();

	int state = STATE_IDLE;
	int state_c = STATE_IDLE;
	vector<int> nb;

	double pz = 0;		//1e-5
	double nz = -0;

	vector<double> sigDetected;

	for(int j = 0; j < vertCount(); j++)		//m_size: size of the mesh
	{
        if (vertOnBoundary[j]) continue;  // ignore boundary vertex
		if (fabs(vSigVal[j]) < lowThresh)				//too small hks discarded
			continue;
		vertRingNeighborVerts(j, ring, nb, false);	
		
		state = STATE_IDLE;
		for (size_t k = 0; k < nb.size(); k++) {	//for each neighbor 
			int ev = nb[k];
			state_c = STATE_IDLE;
			if( vSigVal[j] - vSigVal[ev] < nz)		// low bound
				state_c = STATE_MIN;
			else if( vSigVal[j] - vSigVal[ev] > pz)	// high bound
				state_c = STATE_MAX;

			if(state == STATE_IDLE)				    // two-step change
				state = state_c;
			else if( state != state_c ) {
				state = STATE_IDLE;
				break;
			}
		}
		if(state == STATE_IDLE) continue;

		vFeatures.push_back(make_pair(j, state_c));	//max: 1, min: -1
		sigDetected.push_back(vSigVal[j]);
	}

	sort(sigDetected.begin(), sigDetected.end(), std::greater<double>());
	std::cout << "Max Signature: ";
	for (int i = 0; i < std::min(10, (int)sigDetected.size()); ++i) std::cout << sigDetected[i] << ' ';
	std::cout << std::endl;
}

bool CMesh::hasBoundary()
{
    if (!hasAttr(StrAttrBoundaryVertCount)) calAttrBoundaryVert();        
    return getAttrValue<int>(StrAttrBoundaryVertCount) > 0;
}

MeshCoordinates CMesh::getVertCoordinates() const
{
    MeshCoordinates coords(vertCount());
    ZGeom::VecNd& vx = coords.getCoordFunc(0);
    ZGeom::VecNd& vy = coords.getCoordFunc(1);
    ZGeom::VecNd& vz = coords.getCoordFunc(2);
    for (int i = 0; i < vertCount(); ++i) {
        auto vCoord = m_vVertices[i]->pos();
        vx[i] = vCoord.x;
        vy[i] = vCoord.y;
        vz[i] = vCoord.z;
    }

    return coords;
}

ZGeom::PointCloud3d CMesh::toPointCloud() const
{
    vector<Vec3d> vPoints;
    for (int i = 0; i < vertCount(); ++i) 
        vPoints.push_back(m_vVertices[i]->pos());    
    return ZGeom::PointCloud3d(vPoints);
}

void CMesh::setVertCoordinates( const MeshCoordinates& coords )
{
	ZGeom::logic_assert(coords.size() == vertCount(), "Size of coordinates and mesh not compatible!");
	std::vector<double> vx = coords.getCoordFunc(0).toStdVector(),
	                    vy = coords.getCoordFunc(1).toStdVector(),
	                    vz = coords.getCoordFunc(2).toStdVector();
	
	setVertexCoordinates(vx, vy, vz);
}

void CMesh::setVertexCoordinates( const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord )
{
	assert(vxCoord.size() == vyCoord.size() && vxCoord.size() == vzCoord.size() && vxCoord.size() == vertCount());

	for (int i = 0; i < vertCount(); ++i)	{
		m_vVertices[i]->setPosition(vxCoord[i], vyCoord[i], vzCoord[i]);
	}
}

void CMesh::setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<ZGeom::Vec3d>& vNewPos)
{
	if(vDeformedIdx.size() != vNewPos.size())
		throw std::logic_error("Error: CMesh::setVertexCoordinates; incompatible parameters");

	size_t vsize = vDeformedIdx.size();
	for (size_t i = 0; i < vsize; ++i)
	{
		m_vVertices[vDeformedIdx[i]]->setPosition(vNewPos[i].x, vNewPos[i].y, vNewPos[i].z);
	}
}

double CMesh::getAvgEdgeLength() const
{
	const MeshAttr<double>* attrAvgEdgeLen = getAttr<double>(StrAttrAvgEdgeLength);
	if (attrAvgEdgeLen == NULL) throw std::logic_error("Attribute of average edge length not available!");
	return attrAvgEdgeLen->attrValue();
}

const std::vector<ZGeom::Vec3d>& CMesh::getFaceNormals()
{
	if (!hasAttr(StrAttrFaceNormal)) calAttrFaceNormals();
	const std::vector<ZGeom::Vec3d>& vFaceNormals = getAttrValue<std::vector<ZGeom::Vec3d>>(StrAttrFaceNormal);
	return vFaceNormals;
}

const std::vector<ZGeom::Vec3d>& CMesh::getVertNormals() const
{
	assert(hasAttr(StrAttrVertNormal));
	return getAttrValue<std::vector<ZGeom::Vec3d>>(StrAttrVertNormal);
}

const std::vector<bool>& CMesh::getVertsOnBoundary()
{
	if (!hasAttr(StrAttrVertOnBoundary)) calAttrBoundaryVert();
	return getAttrValue<std::vector<bool>>(StrAttrVertOnBoundary);
}

bool CMesh::isVertOnBoundary(int vi)
{
    return getVertsOnBoundary()[vi];
}

double CMesh::calSurfaceArea() const
{
	double totalSufaceArea(0);
    for (CFace* f : m_vFaces) totalSufaceArea += f->calArea();
	return totalSufaceArea;
}

ZGeom::Vec3d CMesh::calMeshCenter() const
{
	ZGeom::Vec3d center(0, 0, 0);
	for (int i = 0; i < vertCount(); ++i)
		center += m_vVertices[i]->pos();
	center /= vertCount();
	return center;
}

ZGeom::Vec3d CMesh::calBoundingBox( const ZGeom::Vec3d& center ) const
{
	ZGeom::Vec3d boundBox(0.0, 0.0, 0.0);
	for (int i = 0; i < vertCount(); ++i) {
		boundBox.x = (abs(m_vVertices[i]->m_vPosition.x)>abs(boundBox.x)) ? m_vVertices[i]->m_vPosition.x : boundBox.x;
		boundBox.y = (abs(m_vVertices[i]->m_vPosition.y)>abs(boundBox.y)) ? m_vVertices[i]->m_vPosition.y : boundBox.y;
		boundBox.z = (abs(m_vVertices[i]->m_vPosition.z)>abs(boundBox.z)) ? m_vVertices[i]->m_vPosition.z : boundBox.z;
	}
	boundBox.x = abs(boundBox.x - center.x);
	boundBox.y = abs(boundBox.y - center.y);
	boundBox.z = abs(boundBox.z - center.z);
	return boundBox;
}

void CMesh::saveToMetis( const std::string& sFileName ) const
{
	std::ofstream ofs(sFileName.c_str());
	int vertNum = vertCount();
	ofs << vertNum << '\n';
	for (int i = 0; i < vertNum; ++i) {
		std::vector<int> vNbr = getVertNeighborVerts(i, 1, false);
		int adjCount = (int)vNbr.size();
		for (int j = 0; j < adjCount - 1; ++j)
			ofs << vNbr[j] + 1 << ' ';
		ofs << vNbr[adjCount-1] + 1 << '\n';
	}
	ofs.close();
}

void CMesh::partitionToSubMeshes( const std::vector<std::vector<int>*>& vSubMappedIdx, std::vector<CMesh*>& vSubMeshes ) const
{
	assert(vSubMappedIdx.size() == vSubMeshes.size());
	const int nPart = (int)vSubMappedIdx.size();

	for (int partIdx = 0; partIdx < nPart; ++partIdx)
	{
		CMesh& subMesh = *vSubMeshes[partIdx];
		const std::vector<int>& subMappedIdx = *vSubMappedIdx[partIdx];
		const int subVertCount = (int)subMappedIdx.size();

		std::list<ZGeom::Vec3d> VertexList;	//temporary vertex list
		std::list<int> FaceList;		//temporary face list

		for (int i = 0; i < subVertCount; ++i) {
			VertexList.push_back(this->vertPos(subMappedIdx[i]));
		}
		for (CFace* f : m_vFaces) {
			int fv1 = f->getVertexIndex(0), fv2 = f->getVertexIndex(1), fv3 = f->getVertexIndex(2);
			int sfv1 = -1, sfv2 = -1, sfv3 = -1;
			for (int i = 0; i < subVertCount; ++i) {
				if (subMappedIdx[i] == fv1) sfv1 = i;
				else if (subMappedIdx[i] == fv2) sfv2 = i;
				else if (subMappedIdx[i] == fv3) sfv3 = i;
			}
			if (sfv1 >= 0 && sfv2 >= 0 && sfv3 >= 0) {
				FaceList.push_back(sfv1);
				FaceList.push_back(sfv2);
				FaceList.push_back(sfv3);
			}			
		}

		std::set<int> setIsolatedVert;
		for (int i = 0; i < subVertCount; ++i) setIsolatedVert.insert(i);
		for (auto k : FaceList) setIsolatedVert.erase(k);
		assert(setIsolatedVert.empty());

        int nVertex = (int)VertexList.size();
		int nFace = (int)FaceList.size() / 3;
		int nHalfEdge = 3 * nFace;

        vector<CVertex> m_pVertex(nVertex);
        vector<vector<int>> faceVerts(nFace);
		list<ZGeom::Vec3d>::iterator iVertex = VertexList.begin();
		list<int>::iterator iFace = FaceList.begin();

        for (int i = 0; i < nVertex; i++) {
            m_pVertex[i].m_vPosition = *iVertex++;
            m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
        }
        for (int i = 0; i < nFace; i++) {
            faceVerts[i].resize(3);
            for (int j = 0; j < 3; ++j)
                faceVerts[i][j] = *iFace++;
        }

        subMesh.construct(m_pVertex, faceVerts);
	}	
}

void CMesh::getSubMeshFromFaces( const std::vector<int>& vSubFaces, std::string subMeshName, CMesh& submesh )
{
    submesh.clearMesh();
    submesh.setMeshName(subMeshName);

    std::list<ZGeom::Vec3d> VertexList;
    std::list<int> FaceList;
    std::map<int, int> oldidx2newidx;
    for (int fIdx : vSubFaces) {
        CFace* f = m_vFaces[fIdx];
        int fvOld[3] = { f->getVertexIndex(0), f->getVertexIndex(1), f->getVertexIndex(2) };
        for (int k = 0; k < 3; ++k) {
            int oldVIdx = fvOld[k];
            if (oldidx2newidx.find(oldVIdx) == oldidx2newidx.end()) {
                // add new vertex
                oldidx2newidx[oldVIdx] = (int)VertexList.size();
                VertexList.push_back(m_vVertices[oldVIdx]->pos());
            }
            FaceList.push_back(oldidx2newidx[oldVIdx]);
        }
    }

    int nVertex = (int)VertexList.size();
    int nFace = (int)FaceList.size() / 3;
    vector<CVertex> m_pVertex(nVertex);
    vector<vector<int>> faceVerts(nFace);
    list<ZGeom::Vec3d>::iterator iVertex = VertexList.begin();
    list<int>::iterator iFace = FaceList.begin();
    for (int i = 0; i < nVertex; i++) {
        m_pVertex[i].m_vPosition = *iVertex++;
        m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
    }
    for (int i = 0; i < nFace; i++) {
        faceVerts[i].resize(3);
        for (int j = 0; j < 3; ++j)
            faceVerts[i][j] = *iFace++;
    }

    submesh.construct(m_pVertex, faceVerts);
    submesh.addDefaultColorAttr();
}

void CMesh::getSubMesh(const std::vector<int>& subMappedIdx, std::string subMeshName, CMesh& subMesh)
{
    subMesh.clearMesh();
    subMesh.setMeshName(subMeshName);

    const int subVertCount = (int)subMappedIdx.size();

    std::list<ZGeom::Vec3d> VertexList;	//temporary vertex list
    std::list<int> FaceList;		//temporary face list

    for (int i = 0; i < subVertCount; ++i) {
        VertexList.push_back(this->vertPos(subMappedIdx[i]));
    }
    for (CFace* f : m_vFaces) {
        int fv1 = f->getVertexIndex(0), fv2 = f->getVertexIndex(1), fv3 = f->getVertexIndex(2);
        int sfv1 = -1, sfv2 = -1, sfv3 = -1;
        for (int i = 0; i < subVertCount; ++i) {
            if (subMappedIdx[i] == fv1) sfv1 = i;
            else if (subMappedIdx[i] == fv2) sfv2 = i;
            else if (subMappedIdx[i] == fv3) sfv3 = i;
        }
        if (sfv1 >= 0 && sfv2 >= 0 && sfv3 >= 0) {
            FaceList.push_back(sfv1);
            FaceList.push_back(sfv2);
            FaceList.push_back(sfv3);
        }
    }

    std::set<int> setIsolatedVert;
    for (int i = 0; i < subVertCount; ++i) setIsolatedVert.insert(i);
    for (auto k : FaceList) setIsolatedVert.erase(k);
    assert(setIsolatedVert.empty());

    int nVertex = (int)VertexList.size();
    int nFace = (int)FaceList.size() / 3;
    int nHalfEdge = 3 * nFace;

    vector<CVertex> m_pVertex(nVertex);
    vector<vector<int>> faceVerts(nFace);
    list<ZGeom::Vec3d>::iterator iVertex = VertexList.begin();
    list<int>::iterator iFace = FaceList.begin();

    for (int i = 0; i < nVertex; i++) {
        m_pVertex[i].m_vPosition = *iVertex++;
        m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
    }
    for (int i = 0; i < nFace; i++) {
        faceVerts[i].resize(3);
        for (int j = 0; j < 3; ++j)
            faceVerts[i][j] = *iFace++;
    }

    subMesh.construct(m_pVertex, faceVerts);
    subMesh.addDefaultColorAttr();
}

std::vector<ZGeom::Vec3d> CMesh::allVertPos() const
{
    vector<Vec3d> results(vertCount());
    for (int i = 0; i < vertCount(); ++i) {
        results[i] = m_vVertices[i]->pos();
    }
    return results;
}

void CMesh::loadFromOBJ(std::string sFileName)
{
/* -----  format: smf, obj, dat -----
 * vertex:
 *      v x y z,
 * face(triangle):
 *      f v1 v2 v3  (the vertex index is 1-based)
 * ----------------------------------- */
	//open the file
	FILE *f;
	fopen_s(&f, sFileName.c_str(), "r");
	assert(f != NULL);

	std::list<ZGeom::Vec3d> VertexList;	//temporary vertex list
	std::list<int> FaceList;			//temporary face list

	ZGeom::Vec3d vec;
    char ch = 0;
	int l[3];
	short j;

	ch = fgetc(f);
	while(ch > 0)	//save temporary information of vertex and face in list
	{
		switch(ch)
		{
		case 'v':	//vertex
			ch = fgetc(f);
			if((ch!=' ')&&(ch!='\t'))
				break;
			fscanf_s(f, "%lf%lf%lf", &vec.x, &vec.y, &vec.z);
			VertexList.push_back(vec);

		case 'f':	//face
			ch = fgetc(f);
			if((ch!=' ') && (ch!='\t'))
				break;
			fscanf_s(f, "%ld%ld%ld\n", &l[0], &l[1], &l[2]);
			for(j = 0; j < 3; j++)
				FaceList.push_back(l[j] - 1);		// 0-based, vid - 1
			break;
		
		case '#':
			while(ch != '\n' && ch > 0)
				ch = fgetc(f);
			break;
		}
		ch = fgetc(f);
	}
	fclose(f);
	
    int nVertex = (int)VertexList.size();
    int nFace = (int)FaceList.size() / 3;
    vector<CVertex> m_pVertex(nVertex);
    vector<vector<int>> faceVerts(nFace);
	list<ZGeom::Vec3d>::iterator iVertex = VertexList.begin();
	list<int>::iterator iFace = FaceList.begin();
    for (int i = 0; i < nVertex; i++) {
		m_pVertex[i].m_vPosition = *iVertex++;  
		m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
	}    
    for (int i = 0; i < nFace; i++) {
        faceVerts[i].resize(3);
        for (j = 0; j < 3; ++j)
            faceVerts[i][j] = *iFace++;
    }

	construct(m_pVertex, faceVerts);
}

void CMesh::saveToOBJ( std::string sFileName )
{
	// open the file
	FILE *f = NULL;
	fopen_s(&f, sFileName.c_str(),"wb");
	assert(f != NULL);

	// file header
	fprintf(f, "# vertices : %ld\r\n", vertCount());
	fprintf(f, "# faces    : %ld\r\n", faceCount());
	fprintf(f, "\r\n");
	
	// vertices
	for (int i = 0; i < vertCount(); i++)
	{
        ZGeom::Vec3d vt = m_vVertices[i]->pos();
		fprintf(f, "v %lf %lf %lf\r\n", vt.x, vt.y, vt.z);
	}

	// faces
    for (int i = 0; i < faceCount(); i++) {
        fprintf(f, "f %ld %ld %ld\r\n",
                m_vFaces[i]->getVertexIndex(0) + 1, 
                m_vFaces[i]->getVertexIndex(1) + 1, 
                m_vFaces[i]->getVertexIndex(2) + 1);            
    }

	fclose(f);
}

std::vector<bool> CMesh::getVertsOnHoles()
{
    if (!hasAttr(StrAttrBoundaryLoops)) calAttrBoundaryLoops();
    const vector<vector<int>>& boundaryEdges = getAttrValue<vector<vector<int>>>(StrAttrBoundaryLoops);
    vector<bool> result(vertCount(), false);
    for (int i = 0; i < boundaryEdges.size(); ++i) {
        const vector<int>& vEdgeIdx = boundaryEdges[i];
        if (vEdgeIdx.size() > MAX_HOLE_SIZE) continue;
        auto he = getHalfEdge(vEdgeIdx[i]);
        int pv1 = he->getVertIndex(0), pv2 = he->getVertIndex(1);
        result[pv1] = result[pv2] = true;
    }
    return result;
}

void CMesh::setDefaultColor( ZGeom::Colorf color )
{
    m_defaultColor = color;
    vector<Colorf> vDefaultColors(vertCount(), m_defaultColor);
    addColorAttr(StrAttrColorDefault, ZGeom::ColorSignature(vDefaultColors));
}

void CMesh::addDefaultColorAttr()
{
    vector<Colorf> vDefaultColors(vertCount(), m_defaultColor);
    addColorAttr(StrAttrColorDefault, ZGeom::ColorSignature(vDefaultColors));
}

std::vector<std::vector<int>> CMesh::getBoundaryLoopEdges()
{
    if (!hasAttr(StrAttrBoundaryLoops)) calAttrBoundaryLoops();
    return getAttrValue<vector<vector<int>>>(StrAttrBoundaryLoops);
}

std::vector<std::vector<int>> CMesh::getBoundaryLoopVerts()
{
    const vector<vector<int>> &boundaryEdges = getBoundaryLoopEdges();
    vector<vector<int>> boundaryVerts;
    for (vector<int> ve : boundaryEdges) {
        vector<int> vt;
        for (int e : ve) vt.push_back(m_vHalfEdges[e]->getVertIndex(0));
        boundaryVerts.push_back(vt);
    }
    return boundaryVerts;
}

void CMesh::addAttrMeshFeatures( const vector<int>& featureIdx, const std::string& name )
{
    MeshFeatureList fl;
    for (int vi : featureIdx) fl.addFeature(vi, 0);
    addAttrMeshFeatures(fl, name);
}

void CMesh::faceSplit( int fIdx )
{
    faceSplit(m_vFaces[fIdx]);
}

CVertex* CMesh::faceSplit( CFace* face )
{
    auto center = face->calBarycenter();
    CVertex* vc = new CVertex();
    vc->setPosition(center.x, center.y, center.z);

    CHalfEdge *e12 = face->m_HalfEdges[0], *e23 = face->m_HalfEdges[1], *e31 = face->m_HalfEdges[2];
    CVertex* v1 = e12->m_Vertices[0], *v2 = e23->m_Vertices[0], *v3 = e31->m_Vertices[0];
    CFace *f12c = face, *f23c = new CFace(3), *f31c = new CFace(3);
    CHalfEdge *e1c = new CHalfEdge(), *ec1 = new CHalfEdge(),
              *e2c = new CHalfEdge(), *ec2 = new CHalfEdge(),
              *e3c = new CHalfEdge(), *ec3 = new CHalfEdge();
    assoicateVertEdges(v1, vc, e1c, ec1);
    assoicateVertEdges(v2, vc, e2c, ec2);
    assoicateVertEdges(v3, vc, e3c, ec3);
    makeFace(e12, e2c, ec1, f12c);    
    makeFace(e23, e3c, ec2, f23c);
    makeFace(e31, e1c, ec3, f31c);

    addVertex(vc);
    CHalfEdge *arrEdge[6] = { e1c, ec1, e2c, ec2, e3c, ec3 };
    for (int i = 0; i < 6; ++i) addHalfEdge(arrEdge[i]);    
    addFace(f23c); 
    addFace(f31c);    

    return vc;
}

void CMesh::makeFace( CHalfEdge* e1, CHalfEdge* e2, CHalfEdge* e3, CFace *f )
{
    CHalfEdge::makeLoop(e1, e2, e3);
    e1->m_Face = e2->m_Face = e3->m_Face = f;
    f->create(3);
    f->m_HalfEdges = { e1, e2, e3 };
    f->m_Vertices = { e1->vert(0), e2->vert(0), e3->vert(0) };
}

void CMesh::addVertex( CVertex *v )
{
    v->m_vIndex = (int)m_vVertices.size();
    m_vVertices.push_back(v);
}

void CMesh::addHalfEdge( CHalfEdge *e )
{
    e->m_eIndex = (int)m_vHalfEdges.size();
    m_vHalfEdges.push_back(e);
}


void CMesh::addFace( CFace *f )
{
    f->m_fIndex = (int)m_vFaces.size();
    m_vFaces.push_back(f);
}

void CMesh::assoicateVertEdges( CVertex *v1, CVertex *v2, CHalfEdge *e12, CHalfEdge *e21 )
{
    v1->addHalfEdge(e12);
    v2->addHalfEdge(e21);
    e12->setVerts(v1, v2);
    e21->setVerts(v2, v1);
    CHalfEdge::makeTwins(e12, e21);
}

void CMesh::edgeSwap( CHalfEdge* e1 )
{
    CHalfEdge *te1 = e1->m_eTwin;
    if (te1 == NULL) return;
    CHalfEdge *e2 = e1->nextHalfEdge(), *e3 = e1->prevHalfEdge();
    CHalfEdge *te2 = te1->nextHalfEdge(), *te3 = te1->prevHalfEdge();
    CVertex *v1 = e1->vert(0), *v2 = te1->vert(0);
    CVertex *v3 = e2->vert(1), *v4 = te2->vert(1);
    CFace *f1 = e1->getAttachedFace(), *f2 = te1->getAttachedFace();
    assoicateVertEdges(v4, v3, e1, te1);
    v1->removeHalfEdge(e1); v2->removeHalfEdge(te1);
    makeFace(e1, e3, te2, f1);
    makeFace(te1, te3, e2, f2);
}

bool CMesh::relaxEdge(CHalfEdge* e1)
{
    CHalfEdge *te1 = e1->m_eTwin;
    if (te1 == NULL) return false;
    CHalfEdge *e2 = e1->nextHalfEdge(), *e3 = e1->prevHalfEdge();
    CHalfEdge *te2 = te1->nextHalfEdge(), *te3 = te1->prevHalfEdge();
    CVertex *v1 = e1->vert(0), *v2 = te1->vert(0);
    CVertex *v3 = e2->vert(1), *v4 = te2->vert(1);
    pair<Vec3d, double> center1 = ZGeom::circumcenter(v1->pos(), v2->pos(), v3->pos()),
                        center2 = ZGeom::circumcenter(v1->pos(), v2->pos(), v4->pos());

    // test whether non-mutual vertices lie outside of circum-sphere of opposing triangle
    if ((center1.first - (Vec3d)v4->pos()).length() > center1.second &&
        (center2.first - (Vec3d)v3->pos()).length() > center2.second)
        return false;
    else {
        edgeSwap(e1);
        return true;
    }    
}

// std::vector<double> CMesh::calPrincipalCurvature( int k )
// {
// //    assert(k >= 0 && k < = 2); // k == 0 means total curvature; k_total = k1^2 + k2^2 = 4(kh^2 - kg)
//     int N = vertCount();
//     auto curvMean = getMeanCurvature(), curvGauss = getGaussCurvature();
//     std::vector<double> result(N);
//     for (int i = 0; i < N; ++i) {
//         double delta = std::max(0., curvMean[i]*curvMean[i] - curvGauss[i]);        
//         switch (k) {
//         case 0: result[i] = std::max(0., 4.*curvMean[i]*curvMean[i] - 2.*curvGauss[i]); break;
//         case 1: result[i] = curvMean[i] + sqrt(delta); break;
//         case 2: result[i] = curvMean[i] - sqrt(delta); break;
//         }        
//     }
//     return result;
// }

#if 0
void CMesh::loadFromPLY( std::string sFileName )
{
	FILE *f;
	fopen_s(&f, sFileName.c_str(), "rb");
	assert (f != NULL);

	int p1, p2, p3, i, vNr, pNr;
	char nr; 
	float x, y, z;
	char buffer[101];
	vNr = pNr = 0;

	fseek(f, 0, SEEK_SET);

	// READ IN HEADER
	while (fgetc(f) != '\n') // Reads all of the 1st line
		;

	// read the remainder of the header
	// - interested in 3 parts only
	// - format
	// - number of vertices
	// - number of polygons

	fscanf_s(f, "%100s ", buffer);
	while (_stricmp(buffer, "end_header") != 0)
	{
		if (_stricmp(buffer, "format") == 0)
		{
			fscanf_s(f, "%100s ", buffer);
			if (_stricmp(buffer, "ascii") != 0)
			{
				throw runtime_error("PLY file format error: PLY ASCII support only.");
			}
		} else if (_stricmp(buffer, "element") == 0)
		{
			fscanf_s(f, "%100s ", buffer);
			if (_stricmp(buffer, "vertex") == 0)
			{
				fscanf_s(f, "%100s", buffer);
				vNr = atoi(buffer);
			} else if (_stricmp(buffer, "face") == 0)
			{
				fscanf_s(f, "%100s", buffer);
				pNr = atoi(buffer);
			}
		}
		fscanf_s(f, "%100s ", buffer);
	}

	m_nVertex = vNr;
	m_nFace = pNr;
	m_nHalfEdge = 3 * m_nFace;		//number of half-edges

	m_pVertex = new CVertex[m_nVertex];
	m_pFace = new CFace[m_nFace];

	// READ IN VERTICES
	for (i = 0; i < vNr; i++) // Reads the vertices
	{
		if (fscanf_s(f, "%f %f %f", &x, &y, &z) != 3)
			throw runtime_error("PLY file format error: vertex list");
		while (fgetc(f) != '\n'); // Read till end of the line
		// to skip texture/color values

		m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
		m_pVertex[i].m_vPosition = Vector3D(x, y, z);

	}

	// READ IN POLYGONS
	int triCount(0);
	while ((triCount < pNr) && (!(feof(f)))) 
	{
		fscanf_s(f, "%c", &nr);

		switch (nr)
		{
		case '2':
			throw runtime_error("PLY file format error: polygon 2");

		case '3':
			if (fscanf_s(f, "%d %d %d\n", &p1, &p2, &p3) != 3)
				throw runtime_error( "PLY file format error: polygon 3");
			if ((p1 >= vNr) || (p2 >= vNr) || (p3 >= vNr))
				throw runtime_error( "PLY file format error: vertex index out of range");
			else 
			{ 
				m_pFace[triCount].create(3);
				m_pFace[triCount].m_piVertex[0] = p1;
				m_pFace[triCount].m_piVertex[1] = p2;
				m_pFace[triCount].m_piVertex[2] = p3;
			}
			triCount++;
			break;
		case ' ':
		case '\t':
		case '\n':
			// skip leading whitespace characters
			break;
		default:
			// skip any lines that that we are not interested in
			// (i.e. don't begin with the cases above)
			// e.g. N vertex polygons of form "N v1 v2...vN" on the line
			do {
				nr = fgetc(f);
			} while ((!(feof(f))) && nr != '\n');
			break;
		}


	}
	construct(TODO, TODO);
}

void CMesh::loadFromOFF( std::string sFileName )
{
    // -----  format: off
    //                #v, #f, #v+#f
    //vertex:
    //      x y z,
    //face(triangle):
    //      v1 v2 v3  (the vertex index is 1-based)

    //open the file
    FILE *f;
    fopen_s(&f, sFileName.c_str(), "r");

    char *line[100];
    fscanf_s(f, "%s", line);
    int vNum, fNum, vfNum;
    fscanf_s(f, "%ld%ld%ld", &vNum, &fNum, &vfNum);
    assert(vNum + fNum == vfNum);

    cout << "#vertices" << vNum << "    #faces" << fNum << endl;

    list<Vector3D> VertexList;	//temporary vertex list
    list<int> FaceList;			//temporary face list

    Vector3D vec;

    int l[3];

    for (int i = 0; i < vNum; ++i)
    {
        fscanf_s(f, "%lf%lf%lf", &vec.x, &vec.y, &vec.z);
        VertexList.push_back(vec);
    }
    for (int i = 0; i < fNum; ++i)
    {
        fscanf_s(f, "%ld%ld%ld", l, l+1, l+2);
        for (int j = 0; j < 3; ++j)
            FaceList.push_back(l[j]-1);
    }

    m_nVertex = (int)VertexList.size();
    m_nFace = (int)FaceList.size() / 3;
    m_nHalfEdge = 3 * m_nFace;		//number of half-edges

    //read vertices and faces
    m_pVertex = new CVertex[m_nVertex];
    m_pFace = new CFace[m_nFace];

    list<Vector3D>::iterator iVertex = VertexList.begin();
    list<int>::iterator iFace = FaceList.begin();

    for(int i = 0; i < m_nVertex; i++)
    {
        m_pVertex[i].m_vPosition = *iVertex++;  
        m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
    }

    for(int i = 0; i < m_nFace; i++)
    {
        m_pFace[i].create(3);
        for(int j = 0; j < 3; j++)
            m_pFace[i].m_piVertex[j] = *iFace++;
    }

    VertexList.clear();
    FaceList.clear();
    construct(TODO, TODO);
}

void CMesh::loadFromVERT( std::string sFileName )
{
	//open the file
	ifstream infile(sFileName.c_str());
	if( ! infile ) {
		throw std::runtime_error("error: unable to open input vert file!");
	}
	list<Vector3D> VertexList;//temporary vertex list
	list<int> FaceList;//temporary face list

	string sline; // single line
	double dx,dy,dz;
	while( getline(infile, sline) )
	{
		istringstream sin(sline);
		sin >> dx >> dy >> dz;
		Vector3D vec(dx,dy,dz);
		VertexList.push_back(vec);
	}

	string triName = sFileName.substr(0,sFileName.length()-4);
	triName += "tri";

	cout << triName << endl;
	ifstream trifile(triName.c_str());
	if(!trifile) {
		throw std::runtime_error("unable to open input tri file");
	}

	int tx,ty,tz;
	while( getline(trifile, sline) )
	{
		istringstream sin(sline);
		sin >> tx >> ty >> tz;
		FaceList.push_back(tx-1);
		FaceList.push_back(ty-1);
		FaceList.push_back(tz-1);
	}

	
	m_nVertex = (int)VertexList.size();
	m_nFace = (int)FaceList.size()/3;
	m_nHalfEdge = 3*m_nFace;

	//read vertices and faces
	m_pVertex = new CVertex[m_nVertex];
	m_pFace = new CFace[m_nFace];

	int i;
	list<Vector3D>::iterator iVertex = VertexList.begin();
	list<int>::iterator iFace = FaceList.begin();

	for(i=0; i<m_nVertex; i++)
	{
		m_pVertex[i].m_vPosition = *iVertex++;  
		m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
	}
	short j;
	for(i=0; i<m_nFace; i++)
	{
		m_pFace[i].create(3);
		for(j=0;j<3;j++)
			m_pFace[i].m_piVertex[j] = *iFace++;
	}

	VertexList.clear();
	FaceList.clear();
	construct(TODO, TODO);
}

void CMesh::loadFromM( std::string sFileName )
{
	//open the file
	ifstream infile;
	infile.open(sFileName.c_str());

	vector<Vector3D> VertexList;//temporary vertex list
	//vector<int> VertexIDList;//temporary vertex list
	//vector<Vector3D> VertexNormList;//temporary vertex list
	vector<float> VertexColorList;//temporary vertex list
	//vector<Vector2D> VertexParamList;//temporary vertex list
	vector<int> FaceList;//temporary face list

	const int MAX_TOKEN = 1000;
	char str_buf[MAX_TOKEN];
	char token[MAX_TOKEN];

	do{
		infile.getline(str_buf, MAX_TOKEN);

		if(strstr(str_buf, "Vertex") == 0) {
			;
		}
		else
			break;
	}while(!infile.eof());

	int v_index = 1, prev_idx = 0; // start from 1
	// Process vertex data
	while( strstr(str_buf, "Vertex") ) {
		// move the line to the stream buffer
		istringstream strbuf(str_buf/*, strlen(str_buf)*/);

		strbuf >> token; // "Vertex"

		// vertex index
		int v_idx;
		strbuf >> v_idx;
		prev_idx = v_idx;

		//VertexIDList.push_back(v_idx);
		Vector3D vert;
		strbuf >> vert.x >> vert.y >> vert.z;
		VertexList.push_back(vert);

		char c;
		strbuf >> c;
		if( c == '{' && !strbuf.eof()) {
			while(c != '}' && !strbuf.eof()) {
				strbuf.getline(token, MAX_TOKEN);// >> c;

				//read rgb
				char *pstr = strstr(token,"rgb");
				if(pstr != 0) {
					pstr = strstr(pstr, "=");
					pstr = strstr(pstr, "(");
					pstr = &pstr[1];
					istringstream strbuf(pstr/*, strlen(pstr)*/);
					float r,g,b;
					strbuf >> r >> g >> b;

					//cout << ur << " " << ug << " " << ub << endl;
					VertexColorList.push_back(r);
					VertexColorList.push_back(g);
					VertexColorList.push_back(b);
				}

				//read normal
				/*pstr = strstr(token,"normal");
				if(pstr != 0) {
					pstr = strstr(pstr, "=");
					pstr = strstr(pstr, "(");
					pstr = &pstr[1];
					istringstream strbuf(pstr, strlen(pstr));
					Vector3D norm;
					strbuf >> norm.x >> norm.y >> norm.z;
					VertexNormList.push_back(norm);
				}*/

				//read uv
				/*pstr = strstr(token,"uv");
				if(pstr != 0) {
					pstr = strstr(pstr, "=");
					pstr = strstr(pstr, "(");
					pstr = &pstr[1];
					istringstream strbuf(pstr, strlen(pstr));
					Vector2D uv;
					strbuf >> uv.x >> uv.y;
					VertexParamList.push_back(uv);
				}*/
				//read other properties ...
			}
		}
		infile.getline(str_buf, MAX_TOKEN);
	}

	// process face data
	int f_id = 1;
	while( strstr(str_buf, "Face") ) {
		// move the line to the stream buffer
		istringstream strbuf(str_buf/*, strlen(str_buf)*/);

		strbuf >> token; // "Face"

		int f_num;
		strbuf >> f_num;

		int x_id,y_id,z_id;
		strbuf >> x_id >> y_id >> z_id;

		FaceList.push_back(x_id);
		FaceList.push_back(y_id);
		FaceList.push_back(z_id);

		// face properties
		char c;
		strbuf >> c;
		if( c == '{' && !strbuf.eof()) {			
			while(c != '}' && !strbuf.eof()) {
				strbuf.getline(token, MAX_TOKEN);// >> c;
			}
		}	

		infile.getline(str_buf, MAX_TOKEN);
	}


	m_nVertex = (int)VertexList.size();
	m_nFace = (int)FaceList.size()/3;
	m_nHalfEdge = 3 * m_nFace;

	// read vertices and faces
	m_pVertex = new CVertex[m_nVertex];
	m_pFace = new CFace[m_nFace];

	std::vector<ZGeom::Colorf> vVertColors(m_nVertex);
	for(int i = 0; i < m_nVertex; i++)
	{
		m_pVertex[i].m_vPosition= VertexList[i];  
		m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;

		vVertColors[i][0] = VertexColorList[3*i];
		vVertColors[i][1] = VertexColorList[3*i+1];
		vVertColors[i][2] = VertexColorList[3*i+2];
	}
	
    addColorAttr(StrAttrVertColors, vVertColors);

	for(int i = 0; i < m_nFace; i++)
	{
		m_pFace[i].create(3);
		m_pFace[i].m_piVertex[0] = FaceList[i*3]-1;
		m_pFace[i].m_piVertex[1] = FaceList[i*3+1]-1;
		m_pFace[i].m_piVertex[2] = FaceList[i*3+2]-1;
	}

	construct(TODO, TODO);
}


// Gu Mesh is used as a possible way to keep parameterizations
// and the mesh together
void CMesh::saveToM( const std::string& sFileName )
{
	FILE* fp;
	fopen_s(&fp, sFileName.c_str(), "w+");
	
	const std::vector<Vector3D>& vVertNormals = getVertNormals();
	fprintf(fp, "# vertex=%ld, face=%ld\n", m_nVertex, m_nFace);
	
	std::vector<Colorf>* vColors = NULL;
    if (hasAttr(StrAttrVertColors)) vColors = &getVertColors(StrAttrVertColors);

	for (int i = 0; i < m_nVertex; i++)
	{
		const Vector3D& pos = m_pVertex[i].m_vPosition;
		const Vector3D& normal = vVertNormals[i];
		float r(.5), g(.5), b(.5);
		if (vColors) {
			r = (*vColors)[i].r();
			g = (*vColors)[i].g();
			b = (*vColors)[i].b();
		}
		
		fprintf(fp, "Vertex %ld  %.6f %.6f %.6f {rgb=(%.6f %.6f %.6f) normal=(%.6f %.6f %.6f)",
			    i + 1, pos.x, pos.y, pos.z, r, g, b, normal.x, normal.y, normal.z);
		fprintf(fp, "}\n");
	}

	for (int i = 0; i < m_nFace; i++)
	{
		// assert(m_pFace[i].m_nType == 3);
		fprintf(fp, "Face %ld", i+ 1);

		for (int j = 0; j < m_pFace[i].m_nType; j++)
		{
			int edge = m_pFace[i].m_piEdge[j];
			int _vv_j = m_pHalfEdge[edge].m_iVertex[0] + 1;
			fprintf(fp, " %ld", _vv_j);

		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}
#endif
