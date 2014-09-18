#include "Mesh.h"
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <numeric>
#include "zassert.h"
#include "arithmetic.h"

using namespace std;
using ZGeom::Colorf;

const std::string CMesh::StrAttrBoundaryVertCount		= "mesh_boundary_vert_count";
const std::string CMesh::StrAttrBoundaryCount			= "mesh_boundary_count";
const std::string CMesh::StrAttrAvgEdgeLength			= "mesh_average_edge_length";
const std::string CMesh::StrAttrMeshBBox				= "mesh_bounding_box";
const std::string CMesh::StrAttrMeshCenter				= "mesh_center";
const std::string CMesh::StrAttrVertColors				= "vert_color";
const std::string CMesh::StrAttrVertGaussCurvatures		= "vert_gauss_curvature";
const std::string CMesh::StrAttrVertMeanCurvatures		= "vert_mean_curvature";
const std::string CMesh::StrAttrVertNormal				= "vert_normal";
const std::string CMesh::StrAttrVertOnHole				= "vert_on_hole";
const std::string CMesh::StrAttrVertOnBoundary			= "vert_on_boundary";
const std::string CMesh::StrAttrFaceNormal				= "face_normal";
const std::string CMesh::StrAttrVertMixedArea			= "vert_scalar_mixed_area";

CVertex::CVertex()
{
	init();
}

CVertex::CVertex( double x, double y, double z )
{
	init();
	m_vPosition = Vector3D(x,y,z);
}

CVertex::CVertex( const Vector3D& v )
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

//////////////////////////////////////////////////////
//						CVertex						//
//////////////////////////////////////////////////////

CVertex::~CVertex()
{
	delete[] m_piEdge;
}

void CVertex::init()
{
	m_piEdge = NULL;
	mOutValence = 0;
	m_bIsBoundary = false;
	m_mark = -1;
	m_LocalGeodesic = -1.0;
	m_inheap = false;
	m_bIsValid = true;
}

void CVertex::clone(const CVertex& v)
{
	if (this == &v) return;

	m_vIndex			= v.m_vIndex;
	m_vid				= v.m_vid;
	m_vPosition			= v.m_vPosition;
	mOutValence			= v.mOutValence;
	m_bIsBoundary		= v.m_bIsBoundary;
	m_bIsValid			= v.m_bIsValid;

	m_mark = -1;
	m_LocalGeodesic = -1.0;
	m_inheap = false;

	if (v.m_piEdge != NULL) {
		m_piEdge = new int[mOutValence];		// starting half-edge index array
		for (int i = 0; i < mOutValence; ++i)
			m_piEdge[i] = v.m_piEdge[i];
	}
	else m_piEdge = NULL;
}

std::vector<const CFace*> CVertex::getAdjacentFaces() const
{
	vector<const CFace*> pFaces;
	for (CHalfEdge* he : m_HalfEdges) {
		pFaces.push_back(he->getAttachedFace());
	}

	return pFaces;
}

bool CVertex::judgeOnBoundary()
{
	for (auto iter = m_HalfEdges.begin(); iter != m_HalfEdges.end(); ++iter) {
		if ((*iter)->twinHalfEdge() == NULL || false == (*iter)->twinHalfEdge()->isValid()) {
			m_bIsBoundary = true;
			return true;
		}
	}
	m_bIsBoundary = false;
	return false;
}

void CVertex::translateAndScale( const Vector3D& translation, double s )
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

CHalfEdge::CHalfEdge()
{
	m_iVertex[0] = m_iVertex[1] = -1;
	m_iTwinEdge = m_iNextEdge = m_iPrevEdge = m_iFace = -1;
	
	m_Vertices[0] = m_Vertices[1] = NULL;
	m_eTwin= m_eNext = m_ePrev = NULL;
	m_Face = NULL;

	m_bIsValid = true;
}

CHalfEdge::CHalfEdge( int iV0, int iV1 )
{
	m_iVertex[0] = iV0; m_iVertex[1] = iV1;
	m_iTwinEdge = m_iNextEdge = m_iPrevEdge = m_iFace = -1;

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

//////////////////////////////////////////////////////
//						CHalfEdge					//
//////////////////////////////////////////////////////

CHalfEdge::~CHalfEdge()
{
}

void CHalfEdge::clone(const CHalfEdge& e)
{
	if (this == &e) return;

	m_eIndex = e.m_eIndex;

	m_Vertices[0] = m_Vertices[1] = NULL;
	m_eTwin = m_eNext = m_ePrev = NULL;
	m_Face = NULL;

	m_iVertex[0]	= e.m_iVertex[0];		// starting and ending vertex, Vertex0 -­> Vertex1
	m_iVertex[1]	= e.m_iVertex[1];       //
	m_iTwinEdge		= e.m_iTwinEdge;        // reverse half-edge index, -1 if boundary half edge
	m_iNextEdge		= e.m_iNextEdge;		// next half-edge index ( counter-clock wise )
	m_iPrevEdge		= e.m_iPrevEdge;		// previous half-edge index (cc)
	m_iFace			= e.m_iFace;			// attaching face index ( on the left side )

	m_bIsValid      = e.m_bIsValid;
}


double CHalfEdge::getLength() const
{
	Vector3D v = m_Vertices[0]->getPosition() - m_Vertices[1]->getPosition();
	return v.length();
}


CFace::CFace() : m_nType(0), m_piVertex(NULL), m_piEdge(NULL)
{
	m_bIsValid = true;
}

//////////////////////////////////////////////////////
//						CFace						//
//////////////////////////////////////////////////////

CFace::CFace(int s) : m_nType(s), m_piEdge(NULL), m_piVertex(NULL), m_bIsValid(true)
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
	delete []m_piEdge;
	delete []m_piVertex;
}

void CFace::clone( const CFace& f )
{
	if (this == &f) return;

	m_fIndex	= f.m_fIndex;
	m_nType		= f.m_nType;

	m_piVertex = NULL;
	m_piEdge = NULL;

	if (f.m_piVertex && f.m_piEdge) {
		m_piVertex = new int[m_nType];      // polygon vertex index
		m_piEdge = new int[m_nType];        // polygon edge index

		for(int i = 0; i < m_nType; i++) {
			m_piVertex[i] = f.m_piVertex[i];
			m_piEdge[i] = f.m_piEdge[i];
		}
	}

	m_bIsValid = f.m_bIsValid;
}

void CFace::create(int s)
{
	m_nType = s;
	m_piEdge = new int[s];
	m_piVertex = new int[s];
}

bool CFace::hasVertex( int vidx ) const
{
	for (CVertex* pv : m_Vertices) {
		if (pv->getIndex() == vidx)	return true;
	}

	return false;
}

double CFace::distanceToVertex( const CVertex* vp, std::vector<double>& baryCoord )
{	
	//baryCoord.resize(3, 0);
	
	/**** adapted from WildMagic ****/
	Vector3D V[3] = {m_Vertices[0]->getPosition(), m_Vertices[1]->getPosition(), m_Vertices[2]->getPosition()};
	Vector3D p = vp->getPosition();

	Vector3D diff = V[0] - p;
	Vector3D edge0 = V[1] - V[0];
	Vector3D edge1 = V[2] - V[0];
	double   a = edge0.length2();
	double   b = edge0 * edge1;
	double   c = edge1.length2();
	double   d = diff * edge0;
	double   e = diff * edge1;
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
	return TriArea(m_Vertices[0]->getPosition(), m_Vertices[1]->getPosition(), m_Vertices[2]->getPosition());
}

Vector3D CFace::calcNormal() const
{
	Vector3D v[2];
	v[0] = m_Vertices[2]->getPosition() - m_Vertices[0]->getPosition();
	v[1] = m_Vertices[2]->getPosition() - m_Vertices[1]->getPosition();
	Vector3D vNormal= v[0] ^ v[1];
	vNormal.normalize();	
	return vNormal;
}

Vector3D CFace::calBarycenter() const
{
	return (m_Vertices[0]->getPosition() + m_Vertices[1]->getPosition() + m_Vertices[2]->getPosition()) / 3.0;
}


std::vector<double> CFace::getPlaneFunction()
{
	vector<double> para(4);
	Vector3D vNormal = calcNormal();
	para[0] = vNormal[0];
	para[1] = vNormal[1];
	para[2] = vNormal[2];
	double d = vNormal * m_Vertices[0]->getPosition();
	para[3] = -d;
	return para;
}

//////////////////////////////////////////////////////
//						CMesh						//
//////////////////////////////////////////////////////

CMesh::CMesh() : m_meshName(""), 
				 m_nVertex(0), m_nHalfEdge(0), m_nFace(0), 
				 m_pVertex(NULL), m_pHalfEdge(NULL), m_pFace(NULL),
				 m_bIsPointerVectorExist(false), m_bIsIndexArrayExist(false),
				 m_verbose(true)
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

	m_meshName = oldMesh.m_meshName + nameSuffix;
	m_nVertex = oldMesh.m_nVertex;
	m_nFace = oldMesh.m_nFace;
	m_nHalfEdge = oldMesh.m_nHalfEdge;

	m_bIsPointerVectorExist = oldMesh.m_bIsPointerVectorExist;
	m_bIsIndexArrayExist = oldMesh.m_bIsIndexArrayExist;

	copyAttributes(oldMesh.mAttributes);

	if (m_bIsPointerVectorExist)
	{
		for (int i = 0; i < m_nVertex; ++i) {
			this->m_vVertices.push_back(new CVertex(*oldMesh.m_vVertices[i]));
		}
		for (int i = 0; i < m_nFace; ++i) {
			this->m_vFaces.push_back(new CFace(*oldMesh.m_vFaces[i]));
		}
		for (int i = 0; i < m_nHalfEdge; ++i) {
			this->m_vHalfEdges.push_back(new CHalfEdge(*oldMesh.m_vHalfEdges[i]));
		}

		for (int i = 0; i < m_nVertex; ++i) {
			CVertex* curV = this->m_vVertices[i];
			const CVertex* oldV = oldMesh.m_vVertices[i];
			for (int j = 0; j < oldV->mOutValence; ++j) {
				int eidx = oldV->m_HalfEdges[j]->m_eIndex;
				curV->m_HalfEdges.push_back(this->m_vHalfEdges[eidx]);
			}
		}
		for (int i = 0; i < m_nFace; ++i) {
			CFace* curF = this->m_vFaces[i];
			const CFace* oldF = oldMesh.m_vFaces[i];
			for (int j = 0; j < oldF->m_nType; ++j) {
				int vidx = oldF->m_Vertices[j]->m_vIndex;
				int eidx = oldF->m_HalfEdges[j]->m_eIndex;
				curF->m_Vertices.push_back(this->m_vVertices[vidx]);
				curF->m_HalfEdges.push_back(this->m_vHalfEdges[eidx]);
			}
		}
		for (int i = 0; i < m_nHalfEdge; ++i) {
			CHalfEdge* curE = this->m_vHalfEdges[i];
			const CHalfEdge* oldE = oldMesh.m_vHalfEdges[i];
			int vidx0 = oldE->m_Vertices[0]->m_vIndex,
				vidx1 = oldE->m_Vertices[1]->m_vIndex,
				neidx = oldE->m_eNext->m_eIndex,
				peidx = oldE->m_ePrev->m_eIndex,
				fidx = oldE->m_Face->m_fIndex;

			curE->m_Vertices[0] = this->m_vVertices[vidx0];
			curE->m_Vertices[1] = this->m_vVertices[vidx1];
			curE->m_eNext = this->m_vHalfEdges[neidx];
			curE->m_ePrev = this->m_vHalfEdges[peidx];
			curE->m_Face = this->m_vFaces[fidx];
			if (oldE->m_eTwin != NULL) {
				int teidx = oldE->m_eTwin->m_eIndex;
				curE->m_eTwin = this->m_vHalfEdges[teidx];
			}
			else curE->m_eTwin = NULL;
		}

		buildIndexArrays();
	}
	else if (m_bIsIndexArrayExist)
	{
		m_pVertex = new CVertex[m_nVertex];			//array of vertices
		for(int i = 0; i < m_nVertex; i++)
			m_pVertex[i] = oldMesh.m_pVertex[i];

		m_pHalfEdge = new CHalfEdge[m_nHalfEdge]; 	//array of half-edges
		for(int i = 0; i < m_nHalfEdge; i++)
			m_pHalfEdge[i] = oldMesh.m_pHalfEdge[i];

		m_pFace = new CFace[m_nFace];				//array of faces
		for(int i = 0; i < m_nFace; i++)
			m_pFace[i] = oldMesh.m_pFace[i];
		
		buildPointerVectors();
	}
}

void CMesh::clearMesh()
{
	for (auto iter = mAttributes.begin(); iter != mAttributes.end(); ++iter)
		delete iter->second;
	mAttributes.clear();	

	delete[] m_pVertex;   m_pVertex = NULL;
	delete[] m_pHalfEdge; m_pHalfEdge = NULL;
	delete[] m_pFace;     m_pFace = NULL;	
		
	if (m_bIsPointerVectorExist && m_bSeparateStorage) {
		for (CVertex* v : m_vVertices)	delete v;
		for (CHalfEdge* e : m_vHalfEdges) delete e;
		for (CFace* f : m_vFaces) delete f;
	}
	m_vVertices.clear();
	m_vHalfEdges.clear();
	m_vFaces.clear();

	m_nVertex = m_nHalfEdge = m_nFace = 0;
	m_bIsIndexArrayExist = m_bIsPointerVectorExist = false;
	m_meshName = "";
}

void CMesh::load( const std::string& sFileName )
{
	clearMesh();	
	size_t dotPos = sFileName.rfind('.'), slashPos = sFileName.rfind('/');
	m_meshName = sFileName.substr(slashPos+1, dotPos-slashPos-1);
	std::string ext = sFileName.substr(dotPos, sFileName.size() - dotPos);

	if (ext == ".obj" || ext == ".OBJ" || ext == ".Obj") 
        loadFromOBJ(sFileName);
	else if (ext == ".m" || ext == ".M") 
        loadFromM(sFileName);
	else if (ext == ".ply" || ext == ".PLY" || ext == ".Ply") 
        loadFromPLY(sFileName);
    else if (ext == ".vert" || ext == ".VERT")
        loadFromVERT(sFileName);
    else if (ext == ".off" || ext == ".OFF" || ext == ".Off")
        loadFromOFF(sFileName);
	else throw runtime_error("Unrecognizable file extension!");
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

	char ch = 0;
	
	std::list<Vector3D> VertexList;	//temporary vertex list
	std::list<int> FaceList;			//temporary face list

	Vector3D vec;

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
				FaceList.push_back(l[j]-1);		//vid - 1
			break;
		
		case '#':
			while(ch != '\n' && ch > 0)
				ch = fgetc(f);
			break;
		}
		ch = fgetc(f);
	}
	fclose(f);
		
	m_nVertex = (int)VertexList.size();
	m_nFace = (int)FaceList.size() / 3;
	m_nHalfEdge = 3 * m_nFace;		//number of half-edges

	//read vertices and faces
	m_pVertex = new CVertex[m_nVertex];
	m_pFace = new CFace[m_nFace];

	list<Vector3D>::iterator iVertex = VertexList.begin();
	list<int>::iterator iFace = FaceList.begin();

	for(int i = 0; i < m_nVertex; i++) {
		m_pVertex[i].m_vPosition = *iVertex++;  
		m_pVertex[i].m_vIndex = m_pVertex[i].m_vid = i;
	}

	for(int i = 0; i < m_nFace; i++) {
		m_pFace[i].create(3);
		for(j = 0; j < 3; j++)
			m_pFace[i].m_piVertex[j] = *iFace++;
	}

	construct();
}

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
	construct();
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
	construct();
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

	construct();
}

void CMesh::save(string sFileName)
{
	string sExt = sFileName.substr(sFileName.length()-4);
	if (sExt == ".obj" || sExt==".OBJ")
		saveToOBJ(sFileName);
	sExt = sFileName.substr(sFileName.length()-2);
	if (sExt == ".m" || sExt==".M")
		saveToM(sFileName);
}

void CMesh::saveToOBJ( std::string sFileName )
{
	assert( m_pVertex != NULL && m_pFace != NULL );

	// open the file
	FILE *f = NULL;
	fopen_s(&f, sFileName.c_str(),"wb");
	assert(f != NULL);

	// file header
	fprintf(f, "# vertices : %ld\r\n", m_nVertex);
	fprintf(f, "# faces    : %ld\r\n", m_nFace);
	fprintf(f, "\r\n");

	
	// vertices
	Vector3D center = Vector3D(0, 0, 0); /*getCenter();*/
	for (int i = 0; i < m_nVertex; i++)
	{
		Vector3D vt = m_pVertex[i].m_vPosition - center;
		fprintf(f, "v %lf %lf %lf\r\n", vt.x, vt.y, vt.z);
	}

	// faces
	for (int i = 0; i < m_nFace; i++)
		fprintf(f, "f %ld %ld %ld\r\n", m_pFace[i].m_piVertex[0] + 1, m_pFace[i].m_piVertex[1] + 1, m_pFace[i].m_piVertex[2] + 1);

	fclose(f);
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

void CMesh::findHoles()
{
	std::unordered_set<int> markedVertex;
	std::vector<bool> vIsOnHole(m_nVertex, false);

	for(int i = 0; i < m_nVertex; i++)
	{
		if(!m_vVertices[i]->m_bIsBoundary || markedVertex.find(i) == markedVertex.end()) continue;	// pass if vertex is not boundary or is already visited
		int vi = i;
		const CVertex* pvi =  m_vVertices[i];
		const CHalfEdge* peout = pvi->m_HalfEdges.back();
		vector<int> vTmp;
		markedVertex.insert(vi);
		vTmp.push_back(vi);
		while(peout->m_Vertices[1]->getIndex() != i && vTmp.size() <= MAX_HOLE_SIZE)
		{
			pvi = peout->m_Vertices[1];
			vi = pvi->getIndex();
			peout = pvi->m_HalfEdges.back();
			markedVertex.insert(vi);
			vTmp.push_back(vi);
		}
		if((int)vTmp.size() > MAX_HOLE_SIZE) continue; // boundary	
		
		for (int it_v : vTmp) 
			vIsOnHole[it_v] = true;
	}

	addAttr<std::vector<bool>>(vIsOnHole, StrAttrVertOnHole, AR_VERTEX);

#if 0
	for(int i = 0; i < m_nVertex; i++)
	{
		if(!m_pVertex[i].m_bIsBoundary || m_pVertex[i].m_mark > 0) continue;	// not boundary or visited
		int vi = i;
		int eout = m_pVertex[i].m_piEdge[m_pVertex[i].m_nValence - 1];
		vector<int> v_temp;
		while(m_pHalfEdge[eout].m_iVertex[1] != i)
		{
			v_temp.push_back(vi);
			m_pVertex[vi].m_mark = 1;
			clear_list.push_back(vi);
			vi = m_pHalfEdge[eout].m_iVertex[1];
			eout = m_pVertex[vi].m_piEdge[m_pVertex[vi].m_nValence - 1];
		}
		if((int)v_temp.size() > MAX_HOLE_SIZE) 
			continue; // boundary	
		for(vector<int>::iterator it_v = v_temp.begin(); it_v != v_temp.end(); it_v++) 
			m_pVertex[*it_v].m_bIsHole = true;
	}
	vector<int>::iterator it_vc;	
	for(it_vc = clear_list.begin(); it_vc != clear_list.end(); it_vc++) 
		m_pVertex[*it_vc].m_mark = -1;
#endif
}

void CMesh::construct()
{
	// face normal and area, boundary vertices are computed in the process
	assert(m_pVertex != NULL && m_pFace != NULL);

	//delete old edge list
	delete[] m_pHalfEdge;
	m_pHalfEdge = NULL;	 

///////////////////////////////////////////////////////////////////////////////
///////////////	code to construct primitives pointer vectors	///////////////
	m_vVertices.reserve(m_nVertex);
	m_vFaces.reserve(m_nFace);
	m_vHalfEdges.reserve(m_nHalfEdge);

	for (int i = 0; i < m_nVertex; ++i) {
		CVertex* newVertex = new CVertex(m_pVertex[i]);
		newVertex->m_bIsValid = true;
		newVertex->mOutValence = 0;
		newVertex->m_HalfEdges.clear();
		m_vVertices.push_back(newVertex);
	}

	for (int i = 0; i < m_nFace; ++i)
	{
		CFace* newFace = new CFace(m_pFace[i]);
		this->m_vFaces.push_back(newFace);
		newFace->m_nType = 3;
		newFace->m_bIsValid = true;

		for (int j = 0; j < 3; ++j)
			newFace->m_Vertices.push_back( m_vVertices[m_pFace[i].m_piVertex[j]] );

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
			newFace->m_Vertices[j]->mOutValence++;

			this->m_vHalfEdges.push_back(heInFace[j]); 
		}
		
		he1->m_Vertices[0] = newFace->m_Vertices[0];
		he1->m_Vertices[1] = newFace->m_Vertices[1];
		he2->m_Vertices[0] = newFace->m_Vertices[1];
		he2->m_Vertices[1] = newFace->m_Vertices[2];
		he3->m_Vertices[0] = newFace->m_Vertices[2];
		he3->m_Vertices[1] = newFace->m_Vertices[0];

		for (int j = 0; j < newFace->m_Vertices[1]->mOutValence; ++j)
		{
			if (newFace->m_Vertices[1]->m_HalfEdges[j]->m_Vertices[1] == he1->m_Vertices[0])
			{
				newFace->m_Vertices[1]->m_HalfEdges[j]->m_eTwin = he1;
				he1->m_eTwin = newFace->m_Vertices[1]->m_HalfEdges[j];
				break;
			}
		}
		for (int j = 0; j < newFace->m_Vertices[2]->mOutValence; ++j)
		{
			if (newFace->m_Vertices[2]->m_HalfEdges[j]->m_Vertices[1] == he2->m_Vertices[0])
			{
				newFace->m_Vertices[2]->m_HalfEdges[j]->m_eTwin = he2;
				he2->m_eTwin = newFace->m_Vertices[2]->m_HalfEdges[j];
				break;
			}
		}
		for (int j = 0; j < newFace->m_Vertices[0]->mOutValence; ++j)
		{
			if (newFace->m_Vertices[0]->m_HalfEdges[j]->m_Vertices[1] == he3->m_Vertices[0])
			{
				newFace->m_Vertices[0]->m_HalfEdges[j]->m_eTwin = he3;
				he3->m_eTwin = newFace->m_Vertices[0]->m_HalfEdges[j];
				break;
			}
		}

	} //for each face

	for (vector<CVertex*>::iterator iter = m_vVertices.begin(); iter != m_vVertices.end();)		//--re-arrange each vertex's half-edges clockwise--
	{
		CVertex* pV = *iter;
		if (pV == NULL) throw std::logic_error("Error: null CVertex pointer encountered!");
		if(pV->mOutValence != pV->m_HalfEdges.size())
			throw logic_error("Error: CMesh::construct; pV->m_nValence != pV->m_HalfEdges.size()");

		if (pV->mOutValence == 0) {
			delete pV;
			iter = m_vVertices.erase(iter);
			continue;
		}		
		//pV->judgeOnBoundary();
		++iter;
	} //for each vertex
	
	assignElementsIndex();
	this->m_bIsPointerVectorExist = true;

	buildIndexArrays();
	this->m_bIsIndexArrayExist = true;
}

void CMesh::calFaceNormals()
{
	int faceNum = faceCount();
	std::vector<Vector3D> vFaceNormals(faceNum);
	Vector3D v[2];
	for (int fIndex = 0; fIndex < faceNum; ++fIndex) {
		CFace* face = m_vFaces[fIndex];
		v[0] = face->getVertex(2)->getPosition() - face->getVertex(0)->getPosition();
		if (face->m_nType == 3) {
			v[1] = face->getVertex(2)->getPosition() - face->getVertex(1)->getPosition();
		} else {
			v[1] = face->getVertex(3)->getPosition() - face->getVertex(1)->getPosition();
		}

		vFaceNormals[fIndex] = v[0] ^ v[1];
		vFaceNormals[fIndex].normalize();
	}
	
    addAttr<std::vector<Vector3D>>(vFaceNormals, StrAttrFaceNormal, AR_FACE, AT_VEC_VEC3);
}
 
void CMesh::calVertNormals()
{
	const int vertNum = vertCount();
	const std::vector<Vector3D>& vFaceNormals = getFaceNormals();
	std::vector<Vector3D> vVertNormals(vertNum);

	for (int vIndex = 0; vIndex < vertNum; ++vIndex) {
		CVertex* vertex = m_vVertices[vIndex];
		std::vector<int> incidentFacecIdx = getVertexAdjacentFaceIdx(vIndex, 1);
				
		if(incidentFacecIdx.empty()) {
			vVertNormals[vIndex] = Vector3D(0,0,0);
			continue;
		}

		Vector3D vNormal(0,0,0);
		for(int fIndex : incidentFacecIdx) {
			CFace* face = m_vFaces[fIndex];
			Vector3D cv = (face->getVertex(0)->getPosition() + face->getVertex(1)->getPosition() + face->getVertex(2)->getPosition()) / 3.0;
			double wt = 1.0 / (cv - vertex->getPosition()).length();
			vNormal += vFaceNormals[fIndex] * wt;
		}
		vNormal.normalize();
		vVertNormals[vIndex] = vNormal;
	}

	addAttr<std::vector<Vector3D>>(vVertNormals, StrAttrVertNormal, AR_VERTEX, AT_VEC_VEC3);
}

double CMesh::getHalfEdgeLen( int iEdge ) const
{
	CHalfEdge* he = m_vHalfEdges[iEdge];
	CVertex* verts[2] = {he->m_Vertices[0], he->m_Vertices[1]};
	Vector3D v01 = verts[0]->m_vPosition - verts[1]->m_vPosition;
	return v01.length();
}

int CMesh::calEdgeCount()
{
	assert (this->m_bIsPointerVectorExist);
	int twinEdgeCount = 0;
	for (CHalfEdge* he :m_vHalfEdges) {
		if (he->m_eTwin && he->m_eTwin->m_bIsValid)
			twinEdgeCount++;
	}

	if(twinEdgeCount % 2 != 0)
		throw logic_error("Error: CMesh::getEdgeNum; twinEdgeNum should be even number");

	return halfEdgeCount() - twinEdgeCount / 2;
}

int CMesh::calBoundaryNum()
{
	int boundaryNum = 0;
	std::set<int> boundaryIndexSet;
	for( int i = 0; i < m_nVertex; i++ ) {
		if( m_vVertices[i]->m_bIsBoundary )
			boundaryIndexSet.insert( i );
	}

	for( int i = 0; i < m_nVertex; i++ ) {
		// find boundary loop from boundary vertex i if it is not in any loop 
		if( m_vVertices[i]->m_bIsBoundary && boundaryIndexSet.find(i) != boundaryIndexSet.end()) {
			int currentIndex = i;
			int nextIndex = i;
			do {
				currentIndex = nextIndex;
				auto it = boundaryIndexSet.find( currentIndex );
				boundaryIndexSet.erase( it );
				
				int edgeIndex = -1;
				for( auto he_iter = m_vVertices[i]->m_HalfEdges.begin(); he_iter != m_vVertices[i]->m_HalfEdges.end(); ++he_iter ) 
				{
					CHalfEdge* he = *he_iter;
					if( he->twinHalfEdge() == NULL ) {
						edgeIndex = he->getIndex();
						break;
					}
				}
				nextIndex = m_vHalfEdges[edgeIndex]->vert(1)->getIndex();
			} while( nextIndex != i );

			boundaryNum++;
		}
	}

	addAttr<int>(boundaryNum, StrAttrBoundaryCount, AR_UNIFORM, AT_INT);
	return boundaryNum;
}

int CMesh::calBoundaryVert()
{
	std::vector<bool> vVertOnBoundary(m_nVertex, false);
	int bNum = 0;
	for(int i = 0; i < m_nVertex; ++i)
	{
		if(m_vVertices[i]->judgeOnBoundary() ) {
			vVertOnBoundary[i] = true;
			bNum++;
		} else vVertOnBoundary[i] = false;
	}
    addAttr<std::vector<bool>>(vVertOnBoundary, StrAttrVertOnBoundary, AR_VERTEX, AT_VEC_BOOL);
	addAttr<int>(bNum, StrAttrBoundaryVertCount, AR_UNIFORM, AT_INT);

	if (m_bIsIndexArrayExist) {
		for (int i = 0; i < m_nVertex; ++i) 
			m_pVertex[i].m_bIsBoundary = vVertOnBoundary[i];
	}

	return bNum;
}

int CMesh::calEulerNum(  )
{
	int edgeCount = calEdgeCount();
	return m_nVertex - edgeCount + m_nFace;
}

int CMesh::calMeshGenus(  )
{
	int b = calBoundaryNum();
	int euler_number = calEulerNum();
	return ( 2 - euler_number - b ) / 2;
}

double CMesh::calAreaMixed(double a, double b, double c, double& cotan_a, double& cotan_c)
{
	double cosa = (b*b+c*c-a*a)/(2.0*b*c);
	double cosc = (b*b+a*a-c*c)/(2.0*b*a);
	cotan_a = cosa / sqrt(1.0 - cosa*cosa);
	cotan_c = cosc / sqrt(1.0 - cosc*cosc);

	if (a*a + c*c < b*b)
	{
		double s = (a+b+c)/2.0;
		return sqrt(s*(s-a)*(s-b)*(s-c))/2.0;
	}
	else if (a*a + b*b < c*c || b*b + c*c < a*a)
	{
		double s = (a+b+c)/2.0;
		return sqrt(s*(s-a)*(s-b)*(s-c))/4.0;
	}
	else 
	{
		return (a*a*cotan_a + c*c*cotan_c)/8.0;
	}
}


// Calculate mean and Gauss curvatures                                                                       */
//
void CMesh::calCurvatures()
{
	const double pi = ZGeom::PI;
	const int vertNum = this->vertCount();
	std::vector<double> vGaussCurvatures(vertNum);
	std::vector<double> vMeanCurvatures(vertNum);

	for (int vIndex = 0; vIndex < vertNum; ++vIndex) 
	{
		CVertex* vi = m_vVertices[vIndex];
		double sum = 0.0;		// sum of attaching corner's angle
		double amix = 0.0;
		Vector3D kh;

		if(vi->isOnBoundary()) {
			// boundary vertex has zero curvature
			vGaussCurvatures[vIndex] = 0.0;	
			vMeanCurvatures[vIndex] = 0.0;
			continue;
		}

		for (auto he = vi->m_HalfEdges.begin(); he != vi->m_HalfEdges.end(); ++he) {
			CHalfEdge* e0 = *he;
			CHalfEdge* e1 = e0->m_eNext;
			CHalfEdge* e2 = e1->m_eNext;
			double len0 = e0->getLength();
			double len1 = e1->getLength();
			double len2 = e2->getLength();

			// compute corner angle by cosine law 
			double corner = std::acos((len0*len0 + len2*len2 - len1*len1) / (2.0*len0*len2));
			sum += corner;
			double cota, cotc;
			amix += calAreaMixed(len0, len1, len2, cota, cotc);

			const CVertex* pt1 = e1->vert(0);
			const CVertex* pt2 = e1->vert(1);
			kh += (vi->getPosition() - pt1->getPosition()) * cota + (vi->getPosition() - pt2->getPosition()) * cotc;
		}

		vGaussCurvatures[vIndex] = (2.0 * pi - sum) / amix;
		kh = kh / (2.0 * amix);
		vMeanCurvatures[vIndex] = kh.length() / 2.0;	//half magnitude of 
	}

	addAttr<std::vector<double>>(vGaussCurvatures, StrAttrVertGaussCurvatures, AR_VERTEX, AT_VEC_DBL);
	addAttr<std::vector<double>>(vMeanCurvatures, StrAttrVertMeanCurvatures, AR_VERTEX, AT_VEC_DBL);
}

double CMesh::calHalfAreaMixed( double a, double b, double c, double& cotan_a )
{
	if ( a*a + c*c < b*b ) {
		double s = (a+b+c)/2.0;
		return sqrt(s*(s-a)*(s-b)*(s-c))/4.0;
	} else if ( a*a + b*b < c*c || b*b + c*c < a*a) {
		double s = (a+b+c) / 2.0;
		return sqrt(s*(s-a)*(s-b)*(s-c)) / 8.0;
	} else {
		double cosa = (b*b + c*c - a*a) / (2.0*b*c);
		cotan_a = cosa / sqrt(1 - cosa*cosa);
		return (a*a*cotan_a) / 8.0;
	}
}

bool CMesh::calVertexArea(vector<double>& Av)
{
	Av.resize(m_nVertex,0.0);
	for(int f = 0; f < m_nFace; f++)
	{
		int i0 = m_pFace[f].m_piVertex[0];
		int i1 = m_pFace[f].m_piVertex[1];
		int i2 = m_pFace[f].m_piVertex[2];
		Vector3D v0 = m_pVertex[i0].m_vPosition;
		Vector3D v1 = m_pVertex[i1].m_vPosition;
		Vector3D v2 = m_pVertex[i2].m_vPosition;
		v1 = v1 - v0;
		v2 = v2 - v0;
		crossProduct3D(v1, v2, v0);		//v0 = cross(v1,v2)
		double ff = v0.length()/(6.0);//*m_edge*m_edge);
		Av[i0] += ff;
		Av[i1] += ff;
		Av[i2] += ff;
	}
	return true;
}

double CMesh::calLocalGeodesic( int ia, int ib, int ic ) const
{
	// ia - vertex with smaller geodesic; ib - with greater geodesic; ic - update
	double la = (m_pVertex[ib].m_vPosition-m_pVertex[ic].m_vPosition).length();
	double lb = (m_pVertex[ia].m_vPosition-m_pVertex[ic].m_vPosition).length();
	double lc = (m_pVertex[ib].m_vPosition-m_pVertex[ia].m_vPosition).length();
	double ctheta = (la*la+lb*lb-lc*lc)/(2*la*lb);
	double stheta = sqrt(1-ctheta*ctheta);
	double u = m_pVertex[ib].m_LocalGeodesic - m_pVertex[ia].m_LocalGeodesic;
	double ld = lb-la*ctheta;
	double le = la*stheta;
	double delta = lc*lc-u*u;
	double tc = m_pVertex[ic].m_LocalGeodesic;
	if (delta >= 0.0) {
		delta = sqrt(delta);
		double t1 = lb*(u*ld+le*delta)/(lc*lc);//(-B+delta)/(2*A);
		if(t1>u && lb*(t1-u)/t1>la*ctheta && lb*(t1-u)/t1<la/abs(ctheta))
		{
			if(tc<0.0) tc = t1+m_pVertex[ia].m_LocalGeodesic;
			else tc = min(tc,t1+m_pVertex[ia].m_LocalGeodesic);
		}
		else
		{
			double minab = min(lb+m_pVertex[ia].m_LocalGeodesic,la+m_pVertex[ib].m_LocalGeodesic);
			if(tc<0.0) tc = minab;
			else tc = min(tc,minab);
		}
	}
	else {
		double minab = min(lb+m_pVertex[ia].m_LocalGeodesic,la+m_pVertex[ib].m_LocalGeodesic);
		if(tc<0.0) tc = minab;
		else tc = min(tc,minab);
	}
	return tc;
}

bool CMesh::vertGeoNeighborVerts(int i, double ring, vector<GeoNote>& nbg)
{
	GeoQueue heapqueue;

	if(!nbg.empty()) nbg.clear();

	CVertex& notei = m_pVertex[i];
	notei.m_mark = i;
	notei.m_LocalGeodesic = 0.0;
	notei.m_inheap = true;
	nbg.push_back(GeoNote(i,0.0));

	int size = notei.mOutValence;
	int j,k,ia,ib,ic;

	bool flag = true;

	for (j = 0; j < size; j++) {
		int ee = notei.m_piEdge[j];
		int endv = m_pHalfEdge[ee].m_iVertex[1];
		m_pVertex[endv].m_inheap = true;
		Vector3D vt = m_pVertex[endv].m_vPosition - m_pVertex[i].m_vPosition;
		double mgeo = vt.length();
		m_pVertex[endv].m_LocalGeodesic = mgeo;   // geodesic in first ring
		m_pVertex[endv].m_mark = i;
		if(mgeo<ring) nbg.push_back(GeoNote(endv,mgeo));
	}

	for (j = 0; j < size; j++) {
		int e1 = notei.m_piEdge[j];
		ia = m_pHalfEdge[e1].m_iVertex[1];
		int e2 = m_pHalfEdge[e1].m_iNextEdge;
		ib = m_pHalfEdge[e2].m_iVertex[1];
		if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic) {
			int it = ia;
			ia = ib;
			ib = it;
		}
		e1 = m_pHalfEdge[e2].m_iTwinEdge;
		if(e1<0) continue;
		e2 = m_pHalfEdge[e1].m_iNextEdge;
		ic = m_pHalfEdge[e2].m_iVertex[1];
		double mgeo = calLocalGeodesic(ia,ib,ic);
		m_pVertex[ic].m_LocalGeodesic = mgeo;
		m_pVertex[ic].m_inheap = true;
		if(mgeo<ring) heapqueue.push(GeoNote(ic,mgeo));
	} // first ring

	int itr = 0;
	while (!heapqueue.empty())// && itr<MAX_NEIGHBOR_NUMBER) 
	{
		itr++;
		GeoNote nt = heapqueue.top();
		heapqueue.pop();

		int sg = nt.m_id;
		double sgd = nt.m_geodesic;
		if(m_pVertex[sg].m_bIsBoundary) flag = false;

		if(m_pVertex[sg].m_mark==i) continue;  // marched already
		//if(m_pOctave[o].m_pNote[sg].m_LocalGeodesic < sgd) continue;
		if(m_pVertex[sg].m_LocalGeodesic > ring) break;   // reach the upper bound
		nbg.push_back(nt);
		m_pVertex[sg].m_mark = i;

		// update adjacent vertices of sg
		for (k=0; k<m_pVertex[sg].mOutValence; k++) {
			ia = sg;
			int e1 = m_pVertex[sg].m_piEdge[k];
			ib = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ib].m_mark != i) continue; // unreached point
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				ia = ib;
				ib = sg;
			}
			ic = m_pHalfEdge[e2].m_iVertex[1];
			if(m_pVertex[ic].m_mark != i) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					if(gg<ring) heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					if(gg<ring) heapqueue.push(GeoNote(ic,gg));
				}
			}
			e2 = m_pHalfEdge[e1].m_iTwinEdge;
			if(e2<0 || e2>=m_nHalfEdge) continue;
			e1 = m_pHalfEdge[e2].m_iNextEdge;
			ic = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ic].m_mark != i) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					if(gg<ring) heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					if(gg<ring) heapqueue.push(GeoNote(ic,gg));
				}
			}
		}
	}

	for (size_t ni=0; ni<nbg.size(); ni++)
	{
		int pos = nbg[ni].m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}
	// clear heap
	while (!heapqueue.empty())
	{
		GeoNote nt = heapqueue.top();
		heapqueue.pop();
		int pos = nt.m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}

	return flag;
}

double CMesh::calGeodesic( int s, int t ) const
{
	if(s == t) return 0.0;

	GeoQueue heapqueue;

	vector<GeoNote> nbg;
	CVertex& notei = m_pVertex[s];
	notei.m_mark = s;
	notei.m_LocalGeodesic = 0.0;
	notei.m_inheap = true;
	nbg.push_back(GeoNote(s,0.0));

	int size = notei.mOutValence;
	int j,k,ia,ib,ic;

	bool stop = false;
	double geo = 0.0;

	for (j=0; j<size; j++)
	{
		int ee = notei.m_piEdge[j];
		int endv = m_pHalfEdge[ee].m_iVertex[1];
		Vector3D vt = m_pVertex[endv].m_vPosition - m_pVertex[s].m_vPosition;
		double mgeo = vt.length();
		
		if(endv == t) {stop=true; geo=mgeo; break;}	// destination reached
		m_pVertex[endv].m_inheap = true;
		m_pVertex[endv].m_LocalGeodesic = mgeo;   // geodesic in first ring
		m_pVertex[endv].m_mark = s;
		nbg.push_back(GeoNote(endv,mgeo));
	}
	if(!stop)
	{
		for (j=0; j<size; j++)
		{
			int e1 = notei.m_piEdge[j];
			ia = m_pHalfEdge[e1].m_iVertex[1];
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			ib = m_pHalfEdge[e2].m_iVertex[1];
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				int it = ia;
				ia = ib;
				ib = it;
			}
			e1 = m_pHalfEdge[e2].m_iTwinEdge;
			if(e1<0) continue;
			e2 = m_pHalfEdge[e1].m_iNextEdge;
			ic = m_pHalfEdge[e2].m_iVertex[1];
			double mgeo = calLocalGeodesic(ia,ib,ic);
			m_pVertex[ic].m_LocalGeodesic = mgeo;
			m_pVertex[ic].m_inheap = true;
			heapqueue.push(GeoNote(ic,mgeo));
		} // first ring
	}
	
	int count = 0;
	while (!stop && !heapqueue.empty())
	{
		//if(++count == m_nVertex) {break;}
		//cout << ++count << endl;
		GeoNote nt = heapqueue.top();
		heapqueue.pop();

		int sg = nt.m_id;
		if(sg == t) 
		{
			stop=true; 
			geo=nt.m_geodesic; 
			nbg.push_back(nt); 
			break;
		}
		
		double sgd = nt.m_geodesic;
//		if(m_pVertex[sg].m_bIsBoundary) continue;
		if(m_pVertex[sg].m_mark==s) continue;  // matched already
		m_pVertex[sg].m_mark = s;
		nbg.push_back(nt);
		// update adjacent vertices of sg
		for (k=0; k<m_pVertex[sg].mOutValence; k++)
		{
			ia = sg;
			int e1 = m_pVertex[sg].m_piEdge[k];
			ib = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ib].m_mark != s) continue; // unreached point
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				ia = ib;
				ib = sg;
			}
			ic = m_pHalfEdge[e2].m_iVertex[1];
			if(m_pVertex[ic].m_mark != s) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
			}
			e2 = m_pHalfEdge[e1].m_iTwinEdge;
			if(e2<0 || e2>=m_nHalfEdge) continue;
			e1 = m_pHalfEdge[e2].m_iNextEdge;
			ic = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ic].m_mark != s) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
			}
		}
	}

	for (size_t ni=0; ni<nbg.size(); ni++)
	{
		int pos = nbg[ni].m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}
	// clear heap
	while (!heapqueue.empty())
	{
		GeoNote nt = heapqueue.top();
		heapqueue.pop();
		int pos = nt.m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}
	nbg.clear();
	return geo;
}

double CMesh::getGeodesicToBoundary(int s) const
{
	GeoQueue heapqueue;

	vector<GeoNote> nbg;
	CVertex& notei = m_pVertex[s];
	notei.m_mark = s;
	notei.m_LocalGeodesic = 0.0;
	notei.m_inheap = true;
	nbg.push_back(GeoNote(s,0.0));

	int size = notei.mOutValence;
	int j, k, ia, ib, ic;

	bool stop = false;
	double geo = 0.0;

	const std::vector<bool>& vVertIsHole = getVertsOnHole_const();

	for (j = 0; j < size; j++)
	{
		int ee = notei.m_piEdge[j];
		int endv = m_pHalfEdge[ee].m_iVertex[1];
		Vector3D vt = m_pVertex[endv].m_vPosition - m_pVertex[s].m_vPosition;
		double mgeo = vt.length();

		if(m_pVertex[endv].m_bIsBoundary && !vVertIsHole[endv]) {stop=true; geo=mgeo; break;}	// destination reached
		m_pVertex[endv].m_inheap = true;
		m_pVertex[endv].m_LocalGeodesic = mgeo;   // geodesic in first ring
		m_pVertex[endv].m_mark = s;
		nbg.push_back(GeoNote(endv,mgeo));
	}
	if(!stop)
	{
		for (j = 0; j < size; j++)
		{
			int e1 = notei.m_piEdge[j];
			ia = m_pHalfEdge[e1].m_iVertex[1];
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			ib = m_pHalfEdge[e2].m_iVertex[1];
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				int it = ia;
				ia = ib;
				ib = it;
			}
			e1 = m_pHalfEdge[e2].m_iTwinEdge;
			if(e1<0) continue;
			e2 = m_pHalfEdge[e1].m_iNextEdge;
			ic = m_pHalfEdge[e2].m_iVertex[1];
			double mgeo = calLocalGeodesic(ia,ib,ic);
			m_pVertex[ic].m_LocalGeodesic = mgeo;
			m_pVertex[ic].m_inheap = true;
			heapqueue.push(GeoNote(ic,mgeo));
		} // first ring
	}

	int count = -m_nVertex;
	while (!stop && !heapqueue.empty())
	{
		if(++count == m_nVertex) {break;}

		GeoNote nt = heapqueue.top();
		heapqueue.pop();

		int sg = nt.m_id;
		if(m_pVertex[sg].m_bIsBoundary && !vVertIsHole[sg]) 
		{
			stop = true; 
			geo = nt.m_geodesic; 
			nbg.push_back(nt); 
			break;
		}

		double sgd = nt.m_geodesic;
		
		if(m_pVertex[sg].m_mark == s) continue;  // marched already
		m_pVertex[sg].m_mark = s;
		if(vVertIsHole[sg]) continue;
		nbg.push_back(nt);
		// update adjacent vertices of sg
		for (k=0; k<m_pVertex[sg].mOutValence; k++)
		{
			ia = sg;
			int e1 = m_pVertex[sg].m_piEdge[k];
			ib = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ib].m_mark != s) continue; // unreached point
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				ia = ib;
				ib = sg;
			}
			ic = m_pHalfEdge[e2].m_iVertex[1];
			if(m_pVertex[ic].m_mark != s) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
			}
			e2 = m_pHalfEdge[e1].m_iTwinEdge;
			if(e2<0 || e2>=m_nHalfEdge) continue;
			e1 = m_pHalfEdge[e2].m_iNextEdge;
			ic = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ic].m_mark != s) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
			}
		}
	}

	for (size_t ni=0; ni<nbg.size(); ni++)
	{
		int pos = nbg[ni].m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}
	// clear heap
	while (!heapqueue.empty())
	{
		GeoNote nt = heapqueue.top();
		heapqueue.pop();
		int pos = nt.m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}

	//if(geo==0.0) cout<<s<<endl;
	return geo;
}

double CMesh::getGeodesicToBoundary(int s, vector<GeoNote>& nbg)
{
	GeoQueue heapqueue;
	if(!nbg.empty()) nbg.clear();
	
	CVertex& notei = m_pVertex[s];
	notei.m_mark = s;
	notei.m_LocalGeodesic = 0.0;
	notei.m_inheap = true;
	nbg.push_back(GeoNote(s,0.0));

	int size = notei.mOutValence;
	int j,k,ia,ib,ic;

	bool stop = false;
	double geo = 0.0;
	const std::vector<bool>& vVertIsHole = getVertsOnHole();

	for (j=0; j<size; j++)
	{
		int ee = notei.m_piEdge[j];
		int endv = m_pHalfEdge[ee].m_iVertex[1];
		Vector3D vt = m_pVertex[endv].m_vPosition - m_pVertex[s].m_vPosition;
		double mgeo = vt.length();

		if(m_pVertex[endv].m_bIsBoundary && !vVertIsHole[endv]) {stop=true; geo=mgeo; break;}	// destination reached
		m_pVertex[endv].m_inheap = true;
		m_pVertex[endv].m_LocalGeodesic = mgeo;   // geodesic in first ring
		m_pVertex[endv].m_mark = s;
		nbg.push_back(GeoNote(endv,mgeo));
	}
	if(!stop)
	{
		for (j=0; j<size; j++)
		{
			int e1 = notei.m_piEdge[j];
			ia = m_pHalfEdge[e1].m_iVertex[1];
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			ib = m_pHalfEdge[e2].m_iVertex[1];
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				int it = ia;
				ia = ib;
				ib = it;
			}
			e1 = m_pHalfEdge[e2].m_iTwinEdge;
			if(e1<0) continue;
			e2 = m_pHalfEdge[e1].m_iNextEdge;
			ic = m_pHalfEdge[e2].m_iVertex[1];
			double mgeo = calLocalGeodesic(ia,ib,ic);
			m_pVertex[ic].m_LocalGeodesic = mgeo;
			m_pVertex[ic].m_inheap = true;
			heapqueue.push(GeoNote(ic,mgeo));
		} // first ring
	}

	int count = 0;
	while (!stop && !heapqueue.empty())
	{
		if(++count == m_nVertex) count = count;
		if(s==2576) 
		{
			cout<<count<<endl;
		}
		GeoNote nt = heapqueue.top();
		heapqueue.pop();

		int sg = nt.m_id;
		if(m_pVertex[sg].m_bIsBoundary && !vVertIsHole[sg]) 
		{
			stop=true; 
			geo=nt.m_geodesic; 
			nbg.push_back(nt); 
			break;
		}

		double sgd = nt.m_geodesic;
		
		if(m_pVertex[sg].m_mark==s) continue;  // marched already
		m_pVertex[sg].m_mark = s;
		if(vVertIsHole[sg]) continue;

		nbg.push_back(nt);
		// update adjacent vertices of sg
		for (k=0; k<m_pVertex[sg].mOutValence; k++)
		{
			ia = sg;
			int e1 = m_pVertex[sg].m_piEdge[k];
			ib = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ib].m_mark != s) continue; // unreached point
			int e2 = m_pHalfEdge[e1].m_iNextEdge;
			if (m_pVertex[ia].m_LocalGeodesic > m_pVertex[ib].m_LocalGeodesic)
			{
				ia = ib;
				ib = sg;
			}
			ic = m_pHalfEdge[e2].m_iVertex[1];
			if(m_pVertex[ic].m_mark != s) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
			}
			e2 = m_pHalfEdge[e1].m_iTwinEdge;
			if(e2<0 || e2>=m_nHalfEdge) continue;
			e1 = m_pHalfEdge[e2].m_iNextEdge;
			ic = m_pHalfEdge[e1].m_iVertex[1];
			if(m_pVertex[ic].m_mark != s) 
			{
				double gg = calLocalGeodesic(ia,ib,ic);   // update geodesic
				if(m_pVertex[ic].m_LocalGeodesic < 0.0f)
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
				else if(gg < m_pVertex[ic].m_LocalGeodesic) // heaped, shorter patch came
				{
					m_pVertex[ic].m_LocalGeodesic = gg;
					m_pVertex[ic].m_inheap = true;
					heapqueue.push(GeoNote(ic,gg));
				}
			}
		}
	}

	for (size_t ni=0; ni<nbg.size(); ni++)
	{
		int pos = nbg[ni].m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}
	// clear heap
	while (!heapqueue.empty())
	{
		GeoNote nt = heapqueue.top();
		heapqueue.pop();
		int pos = nt.m_id;
		m_pVertex[pos].m_mark = -1;
		m_pVertex[pos].m_LocalGeodesic = -1.0;
		m_pVertex[pos].m_inheap = false;
	}

	//if(geo==0.0) cout<<s<<endl;
	return geo;
}

double CMesh::calGaussianCurvatureIntegration()
{
	const std::vector<double>& vGaussCurv = getGaussCurvature();
	double sum = std::accumulate(vGaussCurv.begin(), vGaussCurv.end(), 0.0);
	return sum;
}

double CMesh::calVolume() const
{
	double vol = 0.0;
	for(int fi = 0; fi < m_nFace; fi++) {
		CFace* face = m_vFaces[fi];
		const Vector3D& v1 = face->getVertex(0)->getPosition();
		const Vector3D& v2 = face->getVertex(1)->getPosition();
		const Vector3D& v3 = face->getVertex(2)->getPosition();

		Vector3D vn;
		crossProduct3D(v1, v2, vn);
		vol += dotProduct3D(vn, v3);
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
		const CVertex* pv = getVertex(*iter);
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

void CMesh::buildPointerVectors()
{
	m_vHalfEdges.reserve(m_nHalfEdge);

	for (int i = 0; i < m_nVertex; ++i)
		m_vVertices.push_back(&m_pVertex[i]);
	for (int i = 0; i < m_nHalfEdge; ++i)
		m_vHalfEdges.push_back(&m_pHalfEdge[i]);
	for (int i = 0; i < m_nFace; ++i)
		m_vFaces.push_back(&m_pFace[i]);

	for (int i = 0; i < m_nVertex; ++i) {	
		CVertex* vert = m_vVertices[i];
		for (int j = 0; j < vert->mOutValence; ++j)
			vert->m_HalfEdges.push_back( m_vHalfEdges[vert->m_piEdge[j]] );
	}

	for (int i = 0; i < m_nHalfEdge; ++i) {
		CHalfEdge* edg = m_vHalfEdges[i];
		edg->m_Vertices[0] = m_vVertices[edg->m_iVertex[0]];
		edg->m_Vertices[1] = m_vVertices[edg->m_iVertex[1]];
		if (edg->m_iTwinEdge >= 0)
			edg->m_eTwin = m_vHalfEdges[edg->m_iTwinEdge];
		edg->m_eNext = m_vHalfEdges[edg->m_iNextEdge];
		edg->m_ePrev = m_vHalfEdges[edg->m_iPrevEdge];
		edg->m_Face = m_vFaces[edg->m_iFace];
	}

	for (int i = 0; i < m_nFace; ++i) {
		CFace* fac = m_vFaces[i];
		for (int j = 0; j < fac->m_nType; ++j) {
			fac->m_Vertices.push_back( m_vVertices[fac->m_piVertex[j]] );
			fac->m_HalfEdges.push_back( m_vHalfEdges[fac->m_piEdge[j]] );
		}
	}

	m_bIsPointerVectorExist = true;
	m_bSeparateStorage = false;
}

void CMesh::buildIndexArrays()
{
	// allocate memory for the half-edge array pointer
	delete []m_pVertex;
	delete []m_pFace;
	delete []m_pHalfEdge;

	m_pVertex = new CVertex[m_nVertex];
	m_pFace = new CFace[m_nFace];
	m_pHalfEdge = new CHalfEdge[m_nHalfEdge];

	for (int i = 0; i < m_nVertex; ++i)
	{
		m_pVertex[i] = *m_vVertices[i];

		if (m_pVertex[i].m_piEdge)
			delete []m_pVertex[i].m_piEdge;

		m_pVertex[i].m_piEdge = new int[m_pVertex[i].mOutValence];
		for (int j = 0; j < m_pVertex[i].mOutValence; ++j)
		{
			m_pVertex[i].m_piEdge[j] = m_vVertices[i]->m_HalfEdges[j]->m_eIndex;
		}
	}
	
	for (int i = 0; i < m_nHalfEdge; ++i)
	{
		m_pHalfEdge[i] = *m_vHalfEdges[i];
		CHalfEdge& he = m_pHalfEdge[i];
		he.m_iFace = m_vHalfEdges[i]->m_Face->m_fIndex;
		he.m_iTwinEdge = (m_vHalfEdges[i]->m_eTwin ? m_vHalfEdges[i]->m_eTwin->m_eIndex : -1);
		he.m_iNextEdge = m_vHalfEdges[i]->m_eNext->m_eIndex;
		he.m_iPrevEdge = m_vHalfEdges[i]->m_ePrev->m_eIndex;
		int i1 = m_vHalfEdges[i]->m_Vertices[0]->m_vIndex;
		int i2 = m_vHalfEdges[i]->m_Vertices[1]->m_vIndex;
		
		assert(i1 >= 0 && i1 < m_nVertex && i2 >= 0 && i2 < m_nVertex);
		
		he.m_iVertex[0] = m_vHalfEdges[i]->m_Vertices[0]->m_vIndex;
		he.m_iVertex[1] = m_vHalfEdges[i]->m_Vertices[1]->m_vIndex;
	}

	for (int i = 0; i < m_nFace; ++i)
	{
		m_pFace[i] = *m_vFaces[i];
		CFace& face = m_pFace[i];
		if (face.m_piEdge) delete []face.m_piEdge;
		if (face.m_piVertex) delete []face.m_piVertex;

		face.m_piEdge = new int[face.m_nType];
		face.m_piVertex = new int[face.m_nType];
		for (int j = 0; j < face.m_nType; ++j)
		{
			int idx = m_vFaces[i]->m_Vertices[j]->m_vIndex;
			
			assert(idx >= 0 && idx < m_nVertex);
			
			face.m_piVertex[j] = m_vFaces[i]->m_Vertices[j]->m_vIndex;
			face.m_piEdge[j] = m_vFaces[i]->m_HalfEdges[j]->m_eIndex;
		}
	}

	m_bIsIndexArrayExist = true;
	m_bSeparateStorage = true;
}

void CMesh::assignElementsIndex()
{
	assert(m_nVertex == m_vVertices.size() && m_nHalfEdge == m_vHalfEdges.size() && m_nFace == m_vFaces.size());
	
	for (int i = 0; i < m_nVertex; ++i)
		m_vVertices[i]->m_vIndex = i;
	
	for (int i = 0; i < m_nHalfEdge; ++i)
		m_vHalfEdges[i]->m_eIndex = i;
	
	for (int i = 0; i < m_nFace; ++i)
		m_vFaces[i]->m_fIndex = i;

	for (int i = 0; i < m_nVertex; ++i) {
		CVertex* pv = m_vVertices[i];
		if (pv->m_piEdge)
			delete []pv->m_piEdge;
		pv->m_piEdge = new int[pv->mOutValence];
		for (int j = 0; j < pv->mOutValence; ++j) {
			pv->m_piEdge[j] = pv->m_HalfEdges[j]->m_eIndex;
		}
	}

	for (int i = 0; i < m_nHalfEdge; ++i) {
		CHalfEdge* he = m_vHalfEdges[i];
		he->m_iVertex[0] = he->m_Vertices[0]->m_vIndex;
		he->m_iVertex[1] = he->m_Vertices[1]->m_vIndex;
		he->m_iTwinEdge = he->m_eTwin ? he->m_eTwin->m_eIndex : -1;
		he->m_iNextEdge = he->m_eNext->m_eIndex;
		he->m_iPrevEdge = he->m_ePrev->m_eIndex;
		he->m_iFace = he->m_Face->m_fIndex;
	}

	for (int i = 0; i < m_nFace; ++i) {
		CFace *pf = m_vFaces[i];
		for (int j = 0; j < pf->m_nType; ++j) {
			pf->m_piVertex[j] = pf->m_Vertices[j]->m_vIndex;
			pf->m_piEdge[j] = pf->m_HalfEdges[j]->m_eIndex;
		}
	}
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
	if (halfEdge->m_eTwin && halfEdge->m_eTwin->m_bIsValid)
	{
		vOppo2 = halfEdge->m_eTwin->m_eNext->m_Vertices[1];
	}

	for (list<CVertex*>::iterator iter1 = v1VList.begin(); iter1 != v1VList.end(); ++iter1)
	{
		for (list<CVertex*>::iterator iter2 = v2VList.begin(); iter2 != v2VList.end(); ++iter2)
		{
			if ( (*iter1 == *iter2) && (*iter1 != vOppo1) && (*iter2 != vOppo2))
				return false;
		}
	}

	return true;
}

void CMesh::scaleAreaToVertexNum()
{
	Vector3D center = calMeshCenter();
	double totalSufaceArea = calSurfaceArea();
	double scale = sqrt( double(m_nVertex) / totalSufaceArea );

	for (int i = 0; i < m_nVertex; ++i) {
		m_vVertices[i]->translateAndScale(-center, scale);
		m_pVertex[i].translateAndScale(-center, scale);
	}
}

void CMesh::gatherStatistics()
{
	// collect/compute statistics
	// 1. face normals
	// 2. vertex normals
	// 3. boundary vertex
	// 4. mesh center
	// 5. bounding box
	// 6. average edge length
	// 7. individual vertex curvature value
	
	calFaceNormals();
	calVertNormals();
	calBoundaryVert();
	
	Vector3D center = calMeshCenter();
	Vector3D boundBox = calBoundingBox(center);

	double edgeLength = 0;
	for (int i = 0; i < m_nHalfEdge; ++i) {
		edgeLength += m_vHalfEdges[i]->getLength();
	}
	edgeLength /= m_nHalfEdge;

	addAttr<double>(edgeLength, StrAttrAvgEdgeLength, AR_UNIFORM, AT_DBL);
	addAttr<Vector3D>(center, StrAttrMeshCenter, AR_UNIFORM, AT_VEC3);
	addAttr<Vector3D>(boundBox, StrAttrMeshBBox, AR_UNIFORM, AT_VEC3);

	calCurvatures();
	calVertMixedAreas();
	findHoles();
}

void CMesh::move( const Vector3D& translation )
{
	for (int i = 0; i < m_nVertex; ++i)
	{
		m_vVertices[i]->translateAndScale(translation, 1.0);
		m_pVertex[i].translateAndScale(translation, 1.0);
	}
}

void CMesh::clearVertexMark()
{
	for (int i = 0; i < m_nVertex; i++) 
		m_pVertex[i].m_mark = -1;
}

void CMesh::scaleEdgeLenToUnit()
{
	Vector3D center(0, 0, 0);
	for (int i = 0; i < m_nVertex; ++i)
		center += m_vVertices[i]->getPosition();
	center /= m_nVertex;

	double length = 0.;
	int edgeNum = 0;
	bool *heVisisted = new bool[m_nHalfEdge];
	for (int i = 0; i < m_nHalfEdge; ++i)
		heVisisted[i] = false;

	for (int i = 0; i < m_nHalfEdge; ++i)
	{
		if (heVisisted[i]) continue;
		
		edgeNum++;
		length += m_vHalfEdges[i]->getLength();

		const CHalfEdge* ptwin = m_vHalfEdges[i]->twinHalfEdge();
		if (ptwin != NULL) heVisisted[ptwin->getIndex()] = true;
	}
	delete []heVisisted;

	length /= edgeNum;
	double scale = 1.0 / length;
	
	for (int i = 0; i < m_nVertex; ++i)
	{
		m_vVertices[i]->translateAndScale(-center, scale);
		m_pVertex[i].translateAndScale(-center, scale);
	}

}

void CMesh::scaleAndTranslate( const Vector3D& center, double scale )
{
	for (int i = 0; i < m_nVertex; ++i)	{
		m_vVertices[i]->translateAndScale(-center, scale);
		m_pVertex[i].translateAndScale(-center, scale);
	}
}

std::vector<int> CMesh::getOriginalVertexIndex() const
{
	vector<int> vret;
	for (CVertex* v : m_vVertices) 
		vret.push_back(v->m_vIndex);
	return vret;
}

void CMesh::vertRingNeighborVerts( int vIndex, int ring, std::vector<int>& nbr, bool inclusive /*= false*/ ) const
{	
	std::set<int> snb;
	vertRingNeighborVerts(vIndex, ring, snb, inclusive);
	
	nbr.clear();
	for (int vn : snb) nbr.push_back(vn);
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
	construct();
}

void CMesh::extractExtrema( const std::vector<double>& vSigVal, int ring, double lowThresh, vector<int>& vFeatures ) const
{
	const int STATE_IDLE = 0;
	const int STATE_MIN	= -1;
	const int STATE_MAX	=  1;

	assert(vSigVal.size() == m_nVertex);
	vFeatures.clear();
	
	int state = STATE_IDLE;
	int state_c = STATE_IDLE;
	vector<int> nb;

	double pz = 1e-5;		//1e-5
	double nz = -1e-5;

	for(int j = 0; j < m_nVertex; j++)		//m_size: size of the mesh
	{
		state = STATE_IDLE;
		if (m_pVertex[j].m_bIsBoundary) 
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

void CMesh::extractExtrema( const std::vector<double>& vSigVal, int ring, std::vector<std::pair<int, int> >& vFeatures, double lowThresh, int avoidBoundary/* = 1*/ ) const
{
	const int STATE_IDLE = 0;
	const int STATE_MIN	= -1;
	const int STATE_MAX	=  1;

	assert(vSigVal.size() == m_nVertex);
	vFeatures.clear();

	int state = STATE_IDLE;
	int state_c = STATE_IDLE;
	vector<int> nb;

	double pz = 0;		//1e-5
	double nz = -0;

	vector<double> sigDetected;

	for(int j = 0; j < m_nVertex; j++)		//m_size: size of the mesh
	{
		if (m_vVertices[j]->m_bIsBoundary) continue;  // ignore boundary vertex
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

bool CMesh::hasBoundary() const
{
	assert(hasAttr(StrAttrBoundaryVertCount));
	const MeshAttr<int> *attrBoundary = getAttr<int>(StrAttrBoundaryVertCount);
	return attrBoundary->attrValue() > 0;
}

double CMesh::calFaceArea( int i ) const
{
	const CFace* f = m_vFaces[i];
	return f->calArea();
}

void CMesh::getVertCoordinateFunction( int dim, std::vector<double>& vCoord ) const
{
	vCoord.resize(m_nVertex);
	switch(dim) {
	case 0:
		for (int i = 0; i < m_nVertex; ++i) vCoord[i] = m_vVertices[i]->getPosition().x; 
		break;
	case 1:
		for (int i = 0; i < m_nVertex; ++i) vCoord[i] = m_vVertices[i]->getPosition().y; 
		break;
	case 2:
		for (int i = 0; i < m_nVertex; ++i) vCoord[i] = m_vVertices[i]->getPosition().z; 
		break;
	default:
		throw std::logic_error("Invalid coordinate dimension!");
	}
}

MeshCoordinates CMesh::getVertCoordinates() const
{
    MeshCoordinates coords(m_nVertex);
    ZGeom::VecNd& vx = coords.getCoordFunc(0);
    ZGeom::VecNd& vy = coords.getCoordFunc(1);
    ZGeom::VecNd& vz = coords.getCoordFunc(2);
    for (int i = 0; i < m_nVertex; ++i) {
        auto vCoord = m_vVertices[i]->getPosition();
        vx[i] = vCoord.x;
        vy[i] = vCoord.y;
        vz[i] = vCoord.z;
    }

    return coords;
}

ZGeom::PointCloud3d CMesh::toPointCloud() const
{
    std::vector<ZGeom::Vec3d> vPoints;
    for (int i = 0; i < m_nVertex; ++i) 
        vPoints.push_back(toVec3d(m_vVertices[i]->getPosition()));    
    return ZGeom::PointCloud3d(vPoints);
}



void CMesh::setVertCoordinates( const MeshCoordinates& coords )
{
	ZGeom::logic_assert(coords.size() == m_nVertex, "Size of coordinates and mesh not compatible!");
	std::vector<double> vx = coords.getCoordFunc(0).toStdVector(),
	                    vy = coords.getCoordFunc(1).toStdVector(),
	                    vz = coords.getCoordFunc(2).toStdVector();
	
	setVertexCoordinates(vx, vy, vz);
}

void CMesh::setVertexCoordinates( const std::vector<double>& vxCoord, const std::vector<double>& vyCoord, const std::vector<double>& vzCoord )
{
	assert(vxCoord.size() == vyCoord.size() && vxCoord.size() == vzCoord.size() && vxCoord.size() == m_nVertex);

	for (int i = 0; i < m_nVertex; ++i)	{
		m_pVertex[i].setPosition(vxCoord[i], vyCoord[i], vzCoord[i]);
		
		if (m_bIsPointerVectorExist)
			m_vVertices[i]->setPosition(vxCoord[i], vyCoord[i], vzCoord[i]);
	}

    // after changing coordinates, vert and face normals need recalculation
    calFaceNormals();
    calVertNormals();
}

void CMesh::setVertexCoordinates(const std::vector<int>& vDeformedIdx, const std::vector<Vector3D>& vNewPos)
{
	if(vDeformedIdx.size() != vNewPos.size())
		throw std::logic_error("Error: CMesh::setVertexCoordinates; incompatible parameters");

	size_t vsize = vDeformedIdx.size();
	for (size_t i = 0; i < vsize; ++i)
	{
		m_pVertex[vDeformedIdx[i]].setPosition(vNewPos[i].x, vNewPos[i].y, vNewPos[i].z);
		
		if (m_bIsPointerVectorExist)
			m_vVertices[vDeformedIdx[i]]->setPosition(vNewPos[i].x, vNewPos[i].y, vNewPos[i].z);
	}
}

double CMesh::getAvgEdgeLength() const
{
	const MeshAttr<double>* attrAvgEdgeLen = getAttr<double>(StrAttrAvgEdgeLength);
	if (attrAvgEdgeLen == NULL) throw std::logic_error("Attribute of average edge length not available!");
	return attrAvgEdgeLen->attrValue();
}

const std::vector<double>& CMesh::getMeanCurvature()
{
	if (!hasAttr(StrAttrVertMeanCurvatures)) this->calCurvatures();		
	return getAttrValue<std::vector<double>>(StrAttrVertMeanCurvatures);
}

const std::vector<double>& CMesh::getMeanCurvature() const
{
	assert(hasAttr(StrAttrVertMeanCurvatures));
	return getAttrValue<std::vector<double>>(StrAttrVertMeanCurvatures);
}

const std::vector<double>& CMesh::getGaussCurvature()
{
	if (!hasAttr(StrAttrVertGaussCurvatures)) calCurvatures();		
	return getAttrValue<std::vector<double>>(StrAttrVertGaussCurvatures);
}

const std::vector<Vector3D>& CMesh::getFaceNormals()
{
	if (!hasAttr(StrAttrFaceNormal)) calFaceNormals();
	const std::vector<Vector3D>& vFaceNormals = getAttrValue<std::vector<Vector3D>>(StrAttrFaceNormal);
	return vFaceNormals;
}

const std::vector<Vector3D>& CMesh::getVertNormals()
{
	if (!hasAttr(StrAttrVertNormal)) calVertNormals();
	return getAttrValue<std::vector<Vector3D>>(StrAttrVertNormal);
}

const std::vector<Vector3D>& CMesh::getVertNormals() const
{
	assert(hasAttr(StrAttrVertNormal));
	return getAttrValue< std::vector<Vector3D>>(StrAttrVertNormal);
}

const std::vector<bool>& CMesh::getVertsOnHole()
{
	if(!hasAttr(StrAttrVertOnHole)) findHoles();
	return getAttrValue<std::vector<bool>>(StrAttrVertOnHole);
}

const std::vector<bool>& CMesh::getVertsOnHole_const() const
{
	assert(hasAttr(StrAttrVertOnHole));
	return getAttrValue<std::vector<bool>>(StrAttrVertOnHole);
}

const std::vector<bool>& CMesh::getVertsOnBoundary()
{
	if (!hasAttr(StrAttrVertOnBoundary)) calBoundaryVert();
	return getAttrValue< std::vector<bool>>(StrAttrVertOnBoundary);
}

double CMesh::calSurfaceArea() const
{
	double totalSufaceArea(0);
	for (int i = 0; i < m_nFace; ++i) 
		totalSufaceArea += calFaceArea(i);
	return totalSufaceArea;
}

Vector3D CMesh::calMeshCenter() const
{
	Vector3D center(0, 0, 0);
	for (int i = 0; i < m_nVertex; ++i)
		center += m_vVertices[i]->getPosition();
	center /= m_nVertex;
	return center;
}

Vector3D CMesh::calBoundingBox( const Vector3D& center ) const
{
	Vector3D boundBox(0.0, 0.0, 0.0);
	for (int i = 0; i < m_nVertex; ++i) {
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

void CMesh::getGraphCSR( std::vector<int>& xadj, std::vector<int>& adjncy ) const
{
	int vertNum = vertCount();
	xadj.resize(vertNum + 1);
	adjncy.clear();
	xadj[0] = 0;
	for (int i = 0; i < vertNum; ++i) {
		std::vector<int> vNbr = getVertNeighborVerts(i, 1, false);
		for (int j : vNbr) adjncy.push_back(j);
		xadj[i+1] = (int)adjncy.size();
	}
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

		std::list<Vector3D> VertexList;	//temporary vertex list
		std::list<int> FaceList;		//temporary face list

		for (int i = 0; i < subVertCount; ++i) {
			VertexList.push_back(this->getVertexPosition(subMappedIdx[i]));
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

		subMesh.m_nVertex = (int)VertexList.size();
		subMesh.m_nFace = (int)FaceList.size() / 3;
		subMesh.m_nHalfEdge = 3 * subMesh.m_nFace;

		subMesh.m_pVertex = new CVertex[subMesh.m_nVertex];
		subMesh.m_pFace = new CFace[subMesh.m_nFace];

		list<Vector3D>::iterator iVertex = VertexList.begin();
		list<int>::iterator iFace = FaceList.begin();
		for(int i = 0; i < subMesh.m_nVertex; i++) {
			subMesh.m_pVertex[i].m_vPosition = *iVertex++;  
			subMesh.m_pVertex[i].m_vIndex = subMesh.m_pVertex[i].m_vid = i;
		}
		for(int i = 0; i < subMesh.m_nFace; i++) {
			subMesh.m_pFace[i].create(3);
			for(int j = 0; j < 3; j++)
				subMesh.m_pFace[i].m_piVertex[j] = *iFace++;
		}
		subMesh.construct();
	}	
}

void CMesh::calVertMixedAreas()
{
	const double pi = ZGeom::PI;
	const int vertCount = this->vertCount();
	std::vector<double> vMixedAreas(vertCount, 0); 

	for (int vIndex = 0; vIndex < vertCount; ++vIndex) 
	{
		CVertex* vi = m_vVertices[vIndex];
		double sum = 0.0;		// sum of attaching corner's angle
		double amix = 0.0;

		for (auto he = vi->m_HalfEdges.begin(); he != vi->m_HalfEdges.end(); ++he) {
			CHalfEdge* e0 = *he;
			CHalfEdge* e1 = e0->m_eNext;
			CHalfEdge* e2 = e1->m_eNext;
			double len0 = e0->getLength();
			double len1 = e1->getLength();
			double len2 = e2->getLength();
			double cota, cotc;
			amix += calAreaMixed(len0, len1, len2, cota, cotc);
		}

		vMixedAreas[vIndex] = amix;
	}

	addAttr<std::vector<double>>(vMixedAreas, StrAttrVertMixedArea, AR_VERTEX, AT_VEC_DBL);
}

const std::vector<double>& CMesh::getVertMixedAreas()
{
	if (!hasAttr(StrAttrVertMixedArea)) calVertMixedAreas();
	return getAttrValue<std::vector<double>>(StrAttrVertMixedArea);
}

const std::vector<double>& CMesh::getVertMixedAreas() const
{
	assert(hasAttr(StrAttrVertMixedArea));
	return getAttrValue<std::vector<double>>(StrAttrVertMixedArea);
}

std::vector<ZGeom::Vec3d> CMesh::getAllVertPositions() const
{
    using ZGeom::Vec3d;
    vector<Vec3d> results(m_nVertex);
    for (int i = 0; i < m_nVertex; ++i) {
        results[i] = Vec3d(m_vVertices[i]->getPosition());
    }
    return results;
}
