#include "Mesh2.h"
#include "zassert.h"

namespace ZGeom
{

const Vertex& Mesh::getVert( uint vIndex ) const
{
	runtime_assert(vIndex < mVerts.size(), "invalid subscript for ZGeom::Mesh::getVert");
	return mVerts[vIndex];
}

Vertex& Mesh::getVert( uint vIndex )
{
	runtime_assert(vIndex < mVerts.size(), "invalid subscript for ZGeom::Mesh::getVert");
	return mVerts[vIndex];
}

const HalfEdge& Mesh::getHalfEdge( uint eIndex ) const
{
	runtime_assert(eIndex < mHalfEdges.size(), "invalid subscript for ZGeom::Mesh::getHalfEdge");
	return mHalfEdges[eIndex];
}

HalfEdge& Mesh::getHalfEdge(uint eIndex) 
{
	runtime_assert(eIndex < mHalfEdges.size(), "invalid subscript for ZGeom::Mesh::getHalfEdge");
	return mHalfEdges[eIndex];
}

const Face& Mesh::getFace( uint fIndex ) const
{
	runtime_assert(fIndex < mFaces.size(), "invalid subscript for ZGeom::Mesh::getFace");
	return mFaces[fIndex];
};

Face& Mesh::getFace( uint fIndex )
{
	runtime_assert(fIndex < mFaces.size(), "invalid subscript for ZGeom::Mesh::getFace");
	return mFaces[fIndex];
};

} // end of namespace