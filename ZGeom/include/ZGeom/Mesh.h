#ifndef ZGEOM_MESH_H
#define ZGEOM_MESH_H
#include <vector>
#include <string>
#include "Vec3.h"

namespace ZGeom
{

class Vertex;
class HalfEdge;
class Face;
class Mesh;

class Vertex
{
public:
	friend Mesh;

	Vertex() : mIndex(-1), mPos(0,0,0) {}
	Vertex(int index, const Vec3d& vec) : mIndex(index), mPos(vec) {}
	Vertex(const Vertex& v) : mIndex(v.mIndex), mPos(v.mPos) {}
	Vertex& operator = (const Vertex& v) { mIndex = v.mIndex; mPos = v.mPos; mHalfEdges.clear(); }
	int getIndex() const { return mIndex; }
	void setIndex(int index) { mIndex = index; }
	Vec3d getPos() const { return mPos; }
	void setPos(const Vec3d& vec) { mPos = vec; }
	const std::vector<HalfEdge*>& getHalfEdges() const { return mHalfEdges; }

private:
	int mIndex;
	Vec3d mPos;
	std::vector<HalfEdge*> mHalfEdges;
};

class HalfEdge
{
public:
	friend Mesh;

	HalfEdge(int index = -1) : mIndex(index) { mVerts[0] = mVerts[1] = nullptr; mFace = nullptr; }
	int getIndex() const { return mIndex; }
	void setIndex(int index) {	mIndex = index; }
	const Vertex* getVert(int i) const { 
		if (i != 0 && i != 1) throw std::runtime_error("Invalid subscript for ZGeom::HalfEdge::getVert");
		return mVerts[i]; 
	}
	void setVert(int i, Vertex* vert) { 
		if (i != 0 && i != 1) throw std::runtime_error("Invalid subscript for ZGeom::HalfEdge::getVert");
		mVerts[i] = vert; 
	}
	const Face* getFace() const { return mFace; }
	void setFace(Face* face) { mFace = face; }

private:
	int mIndex;
	Vertex* mVerts[2];
	Face* mFace;
};

class Face
{
public:
	friend Mesh;

	int getIndex() const { return mIndex; }
	void setIndex(int index) { mIndex = index; }
	const std::vector<Vertex*>& getVerts() const { return mVerts; }
	std::vector<Vertex*>& getVerts() { return mVerts; }
	const Vertex* getVert(uint k) const {
		if (k >= mVerts.size()) 
			throw std::runtime_error("Invalid subscript for ZGeom::Face::mVerts");
		return mVerts[k];
	}
	Vertex* getVert(uint k) {
		if (k >= mVerts.size()) 
			throw std::runtime_error("Invalid subscript for ZGeom::Face::mVerts");
		return mVerts[k];
	}

private:
	int mIndex;
	std::vector<Vertex*> mVerts;
};

class Mesh
{
public:
	Mesh();
	void cloneWithoutAttr(const Mesh& mesh);

	const std::string& getName() const { return mName;}
	void setName(const std::string& name) { mName = name; }

	const std::vector<Vertex>& getVerts() const { return mVerts; }
	std::vector<Vertex>& getVerts() { return mVerts; }
	const Vertex& getVert(uint vIndex) const;
	Vertex& getVert(uint vIndex);

	const std::vector<HalfEdge>& getHalfEdges() const { return mHalfEdges; }
	std::vector<HalfEdge>& getHalfEdges() { return mHalfEdges; }
	const HalfEdge& getHalfEdge(uint eIndex) const;
	HalfEdge& getHalfEdge(uint eIndex);

	const std::vector<Face>& getFaces() const { return mFaces; }
	std::vector<Face>& getFaces() { return mFaces; }
	const Face& getFace(uint fIndex) const;
	Face& getFace(uint fIndex);



private:
	std::string mName;
	std::vector<Vertex> mVerts;
	std::vector<HalfEdge> mHalfEdges;
	std::vector<Face> mFaces;
};

} // end of namespace

#endif