#ifndef MESH_PYRAMID_H
#define MESH_PYRAMID_H
#include "Mesh.h"
#include "Quadric.h"
#include <vector>
#include <list>
#include <map>
#include <iostream>

typedef std::map<int, std::list<int> > VerticesCover;

class SparseMatrix
{
public:
	std::vector<int> m_ii, m_jj;
	std::vector<double> m_ss;
	void dump(std::string filename) const;
};

class MeshPyramid
{
public:
	MeshPyramid();
	MeshPyramid(CMesh* mesh);
	~MeshPyramid();

	void			setInitialMesh(CMesh* mesh);
	void			setLevel(int level, double contract_ratio = 0.5) { m_levels = level; m_contract_ratio = contract_ratio; }
	int				numOfLevels() const { return m_levels; }
	void			construct(std::ostream& ostr);
	void			clear();
	CMesh*			getMesh(int level) const;
	std::list<int>	getCoveredVertexList(int level, int idx) const;	//return index list of vertices that collapsed to a specific vertex on given level 
	void			dumpVertexValence(int level, std::string filename);
	int				**m_Id2IndexMap;

private:
	struct VertexPair
	{
		CVertex* pV[2];
		CHalfEdge* pHE;
		int vToKeep;		//0 or 1
		double cost;
		CHalfEdge* correspondingHalfEdge();
	};

	struct MeshLevel
	{
		CMesh* mesh;
		VerticesCover verticesCover, verticesCoverFull;
	};

	int			m_levels;
	int		    m_nVertices;		// number of vertices
	int			m_nFace;	 		// number of faces
	double      m_contract_ratio;
	CMesh*					originalMesh;
	std::vector<MeshLevel>	m_vMeshes;
	std::vector<Quadric>	m_FaceQuadrics;
	std::vector<Quadric>	m_VertexQuadrics;
};

#endif