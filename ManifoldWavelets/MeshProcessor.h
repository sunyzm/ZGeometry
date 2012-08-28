#pragma once
#include <mesh/Mesh.h>
#include <engine.h>
#include <mesh/Laplacian.h>
#include <string>
#include <vector>

class MeshProcessor
{
public:
	MeshProcessor(void);
	~MeshProcessor(void);
	void init(CMesh* tm, Engine* e);
	void decomposeLaplacian(int nEigFunc);
	void readMHB(std::string path);
	void writeMHB(std::string path);
	void normalizeFrom(const std::vector<double>& vFrom);

	CMesh* mesh;
	Engine *ep;
	ManifoldHarmonics mhb;
	std::vector<double> vDisplaySignature;
	int pRef;
	bool isMHBuilt;
	int m_size;
};

