#pragma once
#include <mesh/Mesh.h>
#include <engine.h>
#include <mesh/Laplacian.h>
#include <string>

class MeshProcessor
{
public:
	MeshProcessor(void);
	~MeshProcessor(void);
	void decomposeLaplacian(int nEigFunc);
	void readMHB(std::string path);
	void writeMHB(std::string path);

	CMesh* mesh;
	Engine *ep;
	ManifoldHarmonics mhb;

	bool isMHBuilt;
	int size;
};

