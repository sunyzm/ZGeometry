#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <ZGeom/ZGeom.h>
#include <ZGeom/util.h>

//////////////////////////////////////////////////////////////////////////
// format: zmesh_info [mesh_file_name] 
//
int main(int argc, char *argv[])
{
	if (argc < 2) {
		std::cerr << "Lack argument!" << std::endl;
		std::exit(-1);
	}
	
	/// preparation of option and parameters
	std::string meshFilename = argv[1];
	
	if (!fileExist(meshFilename)) {
		std::cerr << "Mesh file not existent!" << std::endl;
		std::exit(-1);
	}

	/// load mesh
	CMesh oriMesh;
	oriMesh.load(meshFilename);
	oriMesh.gatherStatistics();

	std::string meshName = oriMesh.getMeshName();
	int totalVertCount = oriMesh.vertCount();
	auto center = oriMesh.getCenter();
	auto vVertesOnBoundary = oriMesh.getVertsOnBoundary();
	int countVertsOnBoundary = (int)std::count_if(vVertesOnBoundary.begin(), vVertesOnBoundary.end(), [](bool bv){ return bv; });

	std::cout << "** info of mesh " << meshName << " **\n";
	std::cout << "- vertex count: " << totalVertCount << '\n';
	std::cout << "- face count: " << oriMesh.faceCount() << '\n';
	std::cout << "- center: (" << center.x << ", " << center.y << ", " << center.z << ")\n";
	std::cout << "- boundary vert num: " << countVertsOnBoundary << '\n';

#ifdef _DEBUG
	std::system("PAUSE");
#endif
	return 0;
}