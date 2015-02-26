#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <metis.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/util.h>

std::vector<int> MetisMeshPartition(const CMesh* mesh, int nPart)
{
	int vertCount = mesh->vertCount();
	std::vector<int> vPart(vertCount);
	int ncon = 1;
	std::vector<int> xadj, adjncy;
    ZGeom::getMeshGraphCSR(*mesh, xadj, adjncy);
	int objval;
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_CONTIG] = 1;
	options[METIS_OPTION_NUMBERING] = 0;

	METIS_PartGraphKway(&vertCount, &ncon, &xadj[0], &adjncy[0], NULL, NULL, NULL, &nPart, NULL, NULL, NULL, &objval, &vPart[0]);

	return vPart;
}

//////////////////////////////////////////////////////////////////////////
// usage: zmesh_partition [mesh-filename] -n/-m [partitionCount/maxPatchSize] [output-path]
//
int main(int argc, char *argv[])
{
	if (argc == 2 && argv[1] == std::string("--help")) {
		std::cout << "usage: zmesh_partition [mesh-filename] -n/-m [partitionCount/maxPatchSize] [output-path]" << '\n';
		exit(0);
	}

	if (argc < 4) {
		std::cerr << "Lack " << 4 -argc << " arguments!" << std::endl;
		exit(-1);
	}
	else if (argc > 5) {
		std::cerr << "Excessive arguments ignored!" << std::endl;
	}

	/// preparation of option and parameters
	std::string meshFilename = argv[1];
	std::string segOption = argv[2];
	std::string segPara = argv[3];
	std::string outputPath = "./";
	if (argc == 5) {
		outputPath = argv[4];
		if (outputPath.back() != '/') outputPath.push_back('/');	// add trailing slash if missing
	}

	if (!fileExist(meshFilename)) {
		std::cerr << "Mesh file not existent!" << std::endl;
		std::exit(-1);
	}
	bool segNumSpecified = false;
	if (segOption == "-n") segNumSpecified = true;
	else if (segOption == "-m") segNumSpecified = false;
	else {
		std::cerr << "Unrecognized option '" << segOption << " '" << std::endl;
		std::exit(-1);
	}
	int nPara = std::stoi(segPara);

	/// load mesh
	CMesh oriMesh;
	oriMesh.load(meshFilename);
    ZGeom::gatherMeshStatistics(oriMesh);

	std::string meshName = oriMesh.getMeshName();
	int totalVertCount = oriMesh.vertCount();
	int nPart;
	if (segNumSpecified) nPart = nPara;
	else nPart = totalVertCount / nPara + 1;
	if (nPart <= 0) {
		std::cerr << "Number of segmentations must be greater than 0" << std::endl;
		std::exit(-1);
	}

	std::vector<CMesh*> vSubMeshes;
	std::vector<std::vector<int>*> vMappedIdx;

	for (int i = 0; i < nPart; ++i) {
		vSubMeshes.push_back(new CMesh());
		vMappedIdx.push_back(new std::vector<int>());
	}
	std::vector<int> vPartIdx = MetisMeshPartition(&oriMesh, nPart);
	for (int vIdx = 0; vIdx < totalVertCount; ++vIdx) {
		vMappedIdx[vPartIdx[vIdx]]->push_back(vIdx);
	}

	double paritionTime = time_call_sec([&]() {
		oriMesh.partitionToSubMeshes(vMappedIdx, vSubMeshes);
	});
	std::cout << "Partition finished! Time consumed: " << paritionTime << "s\n";

	char numBuf[10];
	for (int i = 0; i < nPart; ++i) {
		_itoa_s(i + 1, numBuf, 10);
		std::string subMeshName = meshName + ".sub" + std::string(numBuf);
		vSubMeshes[i]->setMeshName(subMeshName);
	}

	std::cout << "- number of partitions: " << nPart << '\n';
	std::cout << "- sub-mesh sizes: ";
	for (int i = 0; i < nPart; ++i) {
		std::cout << vSubMeshes[i]->getMeshName() << ": " << vSubMeshes[i]->vertCount() << " | ";
		vSubMeshes[i]->save(outputPath + vSubMeshes[i]->getMeshName() + ".obj");
	}
	std::cout << '\n';

	for (auto p : vSubMeshes) delete p;
	for (auto p : vMappedIdx) delete p;

#ifdef _DEBUG
	std::system("PAUSE");
#endif
	exit(0);
}