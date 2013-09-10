#include <ZMesh/Mesh.h>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
	if (argc < 3) {
		std::cerr << "Arguments should include 2 mesh files!" << std::endl;
		exit(-1);
	}
	CMesh mesh1, mesh2;
	if (!mesh1.Load(argv[1]) || !mesh2.Load(argv[2])) {
		std::cerr << "Fail to load mesh files, please check path or file contents!" << std::endl;
		exit(-2);
	}

	char *outfile;
	if (argc >= 4) outfile = argv[3];
	else outfile = "correspondence.txt";

	std::ofstream fout(outfile);

	std::vector<std::pair<int,int>> matchPairs;

	int meshSize1 = mesh1.getMeshSize(), meshSize2 = mesh2.getMeshSize();
	for (int i = 0; i < meshSize1; ++i) {
		const Vector3D& v1 = mesh1.getVertex(i)->getPosition();
		for (int j = 0; j < meshSize2; ++j) {
			const Vector3D& v2 = mesh2.getVertex(j)->getPosition();
			if (v1.equals(v2)) {
				matchPairs.push_back(std::make_pair(i,j));
				break;
			}
		}
	}

	fout << matchPairs.size() << std::endl;
	for (auto iter = matchPairs.begin(); iter != matchPairs.end(); ++iter) {
		fout << iter->first << ' ' << iter->second << std::endl;
	}
	fout.close();

	return 0;
}