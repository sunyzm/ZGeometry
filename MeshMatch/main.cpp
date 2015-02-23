#include <ZGeom/Mesh.h>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
	if (argc < 3) {
		std::cerr << "Arguments should include 2 mesh files!" << std::endl;
		exit(-1);
	}
	CMesh mesh1, mesh2;
	try {
		mesh1.load(argv[1]);
		mesh2.load(argv[2]);
	}
	catch (std::exception& e) {
		std::cerr << e.what() << '\n';
		exit(-2);
	}
	
	char *outfile;
	if (argc >= 4) outfile = argv[3];
	else outfile = "correspondence.txt";

	std::ofstream fout(outfile);

	std::vector<std::pair<int,int>> matchPairs;

	int meshSize1 = mesh1.vertCount(), meshSize2 = mesh2.vertCount();
	for (int i = 0; i < meshSize1; ++i) {
		const ZGeom::Vec3d& v1 = mesh1.getVertex(i)->pos();
		for (int j = 0; j < meshSize2; ++j) {
			const ZGeom::Vec3d& v2 = mesh2.getVertex(j)->pos();
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