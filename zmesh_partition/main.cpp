#include <cstdlib>
#include <iostream>
#include <fstream>
#include "ZGeom/ZGeom.h"

//////////////////////////////////////////////////////////////////////////
// format: zmesh_partition [mesh_file_name] -n/-m [partitionCount/maxPatchSize]
//
int main(int argc, char *argv[])
{
	if (argc < 4) {
		std::cerr << "Lack argument!\n";
		exit(-1);
	}
	else if (argc > 4) {
		std::cerr << "Excessive arguments ignored!\n";
	}



	return 0;
}