#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <ZGeom/ZGeom.h>
#include <ZGeom/Laplacian.h>
#include <ZGeom/util.h>
#include <ZGeom/MatlabEngineWrapper.h>

//////////////////////////////////////////////////////////////////////////
// format: zmesh_operator [mesh_file_name] 
// currently, just output three types of Laplacian operator in .mat Matlab file
int main(int argc, char *argv[])
{
    using namespace std;
    using namespace ZGeom;

    if (argc < 2) {
        std::cerr << "Lack argument!" << std::endl;
        std::exit(-1);
    }

    string meshfile = argv[1];
    if (!fileExist(meshfile)) {
        std::cerr << "Mesh file not existent!" << std::endl;
        std::exit(-1);
    }

    CMesh mesh;
    mesh.load(meshfile);
    Laplacian umbrella, goemUmbrella, cotformula;
    umbrella.constructUmbrella(&mesh);
    cotformula.constructCotFormula(&mesh);
    MatlabEngineWrapper eng;
    eng.open();
    eng.addSparseMat(umbrella.getLS(), "mat_umbrella");
    eng.addSparseMat(cotformula.getW(), "mat_weight");
    eng.addSparseMat(cotformula.getLS(), "mat_cot");
    eng.close();
    exit(0);
}