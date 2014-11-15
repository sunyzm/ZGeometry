#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <regex>
#include <filesystem>
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
    using namespace std::tr2::sys;

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
    int N = mesh.vertCount();
    mesh.scaleAndTranslate(-mesh.calMeshCenter(), 1.0);
    Laplacian umbrella, goemUmbrella, cotformula;
    umbrella.constructUmbrella(&mesh);
    cotformula.constructCotFormula(&mesh);
    MatlabEngineWrapper eng;
    eng.open();
    auto pwd = initial_path<path>();
    eng.eval("cd " + pwd.string()); 

    eng.addSparseMat(umbrella.getLS(), "mat_umbrella");
    eng.addSparseMat(cotformula.getW(), "mat_weight");
    eng.addSparseMat(cotformula.getLS(), "mat_cot");

    auto coord = mesh.getVertCoordinates();
    DenseMatrixd matCoord(N, 3);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < 3; ++j)
            matCoord(i, j) = coord[i][j];    
    mesh.calVertNormals();
    auto normals = mesh.getVertNormals();
    DenseMatrixd matNormal(N, 3);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < 3; ++j)
            matNormal(i, j) = normals[i][j];
    eng.addDenseMat(matCoord, "mat_coord");   
    eng.addDenseMat(matNormal, "mat_normal");

    eng.eval("save laplace.mat mat_umbrella mat_weight mat_cot mat_coord mat_normal");
    eng.close();
    cout << "Laplace mat file saved!" << endl;
    exit(0);
}