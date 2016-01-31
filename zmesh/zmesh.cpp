#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <regex>
#include <filesystem>
#include <set>
#include <ZGeom/ZGeom.h>
#include <ZGeom/Laplacian.h>
#include <ZGeom/util.h>
#include <ZGeom/MatlabEngineWrapper.h>

// synopsis:
// 1a) zmesh eigen f:[input mesh filename] [eigen-num] [output path]
// 1b) zmesh eigen l:[input mesh filelist] [eigen-num] [output path]

using namespace std;
using namespace ZGeom;
using namespace std::tr2::sys;

const std::set<std::string> spectrum_commands{ "spectrum", "spectral", "eigen" };

void computeMeshesSpectrum(const std::vector<std::string>& mesh_filenames, int eigen_num, std::string output_path);

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cerr << "Lack argument!" << std::endl;
        exit(-1);
    }

    string task_string(argv[1]);
    if (spectrum_commands.find(task_string) != spectrum_commands.end()) {
        if (argc < 5) {
            std::cerr << "Lack argument!" << std::endl;
            exit(-1);
        }

        try
        {
            string files_str(argv[2]);
            string eigen_num_str(argv[3]);
            string out_path(argv[4]);
            
            vector<string> filenames;
            if (files_str.substr(0, 2) == "l:") {
                filenames = readFileList(files_str.substr(2));
            }
            else if (files_str.substr(0, 2) == "f:") {
                filenames = vector<string>{files_str.substr(2)};
            }

            computeMeshesSpectrum(filenames, atoi(argv[3]), out_path);
        }
        catch (runtime_error* e)
        {
            cerr << e->what() << endl;
        }
    }

#ifdef _DEBUG
    system("PAUSE");
#endif
    return 0;
}

void computeMeshesSpectrum(const std::vector<std::string>& mesh_filenames, int eigen_num, std::string output_path)
{
    ZGeom::MatlabEngineWrapper eng;
    eng.open();
    auto pwd = initial_path<path>();
    eng.eval("cd " + pwd.string());
    eng.eval("mkdir " + output_path);
    eng.eval("cd " + output_path);
    cout << "# Total meshes: " << mesh_filenames.size() << "\n" << endl;

    CStopWatch timer;
    timer.startTimer();

    int count = 0;
    for (string mesh_file : mesh_filenames)
    {
        count++;
        cout << count << ": " << mesh_file << endl;
        try {
            if (!fileExist(mesh_file)) {
                throw runtime_error("File not exist!");
            }
            CMesh mesh;
            mesh.load(mesh_file);
            Laplacian cot_form_lap;
            cot_form_lap.constructCotFormula(&mesh);
            EigenSystem es;
            cot_form_lap.decompose(eigen_num, &eng, es, true);
            eng.addDoubleScalar(es.eigVecCount(), "eigen_num");
            eng.addDoubleScalar(es.eigVecSize(), "mesh_size");
            eng.addColVec(VecNd(es.getAllEigVals()), "eigenvals");
            eng.addDenseMat(es.getEigenMat(), "eigenvecs", true);
            eng.addSparseMat(es.getInducingMat(), "inducing_mat");

            string mesh_name = mesh.getMeshName();
            string eigen_file_path = mesh_name + ".eigen.mat";
            eng.eval("save " + eigen_file_path + " eigen_num mesh_size eigenvals eigenvecs inducing_mat");
            cout << eigen_file_path << " saved!" << endl;
        }
        catch (runtime_error* e){
            cerr << "Fail to process " + mesh_file << ": " << e->what() << endl;
        }
        catch (exception& e) {
            cerr << "Fail to process " + mesh_file << endl;
        }
    }

    cout << "\n";
    timer.stopTimer("Total time: ");
}