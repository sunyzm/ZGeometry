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
#include <ZGeom/spectral_analysis.h>


// synopsis:
// 1a) zmesh eigen f:[input mesh filename] [eigen-num] [output path]
// 1b) zmesh eigen l:[input mesh filelist] [eigen-num] [output path]
// 2a) zmesh feature f:[input mesh filename] [eigen-path] [output path]
// 2b) zmesh feature l:[input mesh fileliest] [eigen-path] [output path]
// 3   zmesh compress f:[mesh filename]

using namespace std;
using namespace ZGeom;
using namespace std::tr2::sys;

const std::set<std::string> spectrum_commands{ "spectrum", "spectral", "eigen" };
const std::set<std::string> feature_commands{ "feature" };

void computeMeshesSpectrum(const std::vector<std::string>& mesh_filenames, int eigen_num, std::string output_path);
void computeMeshesFeaturePoints(const std::vector<std::string>& mesh_filenames, std::string eigen_path, std::string output_path);

void compressMesh(const std::vector<std::string>& mesh_filenames);

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cerr << "Lack argument!" << std::endl;
        exit(-1);
    }

    string files_str(argv[2]);

    vector<string> filenames;
    if (files_str.substr(0, 2) == "l:") {
        filenames = readFileList(files_str.substr(2));
    }
    else if (files_str.substr(0, 2) == "f:") {
        filenames = vector<string>{files_str.substr(2)};
    }

    string task_string(argv[1]);
    if (spectrum_commands.find(task_string) != spectrum_commands.end()) {
        if (argc < 5) {
            std::cerr << "Lack argument!" << std::endl;
            exit(-1);
        }
        string eigen_num_str(argv[3]);
        string out_path(argv[4]);

        computeMeshesSpectrum(filenames, atoi(argv[3]), out_path);
    }
    else if (feature_commands.find(task_string) != feature_commands.end()) {
        string eigen_path(argv[3]);
        string out_path(argv[4]);
        computeMeshesFeaturePoints(filenames, eigen_path, out_path);
    }
    else if (task_string == "compress") {
        compressMesh(filenames);
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

            string mesh_name = mesh.getMeshName();
            string eigen_file_path = mesh_name + ".eigen.mat";            
            es.saveToMat(&eng, eigen_file_path);
            cout << eigen_file_path << " saved!" << endl;
        }
        catch (runtime_error* e){
            cerr << "Fail to process " + mesh_file << ": " << e->what() << endl;
        }
        catch (exception* e) {
            cerr << "Unknown error in processing " + mesh_file << endl;
        }
    }

    cout << "\n";
    timer.stopTimer("Total time: ");
}

void computeMeshesFeaturePoints(const std::vector<std::string>& mesh_filenames, std::string eigen_path, std::string output_path)
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

            string mesh_name = mesh.getMeshName();
            string eigen_file_path = eigen_path + "/" + mesh_name + ".eigen.mat";

            EigenSystem es;
            es.loadFromMat(&eng, eigen_file_path);

            double t_min = 4 * std::log(10.0) / es.getAllEigVals().back(), t_max = 4 * std::log(10.0) / es.getEigVal(1);
            int nScales = 4;
            const double tMultiplier = std::pow(t_max / t_min, 1.0 / double(nScales - 1));
            std::vector<double> hks_scales(nScales);
            for (int s = 0; s < nScales; ++s) {
                hks_scales[s] = t_min * std::pow(tMultiplier, s);
            }

            set<int> all_features;
            for (int s = 0; s < nScales; ++s) {
                vector<double> values = calHeatKernelSignature(es, hks_scales[s]);
                vector<int> features = extractMeshExtrema(mesh, values, 2);
                for (int fi : features) all_features.insert(fi);
            }

            string out_file_path = output_path + "/" + mesh_name + ".feature";
            ofstream ofs(output_path.c_str());

        }
        catch (runtime_error* e){
            cerr << "Fail to process " + mesh_file << ": " << e->what() << endl;
        }
    } 

    cout << "\n";
    timer.stopTimer("Total time: ");
}


void compressMesh(const std::vector<std::string>& mesh_filenames) {
    string mesh_file = mesh_filenames.front();
    CMesh mesh;
    mesh.load(mesh_file);
    string mesh_name = mesh.getMeshName();    
    ofstream ofs(mesh_name + ".compress.log");

    // what to do:
    // 1. partition mesh with each mesh of size 300, 500, 1000, 1500, 2000
    // 2. For each submesh, do eigendecomposition
    // 3. Compute SGW dictionary
    // 4. S-OMP approximation
    // 5. Output total approximation error

    // 
    // what to output:
    // 1. eigendecomposition time
    // 2. S-OMP time
    // 3. Approximation error

    vector<int> max_sizes{ 300, 600, 1000, 1500, 2000 };
    for (int max_size : max_sizes) {
        int nPart = std::round(mesh.vertCount() / max_size);
        vector<int> part_idx = MetisMeshPartition(&mesh, nPart);
        // TODO

    }

}
