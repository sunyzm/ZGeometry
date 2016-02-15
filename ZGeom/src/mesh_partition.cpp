#include "mesh_partition.h"
#include <queue>
#include <map>
#include <metis.h>
#include "SparseMatrix.h"
#include "Laplacian.h"
#include "geometry_processing.h"

using namespace std;

std::vector<int> ZGeom::MetisMeshPartition(const CMesh* mesh, int nPart)
{
    int vertCount = mesh->vertCount();
    std::vector<int> vPart(vertCount);
    int ncon = 1;
    std::vector<int> xadj, adjncy;
    getMeshGraphCSR(*mesh, xadj, adjncy);
    int objval;
    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_NUMBERING] = 0;

    METIS_PartGraphKway(&vertCount, &ncon, &xadj[0], &adjncy[0], NULL, NULL, NULL, &nPart, NULL, NULL, NULL, &objval, &vPart[0]);

    return vPart;
}

std::vector<int> ZGeom::FiedlerMeshPartition(CMesh* original_mesh, int nParts)
{
    MatlabEngineWrapper eng;
    eng.open();

    struct SubMeshStruct {
        CMesh* submesh;
        map<int, int> map2Original;
        int level;
        SubMeshStruct(CMesh* m, const map<int, int>& vert_map, int l) : submesh(m), map2Original(vert_map), level(l) {}
    };

    int levels = lround(log(double(nParts)) / log(2.));
    vector<int> mesh_vert_partitions(original_mesh->vertCount());
    int cur_par = 0;

    queue<SubMeshStruct> submesh_queue;     

    map<int, int> original_map;
    for (int i = 0; i < original_mesh->vertCount(); ++i)
        original_map.insert(make_pair(i, i));
    submesh_queue.emplace(original_mesh, original_map, 0);

    while (!submesh_queue.empty()) {
        SubMeshStruct sms = submesh_queue.front();
        if (sms.level == levels) {
            for (auto& item : sms.map2Original) {
                mesh_vert_partitions[item.second] = cur_par;
            }
            cur_par++;
        }
        else {
            Laplacian cot_formula_lap;
            cot_formula_lap.constructCotFormula(sms.submesh);
            EigenSystem es;
            cot_formula_lap.decompose(2, &eng, es, true);
            vector<double> fiedler_vec = es.getEigVec(1).toStdVector();
            vector<int> sub_verts_1, sub_verts_2;
            for (int i = 0; i < fiedler_vec.size(); ++i) {
                if (fiedler_vec[i] < 0) sub_verts_1.push_back(i);
                else sub_verts_2.push_back(i);
            }
            CMesh* submesh1 = new CMesh, *submesh2 = new CMesh;
            sms.submesh->getSubMesh(sub_verts_1, "", *submesh1);
            sms.submesh->getSubMesh(sub_verts_2, "", *submesh2);

            map<int, int>& mapToOriginal = sms.map2Original;
            map<int, int> toOriginal_1, toOriginal_2;
            for (int i = 0; i < (int)sub_verts_1.size(); ++i) {
                toOriginal_1.insert(make_pair(i, mapToOriginal[sub_verts_1[i]]));
            }
            for (int i = 0; i < (int)sub_verts_2.size(); ++i) {
                toOriginal_2.insert(make_pair(i, mapToOriginal[sub_verts_2[i]]));
            }

            submesh_queue.emplace(submesh1, toOriginal_1, sms.level + 1);
            submesh_queue.emplace(submesh2, toOriginal_2, sms.level + 1);
        }        

        if (sms.level > 0) delete sms.submesh;
        submesh_queue.pop();
    }

    return mesh_vert_partitions;
}

