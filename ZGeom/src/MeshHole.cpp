#include "MeshHole.h"
#include <ctime>
#include <random>
#include <set>
using namespace std;

MeshHole autoGenerateHole(const CMesh& mesh, int seedVert, int holeSize)
{
    MeshHole result;

    set<int> vertInHole;
    set<int> faceInHole;
    set<int> boundaryVerts;
    default_random_engine generator((unsigned int)time(NULL));

    vertInHole.insert(seedVert);
    for (auto f : mesh.vert(seedVert)->getAdjacentFaces()) {
        faceInHole.insert(f->getFaceIndex());    
        for (int k = 0; k < 3; ++k) {
            int faceVertIdx = f->getVertexIndex(k);
            if (vertInHole.find(faceVertIdx) == vertInHole.end())
                boundaryVerts.insert(faceVertIdx);
        }    
    }
    
    while (vertInHole.size() < holeSize) {
        vector<int> vecBoundaryVerts{ boundaryVerts.begin(), boundaryVerts.end() };
        uniform_int_distribution<int> distr1(0, (int)vecBoundaryVerts.size() - 1);
        int newHoleVert = vecBoundaryVerts[distr1(generator)];
        vertInHole.insert(newHoleVert);
        for (auto f : mesh.vert(newHoleVert)->getAdjacentFaces()) {
            if (faceInHole.find(f->getFaceIndex()) != faceInHole.end())
                continue;       // face already considered     
            faceInHole.insert(f->getFaceIndex());
            for (int k = 0; k < 3; ++k) {
                int faceVertIdx = f->getVertexIndex(k);
                if (vertInHole.find(faceVertIdx) == vertInHole.end())
                    boundaryVerts.insert(faceVertIdx);
            }
        }
        boundaryVerts.erase(newHoleVert);
        for (auto iter = boundaryVerts.begin(); iter != boundaryVerts.end();) {
            int vIdx = *iter;
            bool encompassed = true;
            for (auto f : mesh.vert(vIdx)->getAdjacentFaces()) {
                if (faceInHole.find(f->getFaceIndex()) == faceInHole.end()) {
                    // found a face not in hole, so vIdx is not in hole yet
                    encompassed = false;
                    break;
                }
            }
            if (encompassed) {
                iter = boundaryVerts.erase(iter);
                vertInHole.insert(vIdx);
            }
            else iter++;
        }
    }

    result.mHoleFaces = vector < int > {faceInHole.begin(), faceInHole.end()};
    result.mHoleVerts = vector < int > {vertInHole.begin(), vertInHole.end()};
    result.mHoleBoundaryVerts = vector < int > {boundaryVerts.begin(), boundaryVerts.end()};
    return result;
}

MeshHole autoGenerateHole(const CMesh& mesh, const std::vector<int>& seedVerts, int totalSize)
{
    assert(totalSize >= seedVerts.size());
    MeshHole result;

    set<int> vertInHole;
    set<int> faceInHole;
    set<int> boundaryVerts;
    default_random_engine generator((unsigned int)time(NULL));

    vertInHole = set < int > {seedVerts.begin(), seedVerts.end()};
    for (int vIdx : vertInHole) {
        for (auto f : mesh.vert(vIdx)->getAdjacentFaces()) {
            faceInHole.insert(f->getFaceIndex());
            for (int k = 0; k < 3; ++k) {
                int faceVertIdx = f->getVertexIndex(k);
                if (vertInHole.find(faceVertIdx) == vertInHole.end())
                    boundaryVerts.insert(faceVertIdx);
            }
        }
    }
    

    while (vertInHole.size() < totalSize) {
        uniform_int_distribution<int> distr1(0, (int)boundaryVerts.size() - 1);
        std::set<int>::const_iterator it(boundaryVerts.begin());
        std::advance(it, distr1(generator));
        int newHoleVert = *it;
        vertInHole.insert(newHoleVert);
        for (auto f : mesh.vert(newHoleVert)->getAdjacentFaces()) {
            if (faceInHole.find(f->getFaceIndex()) != faceInHole.end())
                continue;       // face already considered     
            faceInHole.insert(f->getFaceIndex());
            for (int k = 0; k < 3; ++k) {
                int faceVertIdx = f->getVertexIndex(k);
                if (vertInHole.find(faceVertIdx) == vertInHole.end())
                    boundaryVerts.insert(faceVertIdx);
            }
        }
        boundaryVerts.erase(newHoleVert);

        for (auto iter = boundaryVerts.begin(); iter != boundaryVerts.end();) {
            int vIdx = *iter;
            bool encompassed = true;
            for (auto f : mesh.vert(vIdx)->getAdjacentFaces()) {
                if (faceInHole.find(f->getFaceIndex()) == faceInHole.end()) {
                    // found a face not in hole, so vIdx is not in hole yet
                    encompassed = false;
                    break;
                }
            }
            if (encompassed) {
                iter = boundaryVerts.erase(iter);
                vertInHole.insert(vIdx);
            }
            else iter++;
        }
    }

    result.mHoleFaces = vector < int > {faceInHole.begin(), faceInHole.end()};
    result.mHoleVerts = vector < int > {vertInHole.begin(), vertInHole.end()};
    result.mHoleBoundaryVerts = vector < int > {boundaryVerts.begin(), boundaryVerts.end()};
    return result;
}
