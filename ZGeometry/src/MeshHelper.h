#pragma once
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <ZGeom/ZGeom.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include "MeshLaplacian.h"
#include "global.h"

class MeshHelper
{
public:
	MeshHelper();
    MeshHelper(MeshHelper && mh);

	void init(CMesh* tm);
    void addMesh(std::unique_ptr<CMesh> && newMesh, const std::string description = "");
    void nextMesh();
    void prevMesh();
    void revertOriginal();
    CMesh* getMesh() const { return mMeshHistory[currentMeshIdx].get(); }
    CMesh* getMeshByName(const std::string mesh_descript);
    CMesh* getOriginalMesh() const { return mMeshHistory[0].get(); }
    
    void clearMeshRelated();

    /* editing */
    void addNewHandle(int hIdx);
    int getActiveHandle() const { return mActiveHandle; }
    void setActiveHandle(int h) { mActiveHandle = h; }
    std::map<int, ZGeom::Vec3d>& getHandles() { return mHandles; }
    const std::map<int, ZGeom::Vec3d>& getHandles() const { return mHandles; }
    void clearAllHandles() { mHandles.clear(); }
    int  getRefPointIndex() const { return mRefVert; }
    void setRefPointIndex(int i) { mRefVert = i; }
    void setRefPointPosition(int x, int y, int z) { mRefPos = ZGeom::Vec3d(x, y, z); }

	/* Laplacian and MHB related */
	void constructLaplacian(LaplacianType laplacianType = CotFormula);
	void decomposeLaplacian(int nEigFunc, LaplacianType laplacianType = CotFormula);
	const MeshLaplacian& getMeshLaplacian(LaplacianType laplacianType) const { return mMeshLaplacians[laplacianType]; }
	bool hasLaplacian(LaplacianType laplacianType) { return mMeshLaplacians[laplacianType].isLaplacianConstructed(); }
	bool isLaplacianDecomposed(LaplacianType laplacianType) { return !mMHBs[laplacianType].empty(); }
	void loadMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	void saveMHB(const std::string& path, LaplacianType laplacianType = CotFormula);
	const ZGeom::EigenSystem& getMHB(LaplacianType laplacianType) const { return mMHBs[laplacianType]; }
	bool isMHBCacheValid(const std::string& pathMHB, int eigenCount);
	std::string generateMHBPath(const std::string& prefix, LaplacianType laplacianType);
	const ZGeom::EigenSystem& prepareEigenSystem(const MeshLaplacian& laplaceMat, int eigenCount);
	double	calHK(int v1, int v2, double timescale) const;

private:
    MeshHelper(const MeshHelper&);
    MeshHelper& operator = (const MeshHelper&);

public:
    std::vector<std::unique_ptr<CMesh>> mMeshHistory;
    std::map<std::string, CMesh*> mMeshDescriptMap;
    int currentMeshIdx;

	int mRefVert;
	ZGeom::Vec3d mRefPos;
	std::map<int, ZGeom::Vec3d> mHandles;
	int mActiveHandle;

	MeshLaplacian mMeshLaplacians[LaplacianTypeCount];
	ZGeom::EigenSystem mMHBs[LaplacianTypeCount];
};

