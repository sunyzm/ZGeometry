#pragma once
#include <string>
#include <vector>
#include <map>
#include <ZGeom/ZGeom.h>
#include "MeshLaplacian.h"
#include "global.h"

class MeshHelper {
public:
    struct ExtraInfo {
        int ref_vert_idx;
        ZGeom::Vec3d ref_vert_pos;
        std::map<int, ZGeom::Vec3d> handles;
        int active_handle_idx;
        std::map<LaplacianType, MeshLaplacian> laplacians;
        std::map<LaplacianType, ZGeom::EigenSystem> eigensys;
    };

public:
	MeshHelper();
    MeshHelper(MeshHelper && mh);

	void init(CMesh* tm);
    CMesh* getMeshByName(const std::string mesh_descript);    
    void addMesh(std::unique_ptr<CMesh> && newMesh, const std::string description = "");
    void nextMesh();
    void prevMesh();
    void revertOriginal();
    void removeCurrentMesh();
    CMesh* getOriginalMesh() const;
    CMesh* getMesh() const;
    ExtraInfo& getMeshExtraInfo();
    const ExtraInfo& getMeshExtraInfo() const;
    
    /* Laplacian and MHB related */
    bool hasLaplacian(LaplacianType lap_type);
    const MeshLaplacian& getMeshLaplacian(LaplacianType lap_type);
    void constructLaplacian(LaplacianType lap_type = CotFormula);

    bool hasEigenSystem(LaplacianType lap_type);
    ZGeom::EigenSystem& getEigenSystem(LaplacianType lap_type);
    void computeLaplacian(int nEigFunc, LaplacianType laplacianType = CotFormula);
    ZGeom::EigenSystem& prepareEigenSystem(LaplacianType lap_type = CotFormula, int num_eig = -1);
    
    /* editing */
    void addNewHandle(int hIdx);
    int getActiveHandle() const;
    void setActiveHandle(int h);
    std::map<int, ZGeom::Vec3d>& getHandles();
    const std::map<int, ZGeom::Vec3d>& getHandles() const;
    void clearAllHandles();
    int  getRefPointIndex() const;
    void setRefPointIndex(int i);
    const ZGeom::Vec3d& getRefPointPosition() const;
    void setRefPointPosition(int x, int y, int z);

private:
    MeshHelper(const MeshHelper&);
    MeshHelper& operator = (const MeshHelper&);

    bool isMHBCacheValid(const std::string& pathMHB, int eig_num);
    std::string generateMHBPath(const std::string& prefix, LaplacianType lap_type);

public:
    std::vector<std::unique_ptr<CMesh>> mesh_history_;
    std::vector<std::unique_ptr<ExtraInfo>> mesh_extra_infos_;
    int cur_mesh_idx_;

private:
    std::map<LaplacianType, MeshLaplacian>& getAllLaplacians() {
        return getMeshExtraInfo().laplacians;
    }
    std::map<LaplacianType, ZGeom::EigenSystem>& getAllEigensys() {
        return getMeshExtraInfo().eigensys;
    }
};
