#include "MeshHelper.h"
#include <fstream>
#include <stdexcept>
#include <set>
#include <algorithm>
#include <functional>
#include <ppl.h>
#include <mkl.h>
#include <ZGeom/util.h>
#include <ZGeom/MatVecArithmetic.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include "differential_geometry.h"

using namespace std;
using ZGeom::MatlabEngineWrapper;
using ZGeom::EigenSystem;

MeshHelper::MeshHelper()
{
    cur_mesh_idx_ = -1;
}

 MeshHelper::MeshHelper(MeshHelper && mh)
 {
     mesh_history_ = std::move(mh.mesh_history_);
     mesh_extra_infos_ = std::move(mh.mesh_extra_infos_);
     cur_mesh_idx_ = mh.cur_mesh_idx_;
 }

void MeshHelper::init(CMesh* tm)
{
    assert(tm != nullptr);
    mesh_history_.clear();
    mesh_history_.resize(1);
    mesh_history_[0].reset(tm);

    mesh_extra_infos_.clear();
    mesh_extra_infos_.emplace_back(new ExtraInfo());
    mesh_extra_infos_.back()->ref_vert_idx = 0;
    mesh_extra_infos_.back()->ref_vert_pos = tm->vertPos(0);

    cur_mesh_idx_ = 0;
    tm->addAttr<std::string>("original_mesh", CMesh::StrAttrMeshDescription, AR_UNIFORM, AT_STRING);
}

CMesh* MeshHelper::getMesh() const
{
    return mesh_history_[cur_mesh_idx_].get();
}

CMesh* MeshHelper::getMeshByName(const std::string mesh_descript)
{
    for (int k = 0; k < mesh_history_.size(); ++k) {
        if (mesh_history_[k]->getMeshDescription() == mesh_descript)
            return mesh_history_[k].get();
    }

    return nullptr; 
}

CMesh* MeshHelper::getOriginalMesh() const
{
    return mesh_history_[0].get();
}

void MeshHelper::addMesh(std::unique_ptr<CMesh> && newMesh, const std::string description)
{
    mesh_history_.push_back(std::move(newMesh));
    CMesh &mesh = *mesh_history_.back();
    mesh.setMeshDescription(description);

    mesh_extra_infos_.emplace_back(new ExtraInfo());
    mesh_extra_infos_.back()->ref_vert_idx = 0;
    mesh_extra_infos_.back()->ref_vert_pos = mesh.vertPos(0);

    cur_mesh_idx_ = (int)mesh_history_.size() - 1;
}

void MeshHelper::removeCurrentMesh()
{
    if (cur_mesh_idx_ == 0) {
        cout << "Original mesh cannot be removed!" << std::endl;
        return;
    }

    mesh_history_.erase(mesh_history_.begin() + cur_mesh_idx_);
    mesh_extra_infos_.erase(mesh_extra_infos_.begin() + cur_mesh_idx_);

    cur_mesh_idx_ = cur_mesh_idx_ % mesh_history_.size();
    std::cout << "Switch to: " << mesh_history_[cur_mesh_idx_]->getMeshDescription() << std::endl;
}

void MeshHelper::nextMesh()
{
    if (mesh_history_.size() <= 1) return;
    cur_mesh_idx_ = (cur_mesh_idx_ + 1) % mesh_history_.size();    
}

void MeshHelper::prevMesh()
{
    if (mesh_history_.size() <= 1) return;
    cur_mesh_idx_ = (cur_mesh_idx_ + mesh_history_.size() - 1) % mesh_history_.size();
}

void MeshHelper::revertOriginal()
{
    cur_mesh_idx_ = 0;
}

MeshHelper::ExtraInfo& MeshHelper::getMeshExtraInfo()
{
    return *mesh_extra_infos_[cur_mesh_idx_];
}

const MeshHelper::ExtraInfo& MeshHelper::getMeshExtraInfo() const
{
    return *mesh_extra_infos_[cur_mesh_idx_];
}

bool MeshHelper::hasLaplacian(LaplacianType lap_type)
{
    auto& all_laplacians = getMeshExtraInfo().laplacians;
    return all_laplacians.find(lap_type) != all_laplacians.end();
}

const MeshLaplacian& MeshHelper::getMeshLaplacian(LaplacianType lap_type)
{
    if (!hasLaplacian(lap_type)) constructLaplacian(lap_type);
    return getMeshExtraInfo().laplacians[lap_type];
}

void MeshHelper::constructLaplacian( LaplacianType lap_type /*= CotFormula*/ )
{
	if (hasLaplacian(lap_type)) return;
    getAllLaplacians().insert(make_pair(lap_type, MeshLaplacian()));
    MeshLaplacian& laplacian = getAllLaplacians()[lap_type];
    (laplacian.*(laplacian.getConstructFunc(lap_type)))(getMesh());      
}

bool MeshHelper::hasEigenSystem(LaplacianType lap_type)
{
    auto& all_eigensys = getMeshExtraInfo().eigensys;
    return all_eigensys.find(lap_type) != all_eigensys.end();
}

ZGeom::EigenSystem& MeshHelper::getEigenSystem(LaplacianType lap_type)
{
    ZGeom::runtime_assert(hasEigenSystem(lap_type), "Requested MHB doesn't exist!");
    return getMeshExtraInfo().eigensys[lap_type];
}

void MeshHelper::computeLaplacian( int num_eig, LaplacianType lap_type /*= CotFormula*/ )
{
    ZGeom::logic_assert(hasLaplacian(lap_type), "laplacian hasn't been constructed");
    int vert_count = getMesh()->vertCount();
    ZGeom::logic_assert(num_eig < vert_count - 1, "Invalid num_eig!");
    if (num_eig <= 0) num_eig = vert_count - 2;

    std::string path_eigensys = generateMHBPath("cache/", lap_type);
    int using_cache = gSettings.LOAD_MHB_CACHE;
    getAllEigensys().insert(make_pair(lap_type, EigenSystem()));
    EigenSystem& new_eigensys = getEigenSystem(lap_type);
    if (using_cache > 0 && isMHBCacheValid(path_eigensys, num_eig)) {
        new_eigensys.load(path_eigensys);
    }
    else {
        std::cout << "==== Compute eigen-decomposition ====\n";
        getMeshLaplacian(lap_type).meshEigenDecompose(num_eig, &g_engineWrapper, new_eigensys);
        new_eigensys.save(path_eigensys);
    }
}

std::string MeshHelper::generateMHBPath( const std::string& prefix, LaplacianType laplacianType )
{
	std::string s_idx = "0";
	s_idx[0] += (int)laplacianType;
    std::string pathMHB = prefix + getMesh()->getMeshName() + ".mhb." + s_idx;
	return pathMHB;
}

bool MeshHelper::isMHBCacheValid( const std::string& path_mhb, int eigen_num )
{
	if (!fileExist(path_mhb)) return false;

	int nEig, nSize;
	ifstream ifs(path_mhb.c_str(), ios::binary);
	ifs.read((char*)&nEig, sizeof(int));
	ifs.read((char*)&nSize, sizeof(int));
	ifs.close();

    if (nEig != eigen_num || nSize != getMesh()->vertCount()) return false;

	return true;
}

ZGeom::EigenSystem& MeshHelper::prepareEigenSystem(LaplacianType lap_type /*= CotFormula*/, int num_eig /*= -1*/)
{
    int vert_count = getMesh()->vertCount();
    if (num_eig >= vert_count - 1) {
        throw std::runtime_error("Invalid num_eig requested!");
    }
    if (num_eig <= 0) num_eig = vert_count - 2;
    
    if (!hasLaplacian(lap_type)) {
        constructLaplacian(lap_type);
    }
    if (!hasEigenSystem(lap_type) || getEigenSystem(lap_type).eigVecCount() < num_eig) {
        computeLaplacian(num_eig, lap_type);
    }
         
    return getEigenSystem(lap_type);
}

void MeshHelper::addNewHandle( int hIdx )
{
    auto& handles = getMeshExtraInfo().handles;
    auto iter = handles.find(hIdx);
    if (iter != handles.end()) handles.erase(iter);
    else handles.insert(std::make_pair(hIdx, getMesh()->vert(hIdx)->pos()));
}

int MeshHelper::getActiveHandle() const
{
    return getMeshExtraInfo().active_handle_idx;
}

void MeshHelper::clearAllHandles()
{
    getMeshExtraInfo().handles.clear();
}

void MeshHelper::setActiveHandle(int h)
{
    getMeshExtraInfo().active_handle_idx = h;
}

std::map<int, ZGeom::Vec3d>& MeshHelper::getHandles()
{
    return getMeshExtraInfo().handles;
}

const std::map<int, ZGeom::Vec3d>& MeshHelper::getHandles() const
{
    return getMeshExtraInfo().handles;
}

int MeshHelper::getRefPointIndex() const
{
    return getMeshExtraInfo().ref_vert_idx;
}

void MeshHelper::setRefPointIndex(int i)
{
    getMeshExtraInfo().ref_vert_idx = i;
}

const ZGeom::Vec3d& MeshHelper::getRefPointPosition() const
{
    return getMeshExtraInfo().ref_vert_pos;
}

void MeshHelper::setRefPointPosition(int x, int y, int z)
{
    getMeshExtraInfo().ref_vert_pos = ZGeom::Vec3d(x, y, z);
}
