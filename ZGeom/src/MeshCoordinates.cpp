#include "MeshCoordinates.h"

void MeshCoordinates::addWith(double *cx, double *cy, double *cz)
{
    for (int i = 0; i < size(); ++i) {
        this->coord_data[i] += cx[i];
        this->coord_data[size() + i] = cy[i];
        this->coord_data[size() * 2 + i] = cz[i];
    }
}

void MeshCoordinates::addWith(int c, double* raw)
{
    assert(c >= 0 && c < 3);
    for (int i = 0; i < size(); ++i) {
        this->coord_data[i + size() * c] = raw[i];
    }
}

MeshCoordinates MeshCoordinates::add(const MeshCoordinates& mc2) const
{
    assert(size() == mc2.size());
    MeshCoordinates mc3(size());
    for (size_t i = 0; i < coord_data.size(); ++i) {
        mc3.coord_data[i] = this->coord_data[i] + mc2.coord_data[i];
    }
    return mc3;
}

MeshCoordinates MeshCoordinates::substract(const MeshCoordinates& mc2) const
{
    assert(size() == mc2.size());
    MeshCoordinates mc3(size());
    for (size_t i = 0; i < coord_data.size(); ++i) {
        mc3.coord_data[i] = this->coord_data[i] - mc2.coord_data[i];
    }
    return mc3;
}

MeshCoordinates::MeshCoordinates(int mesh_size, double *cx, double *cy, double *cz)
{
    resize(mesh_size);
    std::copy_n(cx, mesh_size, data());
    std::copy_n(cy, mesh_size, data() + mesh_size);
    std::copy_n(cz, mesh_size, data() + mesh_size * 2);
}

MeshCoordinates::MeshCoordinates(int mesh_size, const ZGeom::VecNd& v1, const ZGeom::VecNd& v2, const ZGeom::VecNd& v3)
{
    assert(mesh_size == v1.size() && v1.size() == v2.size() && v2.size() == v3.size());
    resize(mesh_size);
    std::copy_n(v1.c_ptr(), mesh_size, data());
    std::copy_n(v2.c_ptr(), mesh_size, data() + mesh_size);
    std::copy_n(v3.c_ptr(), mesh_size, data() + mesh_size * 2);
}

MeshCoordinates::MeshCoordinates(int mesh_size, const std::vector<ZGeom::VecNd>& vCoords)
{
    assert(vCoords.size() == 3 && mesh_size == vCoords[0].size());
    resize(mesh_size);
    std::copy_n(vCoords[0].c_ptr(), mesh_size, data());
    std::copy_n(vCoords[1].c_ptr(), mesh_size, data() + mesh_size);
    std::copy_n(vCoords[2].c_ptr(), mesh_size, data() + mesh_size * 2);
}

MeshCoordinates::MeshCoordinates(MeshCoordinates&& mc)
{
    this->coord_data = std::move(mc.coord_data);
}

MeshCoordinates& MeshCoordinates::operator=(MeshCoordinates&& mc)
{
    this->coord_data = std::move(mc.coord_data);
    return *this;
}

const ZGeom::VecNd MeshCoordinates::getCoordVec(int c) const
{
    assert(c >= 0 && c < 3);
    return ZGeom::VecNd(data() + size() * c, size());
}

ZGeom::Vec3d MeshCoordinates::getVertCoord(int vi) const
{
    return ZGeom::Vec3d(coord_data[vi], coord_data[vi + size()], coord_data[vi + size() * 2]);
}

void MeshCoordinates::setVertCoord(int vi, ZGeom::Vec3d vec)
{
    coord_data[vi] = vec[0];
    coord_data[vi + size()] = vec[1];
    coord_data[vi + size() * 2] = vec[2];
}

double& MeshCoordinates::operator()(int idx, int c)
{
    assert(c >= 0 && c <= 2);
    return coord_data[size() * c + idx];    
}

std::vector<ZGeom::VecNd> MeshCoordinates::to3Vec() const
{
    return std::vector<ZGeom::VecNd> {getCoordVec(0), getCoordVec(1), getCoordVec(2)};
}

ZGeom::DenseMatrixd MeshCoordinates::toDenseMatrix() const
{
    ZGeom::DenseMatrixd mat(size(), 3);
    for (int i = 0; i < size(); ++i) {
        for (int c = 0; c < 3; ++c)
            mat(i, c) = coord_data[size() * c + i];
    }
    return mat;
}

void MeshCoordinates::fromDenseMatrix(const ZGeom::DenseMatrixd& mat)
{
    ZGeom::logic_assert(mat.colCount() == 3);
    resize(mat.rowCount());
    for (int i = 0; i < size(); ++i) {
        for (int c = 0; c < 3; ++c)
            coord_data[size() * c + i] = mat(i, c);
    }
}

double MeshCoordinates::difference(const MeshCoordinates& mc2) const
{
    assert(this->size() == mc2.size());

    double errorSum(0);
    for (int i = 0; i < size(); ++i) {
        errorSum += std::pow((this->getVertCoord(i) - mc2.getVertCoord(i)).length(), 2);
    }
    errorSum = std::sqrt(errorSum);
    return errorSum;
}

ZGeom::VecNd MeshCoordinates::vertDifference(const MeshCoordinates& mc2) const
{
    assert(this->size() == mc2.size());
    ZGeom::VecNd vDiff(size());
    for (int i = 0; i < size(); ++i) {
        vDiff[i] = (getVertCoord(i) - mc2.getVertCoord(i)).length();
    }
    return vDiff;
}

void MeshCoordinates::setCoord(int c, const ZGeom::VecNd vec)
{
    assert(vec.size() == this->size() && c >= 0 && c < 2);
    std::copy_n(vec.c_ptr(), size(), data() + c * size());
}
