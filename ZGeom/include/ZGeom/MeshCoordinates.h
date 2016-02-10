#ifndef ZGEOM_MESH_COORDINATES_H
#define ZGEOM_MESH_COORDINATES_H
#include <cassert>
#include <vector>
#include "util.h"
#include "Vec3.h"
#include "VecN.h"
#include "DenseMatrix.h"

class MeshCoordinates
{
public:
    MeshCoordinates() {}
    MeshCoordinates(int meshSize) { resize(meshSize); }
    MeshCoordinates(int mesh_size, double *cx, double *cy, double *cz);
    MeshCoordinates(int mesh_size, const ZGeom::VecNd& v1, const ZGeom::VecNd& v2, const ZGeom::VecNd& v3);
    MeshCoordinates(int mesh_size, const std::vector<ZGeom::VecNd>& vCoords);
    MeshCoordinates(const MeshCoordinates& mc) = default;
    MeshCoordinates(MeshCoordinates&& mc);
    MeshCoordinates& operator = (const MeshCoordinates& mc) = default;
    MeshCoordinates& operator = (MeshCoordinates&& mc);

    bool empty() const { return coord_data.empty(); }
    int size() const { return (int)coord_data.size()/3; }
    void resize(int n) { coord_data.resize(n * 3, 0); }
    double* data() { return coord_data.data(); }
    const double* data() const { return coord_data.data(); }
    const ZGeom::VecNd getCoordVec(int c) const;
    const ZGeom::VecNd getXCoord() const { return getCoordVec(0); }
    const ZGeom::VecNd getYCoord() const { return getCoordVec(1); }
    const ZGeom::VecNd getZCoord() const { return getCoordVec(2); }
    double* getCoordData(int c) { return data() + size() * c; }
    double* xCoordData() { return data(); }
    double* yCoordData() { return data() + size(); }
    double* zCoordData() { return data() + size() * 2; }
    void setCoord(int c, const ZGeom::VecNd vec);

    void addWith(double *cx, double *cy, double *cz);
    void addWith(int c, double* raw);
    void addWith(int c, const ZGeom::VecNd& vec_coord);
    MeshCoordinates add(const MeshCoordinates& mc2) const;
    MeshCoordinates substract(const MeshCoordinates& mc2) const;

    ZGeom::Vec3d getVertCoord(int v) const;
    void setVertCoord(int vIdx, ZGeom::Vec3d vec);
    ZGeom::Vec3d operator [] (int v) const { return getVertCoord(v); }
    double& operator() (int idx, int c);
    double& elem(int idx, int c);

    double difference(const MeshCoordinates& mc2) const;
    
    ZGeom::VecNd vertDifference(const MeshCoordinates& mc2) const;
    std::vector<ZGeom::VecNd> to3Vec() const;

    // convert from/to N*3 matrix
    ZGeom::DenseMatrixd toDenseMatrix() const;
    void fromDenseMatrix(const ZGeom::DenseMatrixd& mat);

    void scaleToUnitbox();

private:
    std::vector<double> coord_data;
};

#endif