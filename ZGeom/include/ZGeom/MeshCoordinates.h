#ifndef ZGEOM_MESH_COORDINATES_H
#define ZGEOM_MESH_COORDINATES_H
#include <cassert>
#include <vector>
#include "Vec3.h"
#include "VecN.h"
#include "DenseMatrix.h"

class MeshCoordinates
{
public:
    MeshCoordinates() : mSize(0) {}
    MeshCoordinates(int meshSize) { resize(meshSize); }

    MeshCoordinates(int meshSize, double *cx, double *cy, double *cz)
    {
        resize(meshSize);
        std::copy_n(cx, meshSize, mCoordX.c_ptr());
        std::copy_n(cy, meshSize, mCoordY.c_ptr());
        std::copy_n(cz, meshSize, mCoordZ.c_ptr());
    }

    MeshCoordinates(int meshSize, const ZGeom::VecNd& v1, const ZGeom::VecNd& v2, const ZGeom::VecNd& v3)
    {
        assert(meshSize == v1.size() && v1.size() == v2.size() && v2.size() == v3.size());
        mSize = meshSize;
        mCoordX = v1;
        mCoordY = v2;
        mCoordZ = v3;
    }

    MeshCoordinates(int meshSize, const std::vector<ZGeom::VecNd>& vCoords)
    {
        assert(vCoords.size() >= 3 && meshSize == vCoords[0].size());
        mSize = meshSize;
        mCoordX = vCoords[0];
        mCoordY = vCoords[1];
        mCoordZ = vCoords[2];
    }

    MeshCoordinates(const MeshCoordinates& mc)
    {
        this->mSize = mc.mSize;
        this->mCoordX = mc.mCoordX;
        this->mCoordY = mc.mCoordY;
        this->mCoordZ = mc.mCoordZ;
    }

    MeshCoordinates(MeshCoordinates&& mc)
    {
        this->mSize = mc.mSize;
        mc.mSize = 0;
        this->mCoordX = std::move(mc.mCoordX);
        this->mCoordY = std::move(mc.mCoordY);
        this->mCoordZ = std::move(mc.mCoordZ);
    }

    MeshCoordinates& operator = (const MeshCoordinates& mc)
    {
        this->mSize = mc.mSize;
        this->mCoordX = mc.mCoordX;
        this->mCoordY = mc.mCoordY;
        this->mCoordZ = mc.mCoordZ;
        return *this;
    }

    MeshCoordinates& operator = (MeshCoordinates&& mc)
    {
        this->mSize = mc.mSize;
        mc.mSize = 0;
        this->mCoordX = std::move(mc.mCoordX);
        this->mCoordY = std::move(mc.mCoordY);
        this->mCoordZ = std::move(mc.mCoordZ);
        return *this;
    }

    bool empty() const { return mSize == 0; }
    int size() const { return mSize; }

    void resize(int n)
    {
        mSize = n;
        mCoordX.resize(mSize, 0);
        mCoordY.resize(mSize, 0);
        mCoordZ.resize(mSize, 0);
    }

    void add(double *cx, double *cy, double *cz)
    {
        mCoordX.add(cx);
        mCoordY.add(cy);
        mCoordZ.add(cz);
    }

    MeshCoordinates add(const MeshCoordinates& mc2) const
    {
        assert(this->mSize == mc2.mSize);
        MeshCoordinates mc3(this->mSize);
        mc3.mCoordX = this->mCoordX + mc2.mCoordX;
        mc3.mCoordY = this->mCoordY + mc2.mCoordY;
        mc3.mCoordZ = this->mCoordZ + mc2.mCoordZ;
        return mc3;
    }

    MeshCoordinates substract(const MeshCoordinates& mc2) const
    {
        assert(this->mSize == mc2.mSize);
        MeshCoordinates mc3(this->mSize);
        mc3.mCoordX = this->mCoordX - mc2.mCoordX;
        mc3.mCoordY = this->mCoordY - mc2.mCoordY;
        mc3.mCoordZ = this->mCoordZ - mc2.mCoordZ;
        return mc3;
    }

    const ZGeom::VecNd& getCoordFunc(int i) const
    {
        switch (i)
        {
        case 0: return mCoordX;
        case 1: return mCoordY;
        case 2: return mCoordZ;
        default: throw std::logic_error("Invalid mesh coordinate");
        }
    }

    ZGeom::VecNd& getCoordFunc(int i)
    {
        switch (i)
        {
        case 0: return mCoordX;
        case 1: return mCoordY;
        case 2: return mCoordZ;
        default: throw std::logic_error("Invalid mesh coordinate");
        }
    }

    ZGeom::VecNd& getXCoord() { return mCoordX; }
    ZGeom::VecNd& getYCoord() { return mCoordY; }
    ZGeom::VecNd& getZCoord() { return mCoordZ; }
    const ZGeom::VecNd& getXCoord() const { return mCoordX; }
    const ZGeom::VecNd& getYCoord() const { return mCoordY; }
    const ZGeom::VecNd& getZCoord() const { return mCoordZ; }

    ZGeom::Vec3d getVertCoordinate(int v) const
    {
        assert(v >= 0 && v < mSize);
        return ZGeom::Vec3d(mCoordX[v], mCoordY[v], mCoordZ[v]);
    }

    ZGeom::Vec3d operator [] (int v) const { return getVertCoordinate(v); }
    double& operator() (int idx, int c) {
        if (c == 0) return mCoordX[idx];
        else if (c == 1) return mCoordY[idx];
        else if (c == 2) return mCoordZ[idx];
        else throw std::logic_error("Invalid coordinate inquiry");
    }

    double difference(const MeshCoordinates& mc2) const
    {
        assert(this->mSize == mc2.mSize);

        double errorSum(0);
        for (int i = 0; i < mSize; ++i) {
            errorSum += std::pow((this->getVertCoordinate(i) - mc2.getVertCoordinate(i)).length(), 2);
        }
        errorSum = std::sqrt(errorSum);
        return errorSum;
    }

    ZGeom::VecNd vertDifference(const MeshCoordinates& mc2) const
    {
        assert(mSize == mc2.size());
        ZGeom::VecNd vDiff(mSize);
        for (int i = 0; i < mSize; ++i)
            vDiff[i] = (getVertCoordinate(i) - mc2.getVertCoordinate(i)).length();
        return vDiff;
    }

    std::vector<ZGeom::VecNd> to3Vec() const { return std::vector<ZGeom::VecNd> {mCoordX, mCoordY, mCoordZ}; }
    
    ZGeom::DenseMatrixd toDenseMatrix() const 
    {
        ZGeom::VecNd coords[3] = { mCoordX, mCoordY, mCoordZ };
        ZGeom::DenseMatrixd mat(mSize, 3);
        for (int i = 0; i < mSize; ++i) {
            mat(i, 0) = mCoordX[i];
            mat(i, 1) = mCoordY[i];
            mat(i, 2) = mCoordZ[i];
        }
        return mat;
    }

    void fromDenseMatrix(const ZGeom::DenseMatrixd& mat)
    {
        for (int i = 0; i < mSize; ++i) {
            mCoordX[i] = mat(i, 0);
            mCoordY[i] = mat(i, 1);
            mCoordZ[i] = mat(i, 2);
        }
    }

private:
    int mSize;
    ZGeom::VecNd mCoordX, mCoordY, mCoordZ;
};

#endif