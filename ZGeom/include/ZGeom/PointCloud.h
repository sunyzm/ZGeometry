#ifndef ZGEOM_POINT_CLOUD_H
#define ZGEOM_POINT_CLOUD_H
#include <vector>
#include <string>
#include "Vec3.h"
#include "VecN.h"

namespace ZGeom {

template<typename T>
class PointCloud
{
public:
    PointCloud() = default;
    PointCloud(const PointCloud<T>& pc) = default;
    PointCloud(PointCloud<T>&& pc);
    PointCloud(const std::vector<T>& vp) : vPoints(vp) {}

    void print(const std::string& filename);
    int size() const { return (int)vPoints.size(); }

private:
    std::vector<T> vPoints;
};

typedef PointCloud<Vec3d> PointCloud3d;
typedef PointCloud<VecNd> PointCloudNd;

template<typename T>
PointCloud<T>::PointCloud(PointCloud<T>&& pc)
{
    vPoints = std::move(pc.vPoints);
}

}   // end of namespace


#endif