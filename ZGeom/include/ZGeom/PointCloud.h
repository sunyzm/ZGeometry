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
    PointCloud<T>& operator = (const PointCloud<T>& pc) = default;
    PointCloud<T>& operator = (PointCloud<T>&& pc) { vPoints = std::move(pc.vPoints); return *this; }
    PointCloud(const PointCloud<T>& pc) = default;        
    PointCloud(PointCloud<T>&& pc) { *this = std::move(pc); }
    PointCloud(const std::vector<T>& vp) : vPoints(vp) {}    
    PointCloud(std::vector<T>&& vp) { vPoints = std::move(vp); }
    
    void print(const std::string& filename);
    int size() const { return (int)vPoints.size(); }
    const std::vector<T>& getPoints() const { return vPoints; }
    std::vector<T>& getPoints() { return vPoints; }
    const T& getPoint(int idx) const { return vPoints[idx]; }
    T& getPoint(int idx) { return vPoints[idx]; }

    const T& operator [] (int idx) const { return vPoints[idx]; }
    T& operator[] (int idx) { return vPoints[idx]; }

private:
    std::vector<T> vPoints;
};

typedef PointCloud<Vec3d> PointCloud3d;
typedef PointCloud<VecNd> PointCloudNd;

}   // end of namespace


#endif