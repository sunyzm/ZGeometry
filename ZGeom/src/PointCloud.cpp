#include "PointCloud.h"
#include <fstream>

namespace ZGeom {

template<> 
void PointCloud<Vec3d>::print(const std::string& filename)
{
    std::ofstream ofs(filename.c_str());
    for (auto l : vPoints) {
        ofs << l.x << ' ' << l.y << ' ' << l.z << std::endl;
    }
    ofs.close();
}

template<>
void PointCloud<VecNd>::print(const std::string& filename)
{
    std::ofstream ofs(filename.c_str());
    for (auto pt : vPoints) {
        for (int i = 0; i < (int)pt.size() - 1; ++i)
            ofs << pt[i] << ' ';
        ofs << pt[pt.size() - 1] << std::endl;
    }
}

}   // end of namespace

