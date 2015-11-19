#include <vector>

namespace ZGeom {

struct VertCurvature {
    enum CurvatureType {MEAN, GAUSS, PRINCIPAL_1, PRINCIPAL_2};

    double mean_curv{ 0 };
    double gauss_curv{ 0 };
    double principal_curv_1{ 0 };
    double principal_curv_2{ 0 };
};

}

