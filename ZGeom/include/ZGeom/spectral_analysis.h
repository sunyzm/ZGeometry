#ifndef SPECTRAL_WAVELETS_H
#define SPECTRAL_WAVELETS_H
#include <vector>
#include <functional>
#include "EigenSystem.h"
#include "DenseMatrix.h"

namespace ZGeom
{

std::vector<double>     calSpectralCoeff(const EigenSystem& es, const std::vector<double>& vert_func);

double                  calBiharmonicDist(const EigenSystem& es, int v1, int v2);
std::vector<double>     calAllBiharmonicDist(EigenSystem& es, int vi, int eig_num = -1);

std::vector<double>     calSpectralKernelSignature(const EigenSystem& es, double t, std::function<double(double, double)> gen_func);
DenseMatrixd            calSpectralKernelMatrix(const EigenSystem& hb, double t, std::function<double(double, double)> gen_func);

double                  calHeatKernel(const EigenSystem& es, int i, int j, double t);
std::vector<double>     calAllHeatKernel(const EigenSystem& es, int vi, double t);
DenseMatrixd            calHeatKernelMatrix(const EigenSystem& es, double t);
double                  calHeatTrace(const EigenSystem& es, double t);
std::vector<double>     calHeatKernelSignature(const EigenSystem& es, double t);

std::vector<double>     computeSgwScales(const EigenSystem& es, int num_scales);
std::vector<double>     calAllSgwWavelet(const EigenSystem& es, double t, int vi);
std::vector<double>     calSgwCoeff(const EigenSystem& es, double t, const std::vector<double>& mesh_func);

}
#endif