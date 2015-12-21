#include "spectral_geometry.h"
#include <cmath>
#include <functional>
#include <vector>
#include <ppl.h>
#include "arithmetic.h"

using namespace std;
using namespace concurrency;

std::function<double(double, double)> heat_gen_func = [](double lambda, double tau){
    return std::exp(-lambda*tau);
};

std::function<double(double)> wavelet_gen_1 = [](double x) {
    if (x < 1.0) return x*x;
    else if (x <= 2.0) return (-5 + 11*x - 6*x*x + x*x*x);
    else return 4.0 / x / x;
};

std::vector<double> bivariateKernelRow(const ZGeom::EigenSystem& es, const std::vector<double>& diag_coeff, int vi) {
    int vert_num = es.eigVecSize();
    int eigen_num = es.eigVecCount();

    vector<double> result(vert_num);
    concurrency::parallel_for(0, vert_num, [&](int vj) {
        double sum = 0;
        for (int k = 0; k < eigen_num; ++k) {
            const auto& phi_k = es.getEigVec(k);
            sum += diag_coeff[k] * phi_k[vi] * phi_k[vj];
        }
        result[vj] = sum;
    });

    return result;
}

std::vector<double> bivariateKernelDiagonal(const ZGeom::EigenSystem& es, const std::vector<double>& diag_coeff) {
    int vert_num = es.eigVecSize();
    int eigen_num = es.eigVecCount();

    vector<double> result(vert_num);
    concurrency::parallel_for(0, vert_num, [&](int vi) {
        double sum = 0;
        for (int k = 0; k < eigen_num; ++k) {
            const auto& phi = es.getEigVec(k);
            sum += diag_coeff[k] * phi[vi] * phi[vi];
        }
        result[vi] = sum;
    });

    return result;
}

namespace ZGeom 
{

std::vector<double> calSpectralKernelSignature(const EigenSystem& es, double t, std::function<double(double, double)> gen_func)
{
    const int vert_count = es.eigVecSize();
    const int eigen_num = es.eigVecCount();

    std::vector<double> diag(eigen_num);
    for (int k = 0; k < eigen_num; ++k) {
        diag[k] = gen_func(es.getEigVal(k), t);
    }

    vector<double> result = bivariateKernelDiagonal(es, diag);
    return result;
}

DenseMatrixd calSpectralKernelMatrix(const EigenSystem& es, double t, std::function<double(double, double)> gen_func)
{
    const int vert_num = es.eigVecSize();
    const int eigen_num = es.eigVecCount();

    std::vector<double> diag(eigen_num);
    for (int k = 0; k < eigen_num; ++k) {
        diag[k] = gen_func(es.getEigVal(k), t);
    }
    DenseMatrixd mat_eigen_vec = es.toDenseMatrix();
    DenseMatrixd mat_hk(vert_num, vert_num);

    quadricFormAMP(vert_num, eigen_num, mat_eigen_vec.raw_ptr(), diag.data(), mat_hk.raw_ptr());

    return mat_hk;
}

double calBiharmonicDist(const EigenSystem& es, int v1, int v2)
{
    // BhDist(i,j) = sqrt(sum_{k=1}^{K}(phi_k(i)-phi_k(j))/lambda_k)^2 )
    double sum = 0;
    for (int k = 1; k < es.eigVecCount(); ++k) {
        sum += pow((es.getEigVec(k)[v1] - es.getEigVec(k)[v2]) / es.getEigVal(k), 2);
    }
    return std::sqrt(sum);
}

std::vector<double> calAllBiharmonicDist(EigenSystem& es, int vo, int eig_num)
{
    // BhDist(i,j) = sqrt(sum_{k=1}^{K}(phi_k(i)-phi_k(j))/lambda_k)^2 )
    int vert_count = es.eigVecSize(), eig_count = es.eigVecCount();
    if (eig_num > 0 && eig_num < eig_count) eig_count = eig_num;
    const DenseMatrixd& mat_eig_vec = es.getEigenMat();
    const vector<double>& vec_eig_vals = es.getAllEigVals();

    vector<double> result(vert_count, 0);
    parallel_for(0, vert_count, [&](int vi) {
        double sum(0);
        for (int k = 1; k < eig_count; ++k) {
            sum += ZGeom::sqr((mat_eig_vec(k, vo) - mat_eig_vec(k, vi)) / vec_eig_vals[k]);
        }
        result[vi] = std::sqrt(sum);
    });

    return result;
}

double calHeatKernel(const EigenSystem& es, int i, int j, double t)
{
    double sum = 0;
    for (int k = 0; k < es.eigVecCount(); ++k) {
        const auto& phi = es.getEigVec(k);
        sum += heat_gen_func(es.getEigVal(k), t) * phi[i] * phi[j];
    }
    return sum;
}

std::vector<double> calAllHeatKernel(const EigenSystem& es, int vi, double t)
{
    int vert_num = es.eigVecCount();
    int eigen_num = es.eigVecCount();
    vector<double> coeff_diag(eigen_num);

    for (int k = 0; k < eigen_num; ++k) {
        coeff_diag[k] = heat_gen_func(es.getEigVal(k), t);
    }

    return bivariateKernelRow(es, coeff_diag, vi);
}

std::vector<double> calHeatKernelSignature(const EigenSystem& hb, double t)
{
    return calSpectralKernelSignature(hb, t, heat_gen_func);
}

DenseMatrixd calHeatKernelMatrix(const EigenSystem& hb, double t)
{
    return calSpectralKernelMatrix(hb, t, heat_gen_func);
}

double calHeatTrace(const EigenSystem& es, double t)
{
    int eigen_num = es.eigVecCount();
    double sum(0);
    for (int k = 0; k < eigen_num; ++k) {
        sum += std::exp(-es.getEigVal(k) * t);
    }
    return sum;
}

std::vector<double> calAllSgwWavelet(const EigenSystem& es, double t, int vi)
{
    int vert_num = es.eigVecSize();
    int eigen_num = es.eigVecCount();

    vector<double> diag(eigen_num);
    for (int k = 0; k < eigen_num; ++k) {
        diag[k] = wavelet_gen_1(t * es.getEigVal(k));
    }

    vector<double> result = bivariateKernelRow(es, diag, vi);
    return result;
}

std::vector<double> computeSgwScales(const EigenSystem& es, int num_scales)
{
    if (num_scales <= 1) throw runtime_error("Number of scales must greater than 1!");

    int vert_num = es.eigVecSize();
    int eigen_num = es.eigVecCount();
    const double K = 20;
    const int J = num_scales;
    double max_eigen_val = es.getAllEigVals().back();
    double min_eigen_val = max_eigen_val / K;
    double min_t = 2. / max_eigen_val, max_t = 2. / min_eigen_val;
    double t_multiplier = pow(max_t/min_t, 1.0/(J - 1));
    
    vector<double> result(J);
    for (int s = 0; s < J; ++s) {
        result[s] = min_t * pow(t_multiplier, s);
    }
    return result;
}

} // end of namespace


