#include "arithmetic.h"
#include <stdexcept>
#include <amp.h>

namespace ZGeom {

double triArea(double e1, double e2, double e3)
{    if (!validAsTriangle(e1, e2, e3)) throw std::runtime_error("Invalid triangle side lengths");

    double p = (e1 + e2 + e3) / 2.0;
    return std::sqrt(p*(p - e1)*(p - e2)*(p - e3));
}

std::pair<Vec3d, double> ZGeom::circumcenter(Vec3d p1, Vec3d p2, Vec3d p3)
{
    // conf. http://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates
    
    double d = ((p1 - p2) ^ (p2 - p3)).length();
    double radius = ((p1 - p2).length()*(p2 - p3).length()*(p3 - p1).length()) / (2 * d);
    double alpha = (sqr((p2 - p3).length()) * (p1 - p2).dot(p1 - p3)) / (2 * d*d);
    double beta = (sqr((p1 - p3).length()) * (p2 - p1).dot(p2 - p3)) / (2 * d*d);
    double gamma = (sqr((p1 - p2).length()) * (p3 - p1).dot(p3 - p2)) / (2 * d*d);
    Vec3d center = p1 * alpha + p2 * beta + p3 * gamma;

    return std::make_pair(center, radius);
}

double cosTriSides(double e1, double e2, double e3)
{
    if (!validAsTriangle(e1, e2, e3)) throw std::runtime_error("Invalid triangle side lengths");
    return (e1*e1 + e2*e2 - e3*e3) / (2.0*e1*e2);
}

void triangleCot(double a, double b, double c, double &cotan_a, double &cotan_c)
{
    if (!validAsTriangle(a, b, c)) throw std::runtime_error("Invalid triangle side lengths");
    double cosa = (b*b + c*c - a*a) / (2.0*b*c);
    double cosc = (b*b + a*a - c*c) / (2.0*b*a);
    cotan_a = cosa / std::sqrt(1 - cosa*cosa);
    cotan_c = cosc / std::sqrt(1 - cosc*cosc);
}

double calMixedTriArea(double a, double b, double c)
{
    double cosa = (b*b + c*c - a*a) / (2.0*b*c);
    double cosc = (b*b + a*a - c*c) / (2.0*b*a);
    double cotan_a = cosa / sqrt(1.0 - cosa*cosa);
    double cotan_c = cosc / sqrt(1.0 - cosc*cosc);
    double s = (a + b + c) / 2.0;

    if (a*a + c*c < b*b) {        
        return sqrt(s*(s - a)*(s - b)*(s - c)) / 2.0;
    }
    else if (a*a + b*b < c*c || b*b + c*c < a*a) {
         return sqrt(s*(s - a)*(s - b)*(s - c)) / 4.0;
    }
    else {
        return (a*a*cotan_a + c*c*cotan_c) / 8.0;
    }
}

double calMixedTriArea(double a, double b, double c, double& cotan_a, double& cotan_c)
{
    double cosa = (b*b + c*c - a*a) / (2.0*b*c);
    double cosc = (b*b + a*a - c*c) / (2.0*b*a);
    cotan_a = cosa / sqrt(1.0 - cosa*cosa);
    cotan_c = cosc / sqrt(1.0 - cosc*cosc);
    double s = (a + b + c) / 2.0;

    if (a*a + c*c < b*b) {
        return sqrt(s*(s - a)*(s - b)*(s - c)) / 2.0;
    }
    else if (a*a + b*b < c*c || b*b + c*c < a*a) {
        return sqrt(s*(s - a)*(s - b)*(s - c)) / 4.0;
    }
    else {
        return (a*a*cotan_a + c*c*cotan_c) / 8.0;
    }
}

double calHalfMixedTriArea(double a, double b, double c, double& cotan_a)
{
    if (a*a + c*c < b*b) {
        double s = (a + b + c) / 2.0;
        return sqrt(s*(s - a)*(s - b)*(s - c)) / 4.0;
    }
    else if (a*a + b*b < c*c || b*b + c*c < a*a) {
        double s = (a + b + c) / 2.0;
        return sqrt(s*(s - a)*(s - b)*(s - c)) / 8.0;
    }
    else {
        double cosa = (b*b + c*c - a*a) / (2.0*b*c);
        cotan_a = cosa / sqrt(1 - cosa*cosa);
        return (a*a*cotan_a) / 8.0;
    }
}

void quadricForm(int dim1, int dim2, double* mat1, double* diag, double *matResult)
{
    using namespace Concurrency;
    parallel_for(0, dim1, [=](int i) {
        for (int j = 0; j <= i; ++j) {
            double sum(0);
            for (int k = 0; k < dim2; ++k) {
                sum += diag[k] * mat1[dim1*k + i] * mat1[dim1*k + j];
            }
            matResult[dim1*i + j] = matResult[dim1*j + i] = sum;
        }
    });
}


void quadricFormAMP(int dim1, int dim2, double* mat1, double* diag, double *matResult)
{
    using namespace Concurrency;
    // Y=X'*Q*X
    array_view<double, 2> X(dim2, dim1, mat1);
    array_view<double, 1> Q(dim2, diag);
    int maxDim = 1000;
    int maxBlock = maxDim * maxDim;

    if (dim1 < maxDim) {
        array_view<double, 2> Y(dim1, dim1, matResult);
        Y.discard_data();
        parallel_for_each(Y.extent, [=](index<2> idx) restrict(amp) {
            int row = idx[0], col = idx[1];
            if (row >= col) {
                Y[idx] = 0;
                for (int k = 0; k < dim2; ++k) {
                    Y[idx] += Q(k) * X(k, row) * X(k, col);
                }
                Y(col, row) = Y[idx];
            }
        });
        Y.synchronize();
    }
    else
    {
        int maxIter = dim1 * dim1 / maxBlock + 1;
        for (int l = 0; l < maxIter; ++l) {
            int curPtr = maxBlock * l;
            int blockSize = maxBlock;
            if (l == maxIter - 1) blockSize = dim1*dim1 - curPtr;
            array_view<double, 1> Y(blockSize, matResult + curPtr);
            Y.discard_data();
            parallel_for_each(Y.extent, [=](index<1> idx) restrict(amp) {
                int pos = idx[0] + curPtr;
                int row = pos / dim1, col = pos % dim1;
                if (row >= col) {
                    Y[idx] = 0;
                    for (int k = 0; k < dim2; ++k) {
                        Y[idx] += Q(k) * X(k, row) * X(k, col);
                    }
                }
            });
            Y.synchronize();
        }
        parallel_for(0, dim1, [&](int i) {
            //for (int i = 0; i < dim1; ++i) {
            for (int j = i + 1; j < dim1; ++j) {
                matResult[i*dim1 + j] = matResult[j*dim1 + i];
            }
        }
        );
    }
}

void matVecMulAMP(int dim1, int dim2, double *mat, double *vec, double *vResult)
{
    using namespace Concurrency;
    // Y = M*X
    array_view<double, 2> M(dim1, dim2, mat);
    array_view<double, 1> X(dim2, vec);
    array_view<double, 1> Y(dim1, vResult);
    Y.discard_data();

    parallel_for_each(Y.extent, [=](index<1> idx) restrict(amp) {
        int row = idx[0];
        Y[idx] = 0;
        for (int k = 0; k < dim2; ++k) {
            Y[idx] += M(row, k) * X(k);
        }
    });

    Y.synchronize();
}

std::vector<double> linspace(double x1, double x2, int N)
{
    assert(N >= 1);
    if (N == 1) return std::vector < double > {x2};
    std::vector<double> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = ((N - 1 - i) * x1 + i * x2) / double(N - 1);
    }
    return result;
}

}   // end of namespace