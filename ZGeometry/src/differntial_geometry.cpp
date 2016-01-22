#include "differential_geometry.h"
#include <ZGeom/spectral_analysis.h>
#include <ZGeom/zassert.h>
#define USE_AMP

using namespace std;

void computeSGWMat(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW)
{
    ZGeom::runtime_assert(waveletScaleNum >= 1, "Illegal waveletScaleNum");
    const int vertCount = mhb.eigVecSize();
    const int eigCount = mhb.eigVecCount();
    const int nWaveletScales = waveletScaleNum;
    const double *pEigVals = mhb.getAllEigValsPtr();
    ZGeom::DenseMatrixd matEigVecs = mhb.toDenseMatrix();
    double *pEigVec = matEigVecs.raw_ptr();

    std::function<double(double)> genG = [](double x) {
        if (x < 1) return x*x;
        else if (x <= 2) return (-5. + 11.*x - 6.*x*x + x*x*x);
        else return 4.0 / x / x;
    };

    if (waveletScaleNum >= 2)
    {
        const int nScalingScales = 1;
        const int totalAtomCount = vertCount*(nWaveletScales + 1);
        matSGW.resize(totalAtomCount, vertCount);
        double maxEigVal = mhb.getEigVal(mhb.eigVecCount() - 1);
        const double K = 20.0;
        double minEigVal = maxEigVal / K;
        double minT = 2. / maxEigVal, maxT = 2. / minEigVal;
        const double tMultiplier = std::pow(maxT / minT, 1.0 / double(nWaveletScales - 1));
        std::vector<double> vWaveletScales(nWaveletScales);
        for (int s = 0; s < nWaveletScales; ++s)
            vWaveletScales[s] = minT * std::pow(tMultiplier, s);

        //////////////////////////////////////////////////////////////////////////
        // compute SGW matrix with AMP
        std::vector<double> vDiag(eigCount);
        //// compute wavelet functions
        for (int s = 0; s < nWaveletScales; ++s)
        {
            for (int i = 0; i < eigCount; ++i) {
                vDiag[i] = genG(vWaveletScales[s] * pEigVals[i]);
            }
            double *pResult = matSGW.raw_ptr() + vertCount * vertCount * s;
#ifdef USE_AMP
            ZGeom::quadricFormAMP(
#else		
            ZGeom::quadricForm(
#endif
                vertCount, eigCount, pEigVec, &vDiag[0], pResult);
        }
        //// compute scaling functions 	
        const double gamma = 1.3849;
        std::function<double(double)> genH = [=](double x) {
            return gamma * std::exp(-std::pow(x / (0.6*minEigVal), 4));
        };
        for (int i = 0; i < eigCount; ++i) {
            vDiag[i] = genH(pEigVals[i]);
        }
        double *pResult = matSGW.raw_ptr() + vertCount * vertCount * nWaveletScales;
#ifdef USE_AMP
        ZGeom::quadricFormAMP(
#else		
        ZGeom::quadricForm(
#endif
            vertCount, eigCount, pEigVec, &vDiag[0], pResult);
    }
    /// single-scale wavelets
    else if (waveletScaleNum == 1)
    {
        matSGW.resize(vertCount, vertCount);
        std::vector<double> vDiag(eigCount);
        double maxEigVal = mhb.getEigVal(mhb.eigVecCount() - 1);
        double minT = 2. / maxEigVal;
        for (int i = 0; i < eigCount; ++i) {
            vDiag[i] = genG(minT * pEigVals[i]);
        }
#ifdef USE_AMP
        ZGeom::quadricFormAMP(
#else		
        ZGeom::quadricForm(
#endif
            vertCount, eigCount, pEigVec, &vDiag[0], matSGW.raw_ptr());
    }
}

void computeSGWMat2(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW)
{
    const int vertCount = mhb.eigVecSize();
    const int eigCount = mhb.eigVecCount();

    std::function<double(double)> generator1 = [](double x) { return x * std::exp(-x); };
    std::function<double(double)> generator2 = [](double x) { return std::exp(-x - 1); };

    const int nScales = 5;
    const int totalAtomCount = vertCount*(waveletScaleNum + 1);
    double maxEigVal = mhb.getEigVal(mhb.eigVecCount() - 1);
    double minEigVal = maxEigVal / 20.;//mhb.getEigVal(1);
    double maxT = 1.0 / minEigVal;
    double minT = 1.0 / maxEigVal;
    const double tMultiplier = std::pow(maxT / minT, 1.0 / double(nScales - 1));

    std::vector<double> waveletScales(nScales);
    for (int s = 0; s < nScales; ++s) waveletScales[s] = minT * std::pow(tMultiplier, s);

    matSGW.resize(totalAtomCount, vertCount);

    ZGeom::DenseMatrixd matEigVecs(eigCount, vertCount);
    const double *pEigVals = &(mhb.getAllEigVals()[0]);
    double *pEigVec = matEigVecs.raw_ptr();
    for (int i = 0; i < eigCount; ++i)
        std::copy_n(mhb.getEigVec(i).c_ptr(), vertCount, pEigVec + i*vertCount);


    //////////////////////////////////////////////////////////////////////////
    // compute SGW with AMP
    //
    std::vector<double> vDiag(eigCount);
    for (int s = 0; s < nScales; ++s) {
        for (int i = 0; i < eigCount; ++i) {
            vDiag[i] = generator1(waveletScales[s] * pEigVals[i]);
        }
        double *pResult = matSGW.raw_ptr() + vertCount * vertCount * s;
        ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);
    }

    for (int i = 0; i < eigCount; ++i) {
        vDiag[i] = generator2(pEigVals[i]);
    }
    double *pResult = matSGW.raw_ptr() + vertCount * vertCount * nScales;
    ZGeom::quadricFormAMP(vertCount, eigCount, pEigVec, &vDiag[0], pResult);

}

void computeGeometricLaplacianCoordinate(const CMesh& mesh, const MeshCoordinates& eCoord, MeshCoordinates& lCoord)
{
    using ZGeom::Vec3d;
    assert(mesh.vertCount() == eCoord.size());
    int vertCount = mesh.vertCount();
    lCoord.resize(vertCount);

    for (int i = 0; i < vertCount; ++i)	{
        Vec3d vi = eCoord.getVertCoord(i);
        std::vector<int> neighbors = mesh.getVertNeighborVerts(i, 1, false);
        Vec3d vavg(0, 0, 0);
        double weightSum(0);
        for (int j : neighbors) {
            Vec3d vj = eCoord.getVertCoord(j);
            double weight = 1.0 / (vi - vj).length();
            vavg += vj * weight;
            weightSum += weight;
        }
        vavg /= weightSum;
        lCoord[i] = vi - vavg;
    }
}
