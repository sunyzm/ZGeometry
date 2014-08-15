#pragma once
#include <ZGeom/ZGeom.h>

void computeHKMat(const ZGeom::EigenSystem& mhb, double timescale, ZGeom::DenseMatrixd& matHK);
void computeSGWMat(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW);
void computeSGWMat2(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW);

void computeGeometricLaplacianCoordinate(const CMesh& mesh, const MeshCoordinates& eCoord, MeshCoordinates& lCoord);


