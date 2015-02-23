#pragma once
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <engine.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include <ZGeom/Mesh.h>
#include "MeshLaplacian.h"
#include "global.h"

void computeKernelSignatureFeatures(CMesh& mesh, const ZGeom::EigenSystem& es, const std::vector<double>& timescales, std::function<double(double, double)> genfunc, std::string featureStr);

void computeSGWMat(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW);
void computeSGWMat2(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::DenseMatrixd& matSGW);

void computeGeometricLaplacianCoordinate(const CMesh& mesh, const MeshCoordinates& eCoord, MeshCoordinates& lCoord);



