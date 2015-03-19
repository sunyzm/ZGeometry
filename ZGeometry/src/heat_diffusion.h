#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/SparseSymMatVecSolver.h>
#include "MeshLaplacian.h"

void computeHeatDiffuseMatrix(CMesh& mesh, double t_multiplier, ZGeom::SparseSymMatVecSolver& heat_solver);
std::vector<double> calHeat(CMesh& mesh, int src_idx, ZGeom::SparseSymMatVecSolver& heat_solver);



