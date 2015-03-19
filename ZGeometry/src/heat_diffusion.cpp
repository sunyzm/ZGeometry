#include "heat_diffusion.h"
#include <ZGeom/MatVecArithmetic.h>
#include <algorithm>

void computeHeatDiffuseMatrix(CMesh& mesh, double t_multiplier, ZGeom::SparseSymMatVecSolver& heat_solver)
{
    const int vertCount = mesh.vertCount();
    const double t = std::pow(mesh.getAvgEdgeLength(), 2) * t_multiplier;

    MeshLaplacian laplacian; 
    laplacian.constructCotFormula(&mesh);
    const ZGeom::SparseMatrix<double>& matW = laplacian.getW();
    const ZGeom::SparseMatrix<double>& matLc = laplacian.getLS();	// negative
    
    ZGeom::SparseMatrix<double> heat_diffuse_mat;
    ZGeom::addMatMat(matW, matLc, -t, heat_diffuse_mat);	//A = W - t*Lc
    
    heat_solver.initialize(heat_diffuse_mat, true, true);
}

std::vector<double> calHeat(CMesh& mesh, int src_idx, ZGeom::SparseSymMatVecSolver& heat_solver)
{
    const int vertCount = mesh.vertCount();
    std::vector<double> vHeat(vertCount);
    ZGeom::VecNd vInitHeat(vertCount, 0);
    vInitHeat[src_idx] = 1.0;
    ZGeom::VecNd vSolvedHeat(vertCount, 0);

    heat_solver.solve(1, vInitHeat.c_ptr(), vSolvedHeat.c_ptr());
    std::copy_n(vSolvedHeat.c_ptr(), vertCount, vHeat.begin());
    return vHeat;
}