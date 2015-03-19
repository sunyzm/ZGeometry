#pragma once
#include <ZGeom/ZGeom.h>
#include <vector>

MeshCoordinates least_square_fairing(CMesh& mesh, const ZGeom::MeshRegion& hole_region, int anchor_ring, double anchor_weight);
