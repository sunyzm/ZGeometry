#pragma once
#include <string>
#include <ZGeom/SimpleConfigLoader.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/ColorMap.h>
#include "OutputHelper.h"

extern SimpleConfigLoader g_configMgr;
extern ZGeom::MatlabEngineWrapper g_engineWrapper;
extern OutputHelper qout;

enum GeometryTask {TASK_VIEWING = 0, TASK_REGISTRATION = 1, TASK_EDITING = 2};
extern GeometryTask g_task;

extern double inpainting_error_curving_max;

enum DeformType {
	DEFORM_Simple, DEFORM_Gradient, DEFORM_Laplace, 
	DEFORM_BiLaplace, DEFORM_Shell, DEFORM_SGW,
	DEFORM_CountDeformTypes
};

enum SignatureID {	
	SIGNATURE_ID = 0x0100, SIGNATURE_EIG_FUNC, SIGNATURE_HKS, SIGNATURE_HK, 
	SIGNATURE_MEAN_CURVATURE, SIGNATURE_GAUSS_CURVATURE, SIGNATURE_WKS, 
	SIGNATURE_MHWS, SIGNATURE_MHW, SIGNATURE_SGWS, SIGNATURE_SGW, 
	SIGNATURE_BIHARMONIC_DISTANCE, SIGNATURE_SIMILARITY_MAP, 
	SIGNATURE_ID_COUNT
};

const std::string StrSignatureModeNames[] = {
	"Normalized", "Log-Normalized", "Mark-Neg", 
	"Abs-Normalized", "Band-Curved", "Pos-Neg Plot"
};

enum LaplacianType {
	Tutte = 0, Umbrella = 1, NormalizedUmbrella = 2, Geometric = 3, CotFormula = 4, 
	SymCot = 5, Anisotropic1 = 6, Anisotropic2 = 7, IsoApproximate = 8, 
	LaplacianTypeCount
};

const std::vector<std::string> StrLaplacianTypes {
	    "Umbrella", "CotFormula", "Anisotropic1", "Anisotropic2"
};

struct GlobalSettings {
	GlobalSettings();
	
    ZGeom::ColorMapType ACTIVE_COLOR_MAP_TYPE;
	int DEFAULT_EIGEN_SIZE;
	int DEFAULT_DEFORM_RING; 
	int LOAD_MHB_CACHE;
    int INPUT_SCALE_TO_UNIT;
	double MIN_HK_TIMESCALE;
	double DEFAULT_HK_TIMESCALE;
	double MAX_HK_TIMESCALE;
	double PARAMETER_SLIDER_CENTER;
	double DR_THRESH_INCREMENT;
};
extern GlobalSettings gSettings;

/* mesh property names */
const std::string StrAttrColorGaussCurvature    = "vert_color_gauss_curvature";
const std::string StrAttrColorMeanCurvature     = "vert_color_mean_curvature";
const std::string StrAttrColorEigenFunction     = "vert_color_eig_function";
const std::string StrAttrColorHKS			    = "vert_color_hks";
const std::string StrAttrColorHK			    = "vert_color_hk";
const std::string StrAttrColorBiharmonicField	= "vert_color_biharmonic_field";
const std::string StrAttrColorHeat			    = "vert_color_heat";
const std::string StrAttrColorMHWS			    = "vert_color_mhws";
const std::string StrAttrColorMHW			    = "vert_color_mhw";
const std::string StrAttrColorSGW[]             = {"vert_color_sgw1", "vert_color_sgw2", "vert_color_sgw3", "vert_color_sgw4", "vert_color_sgw5"};
const std::string StrAttrColorGeodesics	        = "vert_color_geodesics";
const std::string StrAttrColorPreset		    =  "vert_color_preset";
const std::string StrAttrColorSimilarity        = "vert_color_similarity";
const std::string StrAttrColorPosDiff           = "vert_color_pos_diff";
const std::string StrAttrColorPartitions        = "vert_color_partitions";
const std::string StrAttrOriginalSignature      = "vert_original_signature";
const std::string StrAttrColorInpaintError      = "vert_color_hole_vert_inpaint_error";

const std::string StrAttrFeatureUnnamed         = "mesh_feature_unnamed";
const std::string StrAttrFeatureHKS			    = "mesh_feature_hks";
const std::string StrAttrFeatureMHWS		    = "mesh_feature_mhks";
const std::string StrAttrFeatureNeighbors       = "mesh_feature_neighbors";
const std::string StrAttrFeatureMultiHKS        = "mesh_feature_multi_hks";
const std::string StrAttrFeatureSparseSGW       = "mesh_feature_sparse_sgw";
const std::string StrAttrFeatureSparseSGW2      = "mesh_feature_sparse_sgw2";
const std::string StrAttrFeatureMCA             = "mesh_feature_mca";
const std::string StrAttrLineFaceNormal         = "mesh_vector_face_normals";
const std::string StrAttrLineVertNormal		    = "mesh_vector_vert_normals";

const std::string StrAttrHoleSurroundingFaces   = "mesh_hole_surrounding_faces";
