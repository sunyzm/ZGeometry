#pragma once
#include <string>
#include <ZUtil/SimpleConfigLoader.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include "OutputHelper.h"

extern SimpleConfigLoader g_configMgr;
extern ZGeom::MatlabEngineWrapper g_engineWrapper;
extern OutputHelper qout;

enum GeometryTask {TASK_VIEWING = 0, TASK_REGISTRATION = 1, TASK_EDITING = 2};
extern GeometryTask g_task;

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
enum SignatureMode {
	SM_Normalized = 0, SM_LogNormalized, SM_MarkNegNormalized, 
	SM_AbsNormalized, SM_BandCurved, SM_PosNegPlot, 
	SM_CountSigModes
};
const std::string StrSignatureModeNames[] = {
	"Normalized", "Log-Normalized", "Mark-Neg", 
	"Abs-Normalized", "Band-Curved", "Pos-Neg Plot"
};

enum LaplacianType {
	Tutte = 0, Umbrella, NormalizedUmbrella, CotFormula, SymCot, 
	Anisotropic1, Anisotropic2, IsoApproximate, 
	LaplacianTypeCount
};

const std::string StrLaplacianTypes[] = {
	"Umbrella", "CotFormula", "Anisotropic1", "Anisotropic2"
};

struct GlobalSettings {
	GlobalSettings();
	
	int DEFAULT_EIGEN_SIZE;
	int DEFAULT_DEFORM_RING; 
	int LOAD_MHB_CACHE;
	double MIN_HK_TIMESCALE;
	double DEFUALT_HK_TIMESCALE;
	double MAX_HK_TIMESCALE;
	double PARAMETER_SLIDER_CENTER;
	double DR_THRESH_INCREMENT;
};
extern GlobalSettings gSettings;

/* mesh property names */
const std::string StrAttrColorUnnamed        = "vert_color_unnamed";
const std::string StrAttrColorGaussCurvature = "vert_color_gauss_curvature";
const std::string StrAttrColorMeanCurvature  = "vert_color_mean_curvature";
const std::string StrAttrColorEigenFunction  = "vert_color_eig_function";
const std::string StrAttrColorWaveletBasis   = "vert_color_wavelet_basis";
const std::string StrAttrColorDictAtom       = "vert_color_dict_atom";
const std::string StrAttrColorHKS			 = "vert_color_hks";
const std::string StrAttrColorHK			 = "vert_color_hk";
const std::string StrAttrColorBiharmonic	 = "vert_color_biharmonic";
const std::string StrAttrColorHeat			 = "vert_color_heat";
const std::string StrAttrColorMHWS			 = "vert_color_mhws";
const std::string StrAttrColorMHW			 = "vert_color_mhw";
const std::string StrAttrColorSGW			 = "vert_color_sgw";
const std::string StrAttrColorGeodesics	     = "vert_color_geodesics";
const std::string StrAttrColorPreset		 = "vert_color_preset";
const std::string StrAttrColorSimilarity     = "vert_color_similarity";
const std::string StrAttrColorPosDiff        = "vert_color_pos_diff";
const std::string StrAttrColorPartitions     = "vert_color_partitions";
const std::string StrAttrOriginalSignature   = "vert_original_signature";

const std::string StrAttrFeatureUnnamed      = "mesh_feature_unnamed";
const std::string StrAttrFeatureHKS			 = "mesh_feature_hks";
const std::string StrAttrFeatureMHWS		 = "mesh_feature_mhks";
const std::string StrAttrFeatureNeighbors    = "mesh_feature_neighbors";
const std::string StrAttrFeatureMultiHKS     = "mesh_feature_multi_hks";
const std::string StrAttrFeatureSparseSGW    = "mesh_feature_sparse_sgw";

const std::string StrAttrVectUnnamed         = "mesh_vector_unnamed";
const std::string StrAttrVecFaceNormal       = "mesh_vector_face_normals";
const std::string StrAttrVecVertNormal		 = "mesh_vector_vert_normals";

