#pragma once
#include <string>
#include <ZUtil/SimpleConfigLoader.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include "OutputHelper.h"

enum GeometryTask {TASK_VIEWING = 0, TASK_REGISTRATION = 1, TASK_EDITING = 2};
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

const std::string StrColorGaussCurvature = "vert_color_gauss_curvature";
const std::string StrColorMeanCurvature  = "vert_color_mean_curvature";
const std::string StrColorEigenFunction  = "vert_color_eig_function";
const std::string StrColorWaveletBasis   = "vert_color_wavelet_basis";
const std::string StrColorDictAtom       = "vert_color_dict_atom";
const std::string StrColorHKS			 = "vert_color_hks";
const std::string StrColorHK			 = "vert_color_hk";
const std::string StrColorBiharmonic	 = "vert_color_biharmonic";
const std::string StrColorHeat			 = "vert_color_heat";
const std::string StrColorMHWS			 = "vert_color_mhws";
const std::string StrColorMHW			 = "vert_color_mhw";
const std::string StrColorSGW			 = "vert_color_sgw";
const std::string StrColorGeodesics	     = "vert_color_geodesics";
const std::string StrColorPreset		 = "vert_color_preset";
const std::string StrColorSimilarity     = "vert_color_similarity";
const std::string StrColorPosDiff        = "vert_color_pos_diff";
const std::string StrColorPartitions     = "vert_color_partitions";
const std::string StrOriginalSignature   = "vert_original_signature";

extern OutputHelper qout;
extern SimpleConfigLoader g_configMgr;
extern GeometryTask g_task;
extern ZGeom::MatlabEngineWrapper g_engineWrapper;

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
