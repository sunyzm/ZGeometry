#pragma once
#include <string>
#include <ZUtil/SimpleConfigLoader.h>
#include "OutputHelper.h"

enum GeometryTask {TASK_VIEWING = 0, TASK_REGISTRATION = 1, TASK_EDITING = 2};

enum DeformType {Simple, Gradient, Laplace, BiLaplace, Shell, SGW};

enum SignatureID {	SIGNATURE_ID = 0x0100, SIGNATURE_EIG_FUNC, SIGNATURE_HKS, SIGNATURE_HK, 
	SIGNATURE_MEAN_CURVATURE, SIGNATURE_GAUSS_CURVATURE, SIGNATURE_WKS, 
	SIGNATURE_MHWS, SIGNATURE_MHW, SIGNATURE_SGWS, SIGNATURE_SGW, 
	SIGNATURE_BIHARMONIC_DISTANCE, SIGNATURE_SIMILARITY_MAP, 
	SIGNATURE_ID_COUNT};

enum SignatureMode {Normalized = 0, LogNormalized, MarkNegNormalized, AbsNormalized, BandCurved, PosNegPlot, CountSigModes};
const std::string StrSignatureModeNames[] = {"Normalized", "Log-Normalized", "Mark-Neg", "Abs-Normalized", "Band-Curved", "Pos-Neg Plot"};

const std::string StrLaplacianTypes[4] = {"Umbrella", "CotFormula", "Anisotropic1", "Anisotropic2"};

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

const std::string StrOriginalSignature   = "vert_original_signature";

extern OutputHelper qout;
extern SimpleConfigLoader g_configMgr;
extern GeometryTask g_task;