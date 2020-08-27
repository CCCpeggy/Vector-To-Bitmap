#include "../../Include/Common.h"

/**
 *
 * Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 * These parameters store values of parameters in the entire system
 * Don't modify unless you know what are they !
 * 
 * Actually some of these parameters can be modified during running time, 
 * please check CVUserInterface.ui by using Qt Designer
 *
*/

// Filename
std::string SystemParams::str_filename = "";


// FilePath
std::string SystemParams::str_Resources							= "./Resources";
std::string SystemParams::str_Resources_Original				= "./Resources/Original";
std::string SystemParams::str_Resources_Binarization			= "./Resources/Binarization";
std::string SystemParams::str_Resources_OCR						= "./Resources/OCR";
std::string SystemParams::str_Resources_CFR						= "./Resources/Cartoon Filter Region";
std::string SystemParams::str_Resources_CFC						= "./Resources/Cartoon Filter Color";
std::string SystemParams::str_Resources_TextureConfig			= "./Resources/TextureConfig/";
std::string SystemParams::str_Resources_Mesh					= "./Resources/Mesh/";
std::string SystemParams::str_Resources_KLMMesh					= "./Resources/KLMMesh/";
std::string SystemParams::str_Resources_Mesh_Parts				= "./Resources/Mesh_Parts/";
/**
 * Segmentation parameters
 * 
 */
int    SystemParams::s_km_n_clusters  = 3;	    // 3
int    SystemParams::s_km_max_iter    = 100;	// Maximum iteration to compute centroid of scribbles' gabor wavelet 100
double SystemParams::s_rescale_factor = 0.25;
double SystemParams::s_min_size_area  = 5000;


// Manjunath - Ma (Not used)
int SystemParams::g_scale		= 4;	// 4
int SystemParams::g_orientation = 6;	// 6
int SystemParams::g_total_dim	= 24;	//
int SystemParams::g_window		= 1;	// 3
double SystemParams::g_ul		= 0.2;	// 0.25
double SystemParams::g_uh		= 0.8;	// 0.75
int SystemParams::g_side		= 60;
int SystemParams::g_flag		= 0;

/**
 * Chan-Vese IPOL (Not used)
 * 
 */
double SystemParams::cv_ip_Tol     = 0.01;
int    SystemParams::cv_ip_MaxIter = 1000;
double SystemParams::cv_ip_Mu      = 0.01;	// mu, more accurately (smaller ? vs. producing a smoother boundary (larger ?.
double SystemParams::cv_ip_Nu      = 0.2;	// nu, sets the penalty
double SystemParams::cv_ip_Lambda1 = 1.0;	// fit weight inside the curve (default 1.0)
double SystemParams::cv_ip_Lambda2 = 1.0;	// fit weight outside the curve (default 1.0)
double SystemParams::cv_ip_dt      = 0.000001;

/**
 * Local binary pattern
 * 
 */
int SystemParams::lbp_neighbors = 16;
int SystemParams::lbp_radius = 4;

/**
 * Triangulation
 * 
 */
double SystemParams::t_smooth_factor = 0.5;				// for bezier interpolation (without reading cubic bezier from potrace)
double SystemParams::t_scale_factor = 1.0;				// scaling down the image for tracing and triangulation
double SystemParams::t_delaunay_aspect_bound = 0.125;	// factor for as-delaunay-as-possible (0.125) //0.0000001

// has issue with Triangulation2
double SystemParams::t_delaunay_max_length = 0 * SystemParams::t_scale_factor;		// set this to zero if you want to make the edges stretch as long as possible (50)

double SystemParams::t_STDis = 5;

/**
 * Tracing
 * 
 */
int    SystemParams::t_opticurve       = 1; //bool
double SystemParams::t_opttolerance    = 0.2; //0.2
double SystemParams::t_subdivide_limit = 0.5 * SystemParams::t_scale_factor;	// not more than 3
double SystemParams::t_rdp_epsilon     = 1.0 * SystemParams::t_scale_factor;	// not more than 3
int	   SystemParams::t_min_poly_area   = 15.0;	// 50
int	   SystemParams::t_rdp_point_min   = 3;

/**
 * Curvature Scale Space
 * 
 */
int SystemParams::css_sigma = 20;

/**
 * Screentone Renderer
 * 
 */
double SystemParams::sr_aa_value		= 0.5;
double SystemParams::sr_tone_frequency	= 1.0;
double SystemParams::sr_blur			= 0.00;
bool   SystemParams::sr_fixed			= false;
int	   SystemParams::sr_blackness		= 0;
double SystemParams::sr_tone_gap		= 1.0;
double SystemParams::sr_zoom_in			= 1.0;
float  SystemParams::zoom_in_value		= 1.0;

// Rendering
double SystemParams::cr_blackness	= 0.0 / 255.0;


/**
 * Cartoon+Texture Filter
 * 
 */
double SystemParams::ct_f_sigma = 5.5;
double SystemParams::ct_v_lim	= 20;

/**
 * Graphics UI
 * 
 */
bool SystemParams::gr_gpu					= false;
bool SystemParams::gr_sc_regions			= true;
bool SystemParams::gr_user_scribbles		= true;
bool SystemParams::gr_show_wireframe		= false;
bool SystemParams::gr_show_bbw_handles	    = true;
bool SystemParams::gr_show_bbw_pseudo_edge	= true;


/**
 * Part extraction
 * 
 */
int	   SystemParams::pe_num_iteration	= 1;
double SystemParams::pe_sample_distance = 2.5;
double SystemParams::pe_search_radius	= 5.0;
double SystemParams::pe_threshold		= 2.5;

//Otsu + Gaussian
int SystemParams::O_length = 11;
int SystemParams::G_block_size = 511;
int SystemParams::G_params = -17;
int SystemParams::OG_sigma = 0;
int SystemParams::OG_Bsigma = 0;
int SystemParams::OG_Wsigma = 0;


//Texture Unit
int SystemParams::TU_maskSize = 0;
int SystemParams::ST_height = 8000;
int SystemParams::ST_width = 5000;
int SystemParams::ST_groupType = 1;
int SystemParams::ST_gapX = 50;
int SystemParams::ST_gapY = 50;
int SystemParams::ST_orientation = 0;
float SystemParams::STG_widthP = 1;
float SystemParams::STG_heightP = 1;

SystemParams::SystemParams() { }
SystemParams::~SystemParams() { }

