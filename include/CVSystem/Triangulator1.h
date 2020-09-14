
/*
 * The most important class for line tracing and constrained triangulation
 * The input is raster image and Potrace is used to perform line tracing.
 * 
 *
 * Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 */

#ifndef __Triangulator_1__
#define __Triangulator_1__

#include "../../Include/Common.h"
#include "CVSystem/CGALTools.h"

namespace CVSystem
{
	struct PotraceInfo;
	struct pathInfo;
	struct curveInfo;

	enum curveType
	{
		CURVE,
		CORNER,
		LINE
	};

	enum PathSign
	{
		UNSET,
		Clockwise,			//順時針
		Counterclockwise	//逆時針
	};

	struct curveInfo
	{
		pathInfo *parent;

		curveInfo *pre;
		curveInfo *next;

		curveType cType;
		pathInfo *oppsite_path;
		int ft_id;
		bool invColor;
		bool covered;

		cv::Point2d points[5];//錨點 控制點 控制點 錨點 交點 normalize
		cv::Point2d oriPoints[5];//錨點 控制點 控制點 錨點 交點
		cv::Point3d KLMs[4];
		cv::Point3d oriKLMs[4];

		cv::Point2d bound_max, bound_min;
		std::vector<cv::Point2d> dPts;
		std::vector<cv::Point2d> crossPos;
		std::vector<int> crossIDs;
		std::vector<double> crossT;

		curveInfo();

		void Copy(curveInfo cInfo);

		void OverrideProperty(curveInfo *cInfo);

		void Reset(cv::Point2d _cp[4]);

		bool Compare(curveInfo CI);

		void DenseSampling();

		void SetBound();

		bool aabbIntersects(curveInfo cInfo);

		void Split( std::vector<curveInfo*> &cInfos );

		void CubicSplit( std::vector<curveInfo*> &cInfos );

		void InsertBack(curveInfo* newNext);

		void InsertFront(curveInfo* newPre);
	};

	struct pathInfo
	{
		PotraceInfo *Graph;
		pathInfo *pre;
		pathInfo *next;
		PathSign pSign;
		int grayScale;
		curveInfo *CI_head;
		curveInfo *CI_tail;
		int curveCount;
		int ColorId; //0 black 1 white 2~X ST

		pathInfo();

		void AddCurve(curveInfo *cInfo);

		int RemoveCurve(curveInfo *cInfo);
	};

	struct PotraceInfo
	{
		int width;
		int height;
		cv::Mat image;

		pathInfo *PI_head;
		pathInfo *PI_tail;
		int pathCount;

		std::vector<curveInfo *> allCurve;
		std::vector<curveInfo *> outerCurve;//最外圍的curve

		PotraceInfo();

		void AddPath(pathInfo *pInfo);

		bool CheckOuterVertex(MyTriangle tri);
	};

	struct stCurveBountLink
	{
		curveInfo *oriCI;
		cv::Point2d cp[4];
		cv::Point2d bound[4];
		int n_bound_point;	// 紀錄使用幾個點做boundary
		stCurveBountLink* nxt;

		stCurveBountLink():oriCI(NULL){}
		stCurveBountLink(curveInfo *CI, cv::Point2d cp0, cv::Point2d cp1);
		stCurveBountLink(curveInfo *CI, cv::Point2d _cp[4]);

		void Reset(cv::Point2d _cp[4]);

		int Intersect( stCurveBountLink* cb ); //return 0:no intersect, 1:intersect, 2:not intersect but connect

		int Intersect( cv::Point2d p1, cv::Point2d p2 ); // intersect with line, return 0:no intersect, 1:intersect, 2:not intersect but connect

		void Intersect( std::vector<double> &cross_para, std::vector<cv::Point2d> &clip_inside, cv::Point2d tri[3] ); // 和三角形的intersection，紀錄和三角形的所有交錯點和切在內部的的線段，參考考I-1

		bool CheckTriangleInside( cv::Point2d tri[3] );
	};

	class Triangulator1
	{
	public:
		Triangulator1();
		~Triangulator1();

		std::vector<double*> debugLines;
		void TraceImage(std::string strFilename, cv::Mat img, int* mask, int* dilatedMask, bool isLog);	// Trace
		
		// void OCalculate();		// Original (Deprecated)
		// void LSCalculate1();	// Least Square (Deprecated)
		void LSCalculate2(cv::Mat img, bool isLog, float offsetX = 0, float offsetY = 0);	// Complete Triangulation
		// void LSCalculate3(cv::Mat img, bool isLog, float offsetX = 0, float offsetY = 0);	// Complete Triangulation
		// void LSCalculateTU(cv::Mat img, bool isLog);	// Complete Triangulation
		
		//KLM
	//	void GetCurveInfo(pathInfo *pInfo, potrace_path_t*  path, int w, int h);
	//	PotraceInfo* parsePngPath(cv::Mat image);
	//	void BackPathAndCoveredTest(CVSystem::PotraceInfo *graph);
	//	int CalculateCubicParam( cv::Point3d klm[4], cv::Point2d cp[4] );
	//	void InsideOutsideColorCheck(PotraceInfo *graph);
	//	void InsideOutsideColorCheck(curveInfo *cInfo);
	//	bool InsideOutsideColorCheck(cv::Point2d cp[4], cv::Point3d klm[4], PathSign pSign);
	//	int PathInsideTest( CVSystem::pathInfo *pInfo , cv::Point2d pos );
	//	void CalTriangleParameter( cv::Vec3d klm[3], cv::Vec4d texCoord[3], CVSystem::curveInfo *cInfo, cv::Point2d cp[4], cv::Point2d tp[3] );
	//	void CubicCubicClipping( CVSystem::curveInfo* cInfo_1, CVSystem::curveInfo* cInfo_2);
	//	bool LineCurveIntersect( double t[3], cv::Point2d cp0, cv::Point2d cp1, cv::Point2d cp2, cv::Point2d cp3, cv::Point2d Lp0, cv::Point2d Lp1 );
	//	void KLMCalculate(cv::Mat img, bool isLog, float offsetX = 0, float offsetY = 0);	// Complete Triangulation

	//	// delete these three functions
	//	std::vector<MyTriangle> GetOTriangles();	// Original
	//	std::vector<MyQuad>		GetOQuads();		// Original
	//	std::vector<MyTriangle> GetLSTriangles();	// LS
	//	
		std::vector<MyQuad> GetLSQuads();	// LS (Least Square Fitting for Cubic Bezier)

	//	// Indexed Format
	//	
		std::vector<CVSystem::MyPoint>	   GetVertexList();
		std::vector<MyIndexedTriangle>	   GetIndexedTriangles();
		std::vector<MyIndexedTriangle>	   GetWTriangles();
		std::vector<MyIndexedTriangle>	   GetSCTriangles();
		std::vector<MyIndexedTriangle>	   GetBTriangles();
		std::vector<MyIndexedTriangle>	   GetBorderWTriangles();
		std::vector<MyIndexedTriangle>	   GetBorderSCTriangles();
		std::vector<MyIndexedBezierCurves> GetIndexedBezierCurves();

	//	std::vector<KLMTriangle>		   GetIndexedKLMTriangles(){return _indexedKLMTriangles;}

	//	// Corner
		std::vector<CVSystem::MyPoint> GetCornerList();	 // Delete this
		std::vector<int> GetCornerIndices();
		std::vector<CVSystem::MyIndexedLine> GetCornerEdges();

	//	std::vector<int> GetPartOffset() { return _offsets; }

	//	// scaled width
	//	int GetWidthScaled() { return _w_scaled; }

	//	// scaled height
	//	int GetHeightScaled(){ return _h_scaled; }

	//	cv::Mat GetTriImage(){ return _TriImage; }

	//public:
	//	// read .myobj file
	//	void ReadFromMYOBJ(std::string filename);

	//	// read .myobj file
	//	void ReadFromMYOBJasSTG(std::string filename);

	//	TriangulationInfo* STGInfo;

	//private:
		// stuff from Potrace
		potrace_bitmap_t* GetBM(cv::Mat img);
		void CurveToBezier(MyPoint p0, MyPoint p1, MyPoint p2, MyPoint p3, MyPoint& cp0, MyPoint& cp1);

		// create CDT		
		CVSystem::CD_Cdt GetCDT(std::vector<MyPoint> pPoly, std::vector<std::vector<MyPoint>> cPoly, bool shouldRefine = true, bool isST = false);
		
		// create CDT
		// CVSystem::CD_Cdt GetCDTWithRectangleBorder(std::vector<std::vector<MyPoint>> polygons);

		// create ST CDT 
		CVSystem::CD_Cdt GetSTCDTWithRectangleBorder(std::vector<std::vector<MyPoint>> polygons, std::vector<std::vector<MyPoint>> STpolygons);

		// Get hole seed
		std::list<CVSystem::CD_Point> GetHoleSeeds(std::vector<MyPoint> pPoly, std::vector<std::vector<MyPoint>> cPoly, bool isST = false);
		
		// Assign vertex index
		void AssignVertexIndex(CVSystem::CD_Cdt& cdt);

		// create CDT		
		// CVSystem::CD_Cdt GetKMLCDT(std::vector<MyPoint> pPoly, std::vector<std::vector<MyPoint>> cPoly, bool shouldRefine = true);

		// create CDT
		// CVSystem::CD_Cdt GetKMLCDT(std::vector<std::vector<MyPoint>> polygons);

		// Original Triangulation (Deprecated)
		std::vector<MyPoint> GetPolygonO(potrace_curve_t curve);
		std::vector<MyQuad>  GetBezierQuad0( potrace_path_t*  path);

		// Get polygon from Potrace
		std::vector<MyPoint> GetPolygonLS(potrace_path_t*  path, bool isDraw = false);

		// Bezier curve fitting
		std::vector<MyQuad>  GetBezierQuadLS(std::vector<CVSystem::MyPoint> poly, std::vector<bool> rdpFlags, bool isInside);
				
		// Uniform sampling
		void Resubdivide(std::vector<CVSystem::MyPoint>& poly);

		// Clean up small polygon
		void CleanUpPoly(std::vector<CVSystem::MyPoint>& poly);

		// RDP partitioning
		void RDPSimplification(std::vector<CVSystem::MyPoint> poly, std::vector<MyPoint>& rdpPoly, std::vector<bool>& bFlag);

		void RDPSimplificationST(std::vector<std::vector<MyPoint>> polygons, std::vector<CVSystem::MyPoint> poly, std::vector<MyPoint>& rdpPoly, std::vector<bool>& bFlag);

		int ReduceSTpoint(std::vector<std::vector<MyPoint>> polygons, std::vector<CVSystem::MyPoint> souPoly, std::vector<CVSystem::MyPoint> &desPoly);
		
		// Build bezier curve and triangle mesh
		void ProcessEdgesFromCDT(CVSystem::CD_Cdt cd_cdt);

		// Corner Detection
		void GetCornersFromPolygon();

		// Corner Detection
		void ProcessPolygonCorners();

		void DiscardVertexOffset(float offsetX, float offsetY);

	private:
		std::string _strFilename;

		std::vector<MyTriangle> _oTriangles;	// deprecated
		std::vector<MyQuad>		_oQuads;		// deprecated				
		std::vector<MyTriangle> _lsTriangles;	// deprecated

		std::vector<MyQuad>		_lsQuads;		// don't delete this, it contains original LS fitting

		// New format !
		std::vector<int>				   _offsets;			 // .myobj offsets
		std::vector<MyPoint>			   _vertexList;			 // The actual vertices
		std::vector<MyIndexedTriangle>     _indexedTriangles;	 // All triangles
		std::vector<MyIndexedTriangle>     _borderSCTriangles;	 // Additional screentone triangles
		std::vector<MyIndexedTriangle>     _borderWTriangles;	 // Additional white triangles
		std::vector<MyIndexedTriangle>     _scTriangles;		 // Screentone triangles
		std::vector<MyIndexedTriangle>     _wTriangles;			 // White triangles
		std::vector<MyIndexedTriangle>     _bTriangles;			 // Black triangles
		std::vector<MyIndexedBezierCurves> _indexedBezierCurves; // Bezier curves (start and end points only)
		
		std::vector<CVSystem::MyPoint>		 _cornerList;		// Delete this
		std::vector<int>					 _cornerIndices;	// this isn't sensible
		std::vector<CVSystem::MyIndexedLine> _cornerEdges;

		//KLM
		std::vector<KLMTriangle>		_indexedKLMTriangles;
		
		cv::Mat _mask;						// screentone mask
		cv::Mat _dilated_mask;				// dilated screentone mask
		cv::Mat _distance_img;				// unused

		cv::Mat _mask_scaled;
		cv::Mat _dilated_mask_scaled;
		cv::Mat _distance_img_scaled;

		cv::Mat _TriImage;					//Image in MyOBJ

		bool _has_tone;	// does mesh have screentone?

		int _w_scaled;	// scaled image width
		int _h_scaled;	// scaled image height

		//// Gaussian refiner
		CSSSmoothing*	 _cssSmoothing;

		//// Potrace
		potrace_state_t* _trace_state;
		potrace_state_t* _trace_stateST;

		cv::Mat STimg;

		cv::Mat testImg;
		int testcount;

		cv::Mat traceImg, traceImgS;
	};
}
#endif

