
/**
 *
 * Some miscellaneous functions
 *
 * Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 */

#ifndef __Utility_Functions__
#define __Utility_Functions__

#define NOMINMAX

#include "../../Include/Common.h"

#define PT_TRI_INSIDE_ERROR 1e-8
#define SAMPLE_DENSE_PATH_LENGTH 5e-4

typedef std::pair<int, double> PairData;

class UtilityFunctions
{
public:

	enum 
	{
		COLINE = -2,
		PARALLEL = -1,
		CROSS_OUTSIDE = 0,
		CROSS_INSIDE = 1
	};

	//// Visualization
	// display image scaled (pointer array)
	static void DisplayImageDebug(double* img, int img_height, int img_width, double scale, std::string title);
	// display image scaled (opencv mat)
	static void DisplayImageDebug(cv::Mat img, double scale, std::string title);

	//// Data normalization
	static double* NormalizePhi(double* img, int img_height, int img_width);

	//// Color
	static void GetInterpolationColors(double d, int* r, int* g, int* b);
	
	//// String
	// split
	static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
	// split
	static std::vector<std::string> split(const std::string &s, char delim);
	
	//// Points
	// line with thickness
	static void UtilityFunctions::GetSegmentPoints(CVSystem::MyLine curLine, CVSystem::MyLine prevLine, CVSystem::MyLine nextLine, double t0, double t1, CVSystem::MyPoint* pA, CVSystem::MyPoint* pB, CVSystem::MyPoint* pC, CVSystem::MyPoint* pD);
	// bresenham rasterization
	static std::vector<CVSystem::MyPoint> Bresenham(float x1, float y1, float x2, float y2);

	//// Angle
	// degree to radian
	static double DegreeToRadian(double deg);
	// radian to degree
	static double RadianToDegree(double rad);
	// unsigned angle
	static double AngleInBetween(CVSystem::MyPoint pt1, CVSystem::MyPoint pt2);
	// signed angle
	static double GetRotation(CVSystem::MyPoint pt1, CVSystem::MyPoint pt2);
	// angle from x axis
	static double GetAngleFromXAxis(CVSystem::MyPoint pt);
	// what quadrant ?
	static int QuadrantFromRad(double rad);
	// rotate a point around center (0, 0)
	static CVSystem::MyPoint Rotate(CVSystem::MyPoint pt, double rad);
	
	//// Maximum value
	static int GetMaxValue(int* arrayData, int numData);

	//// Sorting
	static void QuickSort(std::vector<PairData>& items, int left, int right);	
	static void QuickSortMain(std::vector<PairData>& items);

	//// Curve and resample
	static double CurveLength(std::vector<CVSystem::MyPoint> curves);
	static void UniformResample(CVSystem::MyLine oriLine, std::vector<CVSystem::MyLine>& resampleLines, int N);
	static void UniformResample(CVSystem::MyLine oriLine, std::vector<CVSystem::MyPoint>& resampleCurves, int N);
	static void UniformResample(std::vector<CVSystem::MyPoint>& oriCurve, std::vector<CVSystem::MyPoint>& resampleCurve, int N /*bool isOpen*/ ) ;

	//// Distance to line
	static double DistanceToFiniteLine(CVSystem::MyPoint v, CVSystem::MyPoint w, CVSystem::MyPoint p);

	//// Curve orientation
	static bool IsClockwise(std::vector<CVSystem::MyPoint> polygon);
	static bool IsClockwise(std::vector<int> polygon, std::vector<CVSystem::MyPoint> vertexList);

	//// Data conversion
	// point set to line set
	static std::vector<CVSystem::MyLine>		PointsToLine(std::vector<CVSystem::MyPoint> points);
	// point set to line set
	static std::vector<CVSystem::MyIndexedLine> IndexedPointsToLine(std::vector<int> points);
	// Indexed triangle mesh no non-indexed triangle mesh
	static std::vector<CVSystem::MyTriangle>	ConvertIndexedTriangles(std::vector<CVSystem::MyIndexedTriangle> triangles, std::vector<CVSystem::MyPoint> vertexList);

	//// Indexing
	// check a line exists in data
	static bool DoesExist(std::vector<CVSystem::MyIndexedLine> data, CVSystem::MyIndexedLine line);
	// shift index 
	static void ShiftIndices(std::vector<CVSystem::MyIndexedTriangle>& idxTriangles, int num);
	// shift index 
	static void ShiftIndices(std::vector<CVSystem::MyIndexedBezierCurves>& bzCurves, int num);
	// convert triangle mesh to index list
	static std::vector<int> GetDistinctIndices(std::vector<CVSystem::MyIndexedTriangle> triangles, int numVertices);

	static cv::Point2f NextPoint(cv::Point2f p, double dis, double orientation);

	static bool InsideMask(cv::Mat mask, cv::Point2f p);

	static bool InsideImage(int imageWidth, int imageHeight, cv::Point2f p);

	static double Dis2f(cv::Point2f p1, cv::Point2f p2);

	static cv::Point2f intersection(cv::Point2f p1s, cv::Point2f p1e, cv::Point2f p2s, cv::Point2f p2e);

	static void ImageIntersectWithLinePoint(cv::Point &cp1, cv::Point &cp2, cv::Mat image, cv::Point p1, cv::Point p2);

	// signed angle
	static double GetRotation(cv::Point2f pt1, cv::Point2f pt2);
	// angle from x axis
	static double GetAngleFromXAxis(cv::Point2f pt);

	static cv::Point2d GetCubicPt(cv::Point2d p1, cv::Point2d p2, cv::Point2d p3, cv::Point2d p4, double t );

	static cv::Point2d GetCubicDerivative( cv::Point2d p1, cv::Point2d p2, cv::Point2d p3, cv::Point2d p4, double t );

	static void splitBezierCubic(cv::Point2d orgCP[4], cv::Point2d newCPA[4], cv::Point2d newCPB[4], double t);

	static void splitBezierCubic( cv::Point2d orgCP[4], cv::Point2d resCP[4], double t1, double t2 );

	static double determinant_3x3( double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21, double m22 );

	static int PointInLine2D( const cv::Point2d& p, const cv::Point2d& Lp1, const cv::Point2d& Lp2 );

	static double getLinePara( const cv::Point2d &p, const cv::Point2d &Lp1, const cv::Point2d &Lp2 );

	static void GetBaryCoord( cv::Point3d &bc, cv::Point2d &p, cv::Point2d &p1, cv::Point2d &p2, cv::Point2d &p3 );

	static inline double signd(const double i) { return i<0 ? -1:1;}

	static inline double Length(cv::Point p) { return sqrt((double)(p.x * p.x + p.y * p.y));}

	static inline double Length(cv::Point p1, cv::Point p2) { return sqrt((double)((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y)));}

	static inline bool pointEqual( const cv::Point2d& p1, const cv::Point2d& p2, double threshold ) { return fabs(p1.x-p2.x)+fabs(p1.y-p2.y) < threshold; }

	static inline bool pointEqual( const double p1x, const double p1y, const double p2x, const double p2y, double threshold ) { return fabs(p1x-p2x)+fabs(p1y-p2y) < threshold; }

	static int LineIntersect( cv::Point2d &cross_p, cv::Point2d &L1_p1, cv::Point2d &L1_p2, cv::Point2d &L2_p1, cv::Point2d &L2_p2 );

	static bool PointTriangleInside(cv::Point3d &w, cv::Point2d &pos, cv::Point2d &tp0, cv::Point2d &tp1, cv::Point2d &tp2);

	static bool PointTriangleInside(cv::Point3d &w, cv::Point2d pos, cv::Point2d tp[3]);

	static bool CubicToLineApx( cv::Point2d cp[4] );

	static void CubicRoots( double t[3], double p[4] );

	static double PointLineDistance( cv::Point2d& p, cv::Point2d& L_p1, cv::Point2d& L_p2 );

	static bool PointLineInside( cv::Point2d point, cv::Point2d Line_p1, cv::Point2d Line_p2 );

	static double GetLinePara( cv::Point2d p, cv::Point2d Lp1, cv::Point2d Lp2 );

	static bool AabbIntersect( const cv::Point2d p11, const cv::Point2d p12, const cv::Point2d p21, const cv::Point2d p22 );

};

#endif
