
/**
 *
 * Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 */

#include "../../Include/Common.h"

double* UtilityFunctions::NormalizePhi(double* img, int img_height, int img_width)
{
	double* normImg = (double*)malloc(img_width * img_height * sizeof(double));

	double minVal = DBL_MAX; double maxVal = DBL_MIN;
	for(int i = 0; i < img_width; i++)
	{
		for(int j = 0; j < img_height; j++)
		{
			double d = img[i + j * img_width];
			if(d < minVal) minVal = d;
			if(d > maxVal) maxVal = d;
		}
	}

	double maxMin = maxVal - minVal;
	for( int j = 0; j < img_height; j++ )
		for( int i = 0; i < img_width; i++ )
			normImg[i + j * img_width] = ( (img[i + j * img_width] - minVal) / maxMin) - 0.5;

	return normImg;
}

void UtilityFunctions::DisplayImageDebug(double* img, int img_height, int img_width, double scale, std::string title)
{
	double minVal = DBL_MAX;
	double maxVal = DBL_MIN;
	for(int i = 0; i < img_width; i++)
	{
		for(int j = 0; j < img_height; j++)
		{
			double d = img[i + j * img_width];
			if(d < minVal) minVal = d;
			if(d > maxVal) maxVal = d;
		}
	}

	cv::Mat drawing = cv::Mat::zeros( cv::Size(img_width, img_height), CV_8UC3 );

	for( int j = 0; j < img_height; j++ )
	{ 
		for( int i = 0; i < img_width; i++ )
		{
			int idx = i + j * img_width;
			double d = img[i + j * img_width];

			int _b, _g, _r;
			GetInterpolationColors(d, &_b, &_g, &_r);
			drawing.at<cv::Vec3b>(j,i) = cv::Vec3b(_b, _g, _r);
		}
	}
		
	//cv::imwrite(title + ".png", drawing);
	cv::resize(drawing, drawing, cv::Size(img_width * scale, img_height * scale));
	cv::imshow(title, drawing);
}

void UtilityFunctions::DisplayImageDebug(cv::Mat img, double scale, std::string title)
{
	int img_height = img.rows * scale;
	int img_width = img.cols * scale;

	cv::Mat dispImg = cv::Mat::zeros(img_height, img_width, CV_8UC3);
	cv::Mat scaledImg = cv::Mat::zeros(img_height, img_width, CV_32FC1);
	cv::resize(img, scaledImg, scaledImg.size());

	cv::cvtColor(dispImg, dispImg, CV_BGR2RGB);
	cv::normalize(scaledImg, scaledImg, 0.0, 1.0, cv::NORM_MINMAX, CV_32FC1);

	for(int a = 0; a < scaledImg.cols; a++)
	{
		for(int b = 0; b < scaledImg.rows; b++)
		{
			int _b, _g, _r;
			GetInterpolationColors(scaledImg.ptr<float>(b, a)[0], &_b, &_g, &_r);
			dispImg.at<cv::Vec3b>(b, a) = cv::Vec3b(_b, _g, _r);
		}
	}
	/*
	cv::Mat inpImg;
	scaledImg.convertTo(inpImg, CV_32FC1);

	double minVal = DBL_MAX;
	double maxVal = DBL_MIN;
	for(int i = 0; i < img_width; i++)
	{
		for(int j = 0; j < img_height; j++)
		{
			double d = inpImg.ptr<float>(j, i)[0];
			if(d < minVal) minVal = d;
			if(d > maxVal) maxVal = d;
		}
	}

	cv::Mat drawing = cv::Mat::zeros( cv::Size(img_width, img_height), CV_8UC3 );
	for( int j = 0; j < img_height; j++ )
	{ 
		for( int i = 0; i < img_width; i++ )
		{
			double d = inpImg.ptr<float>(j, i)[0];
			d = (d - minVal) / (maxVal - minVal);
			drawing.at<cv::Vec3b>(j,i) = cv::Vec3b(0, 0, 255 * d);
		}
	}
	cv::imwrite(title + ".png", drawing);
	*/

	cv::imshow(title, dispImg);
	cv::imwrite(title + ".png", dispImg);
}

void UtilityFunctions::GetInterpolationColors(double d, int* b, int* g, int* r)
{
	int blue = 0;
	int green = 0;
	int red = 0;

	if (d < 0.5)
	{
		blue   = (255 * d * 2.0) +  0   * (0.5 - d) * 2.0;
		green  = (255 * d * 2.0) +  0   * (0.5 - d) * 2.0;
		red    = (0   * d * 2.0) +  255 * (0.5 - d) * 2.0;
	}
	else
	{
		blue   = 255 * (d - 0.5) * 2.0 + 255 * (1.0 - d) * 2.0;
		green  = 0   * (d - 0.5) * 2.0 + 255 * (1.0 - d) * 2.0;  
		red    = 0   * (d - 0.5) * 2.0 + 0   * (1.0 - d) * 2.0;  
	}

	(*b) = blue;
	(*g) = green;
	(*r) = red;
}

int UtilityFunctions::GetMaxValue(int* arrayData, int numData)
{
	int maxVal = INT_MIN;
	for(int a = 0; a < numData; a++)
		if(arrayData[a] > maxVal) maxVal = arrayData[a];
	return maxVal;
}

std::vector<std::string>& UtilityFunctions::split(const std::string &s, char delim, std::vector<std::string> &elems) 
{
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) 
	{
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> UtilityFunctions::split(const std::string &s, char delim) 
{
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}


std::vector<CVSystem::MyPoint> UtilityFunctions::Bresenham(float x1, float y1, float x2, float y2)
{
	std::vector<CVSystem::MyPoint> pts;

	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = abs(delta_y) << 1;

	pts.push_back(CVSystem::MyPoint(x1, y1));

	if (delta_x >= delta_y)
	{	
		int error(delta_y - (delta_x >> 1));	// error may go below zero
		while (x1 != x2)
		{
			if ((error >= 0) && (error || (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}	// else do nothing
			error += delta_y;
			x1 += ix;
			pts.push_back(CVSystem::MyPoint(x1, y1));
		}
	}
	else
	{	
		int error(delta_x - (delta_y >> 1));	// error may go below zero
		while (y1 != y2)
		{
			if ((error >= 0) && (error || (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}	// else do nothing
			error += delta_x;
			y1 += iy;
			pts.push_back(CVSystem::MyPoint(x1, y1));
		}
	}
	return pts;
}


void UtilityFunctions::QuickSort(std::vector<PairData>& items, int left, int right)
{
	int i, j;
	PairData x, y;

	i = left;
	j = right;

	x = items[(left + right) / 2];

	do {
		while((items[i].second < x.second) && (i < right))
			i++;

		while((x.second < items[j].second) && (j > left))
			j--;

		if(i <= j) 
		{
			y = items[i];
			items[i] = items[j];
			items[j]= y;

			i++; j--;
		}
	} while(i <= j);

	if(i < right)
		QuickSort(items, i, right);
	if(left < j)
		QuickSort(items, left, j);

}

// Sorting pair

void UtilityFunctions::QuickSortMain(std::vector<PairData>& items)//min~max
{
	QuickSort(items, 0, items.size() - 1);
}

void UtilityFunctions::GetSegmentPoints(CVSystem::MyLine curLine, CVSystem::MyLine prevLine, CVSystem::MyLine nextLine, double t0, double t1, CVSystem::MyPoint* pA, CVSystem::MyPoint* pB, CVSystem::MyPoint* pC, CVSystem::MyPoint* pD)
{
	// point
	CVSystem::MyPoint p0 = curLine.GetPointA();
	CVSystem::MyPoint p1 = curLine.GetPointB();

	// direction
	CVSystem::MyPoint prevDirNorm  = prevLine.Direction().Norm();
	CVSystem::MyPoint curDirNorm   = curLine.Direction().Norm();
	CVSystem::MyPoint nextDirNorm  = nextLine.Direction().Norm();

	CVSystem::MyPoint d0;
	CVSystem::MyPoint d1;			

	if(!prevLine.Invalid() && !nextLine.Invalid())		// normal drawing
	{
		d0 = (prevDirNorm + curDirNorm) / 2.0;
		d1 = (nextDirNorm + curDirNorm) / 2.0;
	}
	else if(prevLine.Invalid() && nextLine.Invalid())	// line consists of only one line segment
	{
		d0 = curDirNorm;
		d1 = curDirNorm;
	}
	else if(prevLine.Invalid())							// draw starting segment
	{
		d0 = curDirNorm;
		d1 = (nextDirNorm + curDirNorm) / 2.0;
	}
	else if(nextLine.Invalid())							// draw ending segment
	{
		d0 = (prevDirNorm + curDirNorm) / 2.0;
		d1 = curDirNorm;
	}

	// thickness
	float p0HalfThickness = t0;
	float p1HalfThickness = t1;

	if(p0HalfThickness <= 0.0) p0HalfThickness = 0.5;
	if(p1HalfThickness <= 0.0) p1HalfThickness = 0.5;

	d0 = d0.Norm();
	d1 = d1.Norm();

	d0 *= p0HalfThickness;
	d1 *= p1HalfThickness;

	CVSystem::MyPoint d0Left (-d0.y,  d0.x);
	CVSystem::MyPoint d0Right( d0.y, -d0.x);
	CVSystem::MyPoint d1Left (-d1.y,  d1.x);
	CVSystem::MyPoint d1Right( d1.y, -d1.x);

	*pA = p0 + d0Left;
	*pB = p0 + d0Right;
	*pC = p1 + d1Left;
	*pD = p1 + d1Right;
}


double UtilityFunctions::DegreeToRadian(double deg)
{
	return ((double)deg) * M_PI / 180.0; 
}

double UtilityFunctions::RadianToDegree(double rad)
{
	return rad * 180.0 / M_PI;
}

double UtilityFunctions::AngleInBetween(CVSystem::MyPoint pt1, CVSystem::MyPoint pt2)
{
	return acos(pt1.Dot(pt2) / (pt1.Length() * pt2.Length()));
}

double UtilityFunctions::GetRotation(CVSystem::MyPoint pt1, CVSystem::MyPoint pt2)
{
	using namespace CVSystem;
	double perpDot = pt1.x * pt2.y - pt1.y * pt2.x;
	return (float)atan2(perpDot, pt1.Dot(pt2));
}

double UtilityFunctions::GetAngleFromXAxis(CVSystem::MyPoint pt)
{
	using namespace CVSystem;
	return GetRotation(MyPoint(1, 0), pt);
}

int  UtilityFunctions::QuadrantFromRad(double rad)
{
	double half_pi = M_PI / 2.0;
	if(rad >= 0 && rad < half_pi)
	{
		return 0;
	}
	else if(rad >= half_pi && rad < M_PI)
	{
		return 1;
	}
	else if(rad >= M_PI && rad < (M_PI + half_pi))
	{
		return 2;
	}
	return 3;
}

CVSystem::MyPoint UtilityFunctions::Rotate(CVSystem::MyPoint pt, double rad)
{
	double cs = cos(rad);
	double sn = sin(rad);

	double x = pt.x * cs - pt.y * sn;
	double y = pt.x * sn + pt.y * cs;

	return CVSystem::MyPoint(x, y);
}

// Todo: implement closed curves
double UtilityFunctions::CurveLength(std::vector<CVSystem::MyPoint> curves)
{
	double length = 0.0;
	for(size_t a = 1; a < curves.size(); a++) { length += curves[a].Distance(curves[a-1]); }
	return length;
}

void UtilityFunctions::UniformResample(CVSystem::MyLine oriLine, std::vector<CVSystem::MyLine>& resampleLines, int N)
{
	using namespace CVSystem;
	std::vector<MyPoint> oriCurve;
	std::vector<CVSystem::MyPoint> resampleCurve;
	oriCurve.push_back(MyPoint(oriLine.XA, oriLine.YA));
	oriCurve.push_back(MyPoint(oriLine.XB, oriLine.YB));
	UniformResample(oriCurve, resampleCurve, N);

	for(size_t a = 1; a < resampleCurve.size(); a++)
		{ resampleLines.push_back(MyLine(resampleCurve[a-1].x, resampleCurve[a-1].y, resampleCurve[a].x, resampleCurve[a].y)); }
}

void UtilityFunctions::UniformResample(CVSystem::MyLine oriLine, std::vector<CVSystem::MyPoint>& resampleCurves, int N)
{
	using namespace CVSystem;
	std::vector<MyPoint> oriCurve;
	oriCurve.push_back(MyPoint(oriLine.XA, oriLine.YA));
	oriCurve.push_back(MyPoint(oriLine.XB, oriLine.YB));
	UniformResample(oriCurve, resampleCurves, N);
}

void UtilityFunctions::UniformResample(std::vector<CVSystem::MyPoint>& oriCurve, std::vector<CVSystem::MyPoint>& resampleCurve, int N /*bool isOpen*/ ) 
{
	using namespace CVSystem;

	resampleCurve.clear();
	for(int a = 0; a <= N; a++) { resampleCurve.push_back(MyPoint(0,0)); }
	resampleCurve[0].x = oriCurve[0].x;
	resampleCurve[0].y = oriCurve[0].y;
	resampleCurve[N].x = oriCurve[oriCurve.size()-1].x;
	resampleCurve[N].y = oriCurve[oriCurve.size()-1].y;

	double pl_length = CurveLength(oriCurve);

	double resample_size = pl_length / (double) N;
	int curr = 0;
	double dist = 0.0;
	for (int i = 1; i < N; ) 
	{   
		double last_dist = oriCurve[curr].Distance(oriCurve[curr+1]);
		
		dist += last_dist;

		if (dist >= resample_size) 
		{
			//put a point on line
			double _d = last_dist - (dist-resample_size);
			MyPoint cp(oriCurve[curr].x, oriCurve[curr].y); 
			MyPoint cp1(oriCurve[curr+1].x, oriCurve[curr+1].y);

			MyPoint dirv = cp1-cp; 
			dirv = dirv * (1.0 / dirv.Length());
			
			resampleCurve[i] = cp + dirv * _d;
			i++;

			dist = last_dist - _d; //remaining dist         

			//if remaining dist to next point needs more sampling... (within some epsilon)
			while (dist - resample_size > 1e-3) 
			{
				resampleCurve[i] = resampleCurve[i-1] + dirv * resample_size;
				dist -= resample_size;
				i++;
			}
		}
		curr++;
	}

	/*for(size_t a = 1; a < resampleCurve.size(); a++)
	{
		std::cout << resampleCurve[a].Distance(resampleCurve[a-1]) << " ";
	}
	std::cout << "\n\n";*/
}

double UtilityFunctions::DistanceToFiniteLine(CVSystem::MyPoint v, CVSystem::MyPoint w, CVSystem::MyPoint p) 
{
	using namespace CVSystem;

	// Return minimum distance between line segment vw and point p
	double l2 = v.DistanceSquared(w);					   // i.e. |w-v|^2 -  avoid a sqrt
	if (l2 > -M_EPS && l2 < M_EPS) return p.Distance(v);   // v == w case

	// Consider the line extending the segment, parameterized as v + t (w - v).
	// We find projection of point p onto the line. 
	// It falls where t = [(p-v) . (w-v)] / |w-v|^2
	double t = (p - v).Dot(w - v) / l2;

	if (t < 0.0)	  { return  p.Distance(v); }       // Beyond the 'v' end of the segment
	else if (t > 1.0) { return  p.Distance(w); }  // Beyond the 'w' end of the segment
	MyPoint projection = v + (w - v) * t;     // Projection falls on the segment
	return p.Distance(projection);
}

/*
double UtilityFunctions::DistanceToInfinityLine(CVSystem::MyPoint v, CVSystem::MyPoint w, CVSystem::MyPoint p) 
{
	// http://stackoverflow.com/questions/15881718/perpendicular-distance-from-a-point-to-a-line-represented-by-two-points

	MyPoint AX = X - L.A;
	float DL = Dot(L.dir,AX);
	return sqrt(Dot(AX,AX) - DL*DL);

	//Coordinates are (a,b) and (c,d)
	//the point (x,y) is the required point.
	$a=1;
	$b=2;
	$c=3;
	$d=4;

	$m=($d-$b)/($c-$a);
	//echo $m."\n";

	$x=10;
	$y=20;
	//echo $y-($m*$x)-$b+($m*$a)."\n";
	$distance=abs($y-($m*$x)-$b+($m*$a))/sqrt(1+($m*$m));
	echo $distance;
}
*/

bool UtilityFunctions::IsClockwise(std::vector<CVSystem::MyPoint> polygon)
{
	using namespace CVSystem;
	double sumValue = 0;
	for(size_t a = 0; a < polygon.size(); a++)
	{
		MyPoint curPt = polygon[a];
		MyPoint nextPt;
		if(a == polygon.size() - 1) 
			{nextPt = polygon[0];}
		else
			{nextPt = polygon[a+1];}

		// sumValue += ((x2 - x1) * (y2 + y1));
		sumValue += ((nextPt.x - curPt.x) * (nextPt.y + curPt.y));
	}

	if(sumValue >= 0) { return false; }

	return true;
}

bool UtilityFunctions::IsClockwise(std::vector<int> polygon, std::vector<CVSystem::MyPoint> vertexList)
{

	using namespace CVSystem;
	double sumValue = 0;
	for(size_t a = 0; a < polygon.size(); a++)
	{
		MyPoint curPt = vertexList[polygon[a]];
		MyPoint nextPt;
		if(a == polygon.size() - 1) 
			{nextPt = vertexList[polygon[0]];}
		else
			{nextPt = vertexList[polygon[a+1]];}

		// sumValue += ((x2 - x1) * (y2 + y1));
		sumValue += ((nextPt.x - curPt.x) * (nextPt.y + curPt.y));
	}

	if(sumValue >= 0) { return false; }

	return true;
}

std::vector<CVSystem::MyLine> UtilityFunctions::PointsToLine(std::vector<CVSystem::MyPoint> points)
{
	using namespace CVSystem;
	std::vector<MyLine> lines;
	for(size_t a = 1; a < points.size(); a++)
		{ lines.push_back(MyLine(points[a - 1], points[a])); }
	lines.push_back(MyLine(points[points.size() - 1], points[0])); 
	return lines;
}

std::vector<CVSystem::MyIndexedLine> UtilityFunctions::IndexedPointsToLine(std::vector<int> points)
{
	using namespace CVSystem;
	std::vector<MyIndexedLine> lines;
	for(size_t a = 1; a < points.size(); a++)
	{ lines.push_back(MyIndexedLine(points[a - 1], points[a])); }
	lines.push_back(MyIndexedLine(points[points.size() - 1], points[0])); 
	return lines;
}

bool UtilityFunctions::DoesExist(std::vector<CVSystem::MyIndexedLine> data, CVSystem::MyIndexedLine line)
{
	using namespace CVSystem;
	for(size_t a = 0; a < data.size(); a++)
	{
		if(data[a].index0 == line.index0 && data[a].index1 == line.index1) {return true;}
		if(data[a].index1 == line.index0 && data[a].index0 == line.index1) {return true;}
	}

	return false;
}

std::vector<CVSystem::MyTriangle> 
	UtilityFunctions::ConvertIndexedTriangles(std::vector<CVSystem::MyIndexedTriangle> triangles, std::vector<CVSystem::MyPoint> vertexList)
{
	using namespace CVSystem;
	std::vector<MyTriangle>  trTemp;

	for(size_t a = 0; a < triangles.size(); a++)
	{
		MyPoint ptA = vertexList[triangles[a].idx0];
		MyPoint ptB = vertexList[triangles[a].idx1];
		MyPoint ptC = vertexList[triangles[a].idx2];

		MyTriangle tr(ptA, ptB, ptC,triangles[a].tri_type);
		trTemp.push_back(tr);
	}

	return trTemp;
}

void UtilityFunctions::ShiftIndices(std::vector<CVSystem::MyIndexedTriangle>& idxTriangles, int num)
{
	for(size_t a = 0; a < idxTriangles.size(); a++)
	{
		idxTriangles[a].idx0 += num;
		idxTriangles[a].idx1 += num;
		idxTriangles[a].idx2 += num;
	}
}

void UtilityFunctions::ShiftIndices(std::vector<CVSystem::MyIndexedBezierCurves>& bzCurves, int num)
{
	for(size_t a = 0; a < bzCurves.size(); a++)
	{
		for(size_t b = 0; b < bzCurves[a].indices.size(); b++)
			{ bzCurves[a].indices[b] += num; }
	}
}


std::vector<int> UtilityFunctions::GetDistinctIndices(std::vector<CVSystem::MyIndexedTriangle> triangles, int numVertices)
{
	std::vector<int> indices;

	std::vector<bool> flags;
	for(int a = 0; a < numVertices; a++)
		{ flags.push_back(false); }


	for(size_t a = 0; a < triangles.size(); a++)
	{
		int idx0 = triangles[a].idx0;
		int idx1 = triangles[a].idx1;
		int idx2 = triangles[a].idx2;

		if(!flags[idx0])
		{
			indices.push_back(idx0);
			flags[idx0] = true;
		}

		if(!flags[idx1])
		{
			indices.push_back(idx1);
			flags[idx1] = true;
		}

		if(!flags[idx2])
		{
			indices.push_back(idx2);
			flags[idx2] = true;
		}
	}


	return indices;
}

cv::Point2f UtilityFunctions::NextPoint(cv::Point2f p, double dis, double orientation){
	return cv::Point2f(p.x + dis * cos(UtilityFunctions::DegreeToRadian(orientation)), p.y + dis * sin(UtilityFunctions::DegreeToRadian(orientation)));
}

bool UtilityFunctions::InsideMask(cv::Mat mask, cv::Point2f p)
{
	if(mask.ptr<uchar>((int)p.y, (int)p.x)[0] > 127 )
		return true;
	return false;
}

bool UtilityFunctions::InsideImage(int imageWidth, int imageHeight, cv::Point2f p)
{
	if((p.x >= 0 && p.x <= imageWidth - 1) && (p.y >= 0 && p.y <= imageHeight - 1))
		return true;
	return false;
}

double UtilityFunctions::Dis2f(cv::Point2f p1, cv::Point2f p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

cv::Point2f UtilityFunctions::intersection(cv::Point2f p1s, cv::Point2f p1e, cv::Point2f p2s, cv::Point2f p2e)
{
	cv::Point2f v1 = p1e - p1s, v2 = p2e - p2s, s = p2s - p1s;

	// ��u����A���I���s�b�C
	// ��u���|�A���I�L���h�C
	if (v1.cross(v2) == 0) return cv::Point2f(-999, -999);

	// �p����I
	return p1s + v1 * (s.cross(v2) / v1.cross(v2));
}

void UtilityFunctions::ImageIntersectWithLinePoint(cv::Point &cp1, cv::Point &cp2, cv::Mat image, cv::Point p1, cv::Point p2)
{
	cv::Point2f LT(				0,				0);
	cv::Point2f LD(				0, image.rows - 1);
	cv::Point2f RT(image.cols - 1,				0);
	cv::Point2f RD(image.cols - 1, image.rows - 1);

	int Iedge = -1;
	/*
	 1
	  ___
	0|   |2
	 |   |
	  ___
	  3
	*/

	cv::Point2f tp2f = intersection(p1, p2, LT, LD);
	Iedge++;
	if(tp2f.y >= 0 && tp2f.y <= LD.y)
		cp1 = cv::Point(tp2f.x, tp2f.y);
	else
	{
		Iedge++;
		cv::Point2f tp2f = intersection(p1, p2, LT, RT);
		if(tp2f.x >= 0 && tp2f.x <= RT.x)
			cp1 = cv::Point(tp2f.x, tp2f.y);
		else
		{
			Iedge++;
			cv::Point2f tp2f = intersection(p1, p2, RT, RD);
			if(tp2f.y >= 0 && tp2f.y <= RD.y)
				cp1 = cv::Point(tp2f.x, tp2f.y);
			else
			{
				Iedge++;
				cv::Point2f tp2f = intersection(p1, p2, RD, LD);
				if(tp2f.x >= 0 && tp2f.x <= RD.x)
					cp1 = cv::Point(tp2f.x, tp2f.y);
				else
					cp1 = cv::Point(-1, -1);
			}
		}
	}

	tp2f = intersection(p1, p2, LT, LD);
	if(Iedge != 0 && tp2f.y >= 0 && tp2f.y <= LD.y)
		cp2 = cv::Point(tp2f.x, tp2f.y);
	else
	{
		cv::Point2f tp2f = intersection(p1, p2, LT, RT);
		if(Iedge != 1 && tp2f.x >= 0 && tp2f.x <= RT.x)
			cp2 = cv::Point(tp2f.x, tp2f.y);
		else
		{
			cv::Point2f tp2f = intersection(p1, p2, RT, RD);
			if(Iedge != 2 && tp2f.y >= 0 && tp2f.y <= RD.y)
				cp2 = cv::Point(tp2f.x, tp2f.y);
			else
			{
				cv::Point2f tp2f = UtilityFunctions::intersection(p1, p2, RD, LD);
				if(tp2f.x >= 0 && tp2f.x <= RD.x)
					cp2 = cv::Point(tp2f.x, tp2f.y);
				else
					cp2 = cv::Point(-1, -1);
			}
		}
	}
}

double UtilityFunctions::GetRotation(cv::Point2f pt1, cv::Point2f pt2)
{
	float perpDot = pt1.x * pt2.y - pt1.y * pt2.x;
	return (float)atan2(perpDot, pt1.dot(pt2));
}

double UtilityFunctions::GetAngleFromXAxis(cv::Point2f pt)
{
	return GetRotation(cv::Point2f(1, 0), pt);
}

void UtilityFunctions::splitBezierCubic(cv::Point2d orgCP[4], cv::Point2d newCPA[4], cv::Point2d newCPB[4], double t)
{
	newCPA[0] = orgCP[0];
	cv::Point2d orgCPc[4];
	for(int i = 0; i < 4; i++) 
		orgCPc[i] = orgCP[i];

	cv::Point2d cp01 = (orgCPc[1]-orgCPc[0])*t + orgCPc[0];
	cv::Point2d cp12 = (orgCPc[2]-orgCPc[1])*t + orgCPc[1];
	cv::Point2d cp23 = (orgCPc[3]-orgCPc[2])*t + orgCPc[2];
	cv::Point2d cp012 = (cp12 - cp01)*t + cp01;
	cv::Point2d cp123 = (cp23-cp12)*t + cp12;
	cv::Point2d cp0123 = (cp123 - cp012)*t + cp012;

	newCPA[1] = cp01;
	newCPA[2] = cp012;
	newCPA[3] = cp0123;

	newCPB[0] = cp0123;
	newCPB[1] = cp123;
	newCPB[2] = cp23;
	newCPB[3] = orgCPc[3];
}

void UtilityFunctions::splitBezierCubic( cv::Point2d orgCP[4], cv::Point2d resCP[4], double t1, double t2 )
{
	cv::Point2d tmp_cp1[4], tmp_cp2[4];

	splitBezierCubic(orgCP, tmp_cp1, tmp_cp2, t2);
	splitBezierCubic(tmp_cp1, tmp_cp2, resCP, t1/t2);
}

double UtilityFunctions::determinant_3x3( double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21, double m22 )
{
	return m00 * m11 * m22 +
		m01 * m12 * m20 +
		m02 * m10 * m21 -
		m02 * m11 * m20 -
		m01 * m10 * m22 -
		m00 * m12 * m21;
}

int UtilityFunctions::PointInLine2D( const cv::Point2d& p, const cv::Point2d& Lp1, const cv::Point2d& Lp2 )
{
	const double num_err = 1e-10;
	double Lx = Lp2.x - Lp1.x;
	double Ly = Lp2.y - Lp1.y;

	if ( abs(Lx) >= num_err && abs(Ly) >= num_err )
	{
		double Lpx = (p.x - Lp1.x) / Lx;
		double Lpy = (p.y - Lp1.y) / Ly;
		if ( abs( Lpx-Lpy ) < 1e-6 ) //colinear
		{
			if ( Lpx>1 || Lpx<0 ) // outside
				return -1;
			else
				return 1;
		}
		else
			return 0;
	}
	else if ( abs(Lx) <= num_err ) // perpendicular
	{
		if ( abs(p.x - Lp1.x) >= num_err )
		{
			return 0;
		}
		else
		{
			double Lpy = (p.y-Lp1.y) / Ly;
			if ( Lpy>1 || Lpy<0 )
				return -1;
			else
				return 1;
		}
	}
	else if ( abs(Ly) <= num_err ) // horizontal
	{
		double dif_y = abs(p.y - Lp1.y);
		if ( dif_y >= num_err )
		{
			return 0;
		}
		else
		{
			double Lpx = (p.x-Lp1.x) / Lx;
			if ( Lpx>1 || Lpx<0 )
				return -1;
			else
				return 1;
		}
	}
	return 0;
}

double UtilityFunctions::getLinePara( const cv::Point2d &p, const cv::Point2d &Lp1, const cv::Point2d &Lp2 )
{
	if ( Lp1.x != Lp2.x )
		return ( p.x - Lp1.x ) / ( Lp2.x - Lp1.x );
	else if ( Lp1.y != Lp2.y )
		return ( p.y - Lp1.y ) / ( Lp2.y - Lp1.y );

	return 0;
}

void UtilityFunctions::GetBaryCoord( cv::Point3d &bc, cv::Point2d &p, cv::Point2d &p1, cv::Point2d &p2, cv::Point2d &p3 )
{
	// 	double va = (p1 - p2).length();
	// 	double vb = (p2 - p3).length();
	// 	double vc = (p3 - p1).length();
	// 	double totalarea = triarea(va, vb, vc);
	// 	double len0 = (p1 - p).length();
	// 	double len1 = (p2 - p).length();
	// 	double len2 = (p3 - p).length();
	// 	double tarea[3];
	// 	tarea[0] = triarea(vb, len1, len2);
	// 	tarea[1] = triarea(vc, len0, len2);
	// 	bc[0] = tarea[0]/totalarea;
	// 	bc[1] = tarea[1]/totalarea;
	// 	bc[2] = 1 - bc[1] - bc[0];
	double cp = (p2.y-p3.y)*(p1.x-p3.x) + (p3.x-p2.x)*(p1.y-p3.y);
	double den;
	if (abs(cp)<1e-10)
		cp = 1e-10*signd(cp);
	den = 1/cp;

	bc.x = ((p2.y-p3.y)*(p.x-p3.x) + (p3.x-p2.x)*(p.y-p3.y)) * den;
	bc.y = ((p3.y-p1.y)*(p.x-p3.x) + (p1.x-p3.x)*(p.y-p3.y)) * den;
	bc.z = 1 - bc.x - bc.y;
}

int UtilityFunctions::LineIntersect( cv::Point2d &cross_p, cv::Point2d &L1_p1, cv::Point2d &L1_p2, cv::Point2d &L2_p1, cv::Point2d &L2_p2 )
{
	double parallelCheck = (L1_p1.x-L1_p2.x)*(L2_p1.y-L2_p2.y) - (L1_p1.y-L1_p2.y)*(L2_p1.x-L2_p2.x);
	if ( fabs(parallelCheck)<1e-20 )
		return PARALLEL;
	cross_p.x = ( (L1_p1.x*L1_p2.y-L1_p1.y*L1_p2.x)*(L2_p1.x-L2_p2.x) - (L2_p1.x*L2_p2.y-L2_p1.y*L2_p2.x)*(L1_p1.x-L1_p2.x) )/parallelCheck;
	cross_p.y = ( (L1_p1.x*L1_p2.y-L1_p1.y*L1_p2.x)*(L2_p1.y-L2_p2.y) - (L2_p1.x*L2_p2.y-L2_p1.y*L2_p2.x)*(L1_p1.y-L1_p2.y) )/parallelCheck;

	if ( cross_p.x>L1_p1.x && cross_p.x>L1_p2.x )
		return CROSS_OUTSIDE;
	if ( cross_p.x<L1_p1.x && cross_p.x<L1_p2.x )
		return CROSS_OUTSIDE;
	if ( cross_p.y>L1_p1.y && cross_p.y>L1_p2.y )
		return CROSS_OUTSIDE;
	if ( cross_p.y<L1_p1.y && cross_p.y<L1_p2.y )
		return CROSS_OUTSIDE;

	if ( cross_p.x>L2_p1.x+1e-10 && cross_p.x>L2_p2.x+1e-10 )
		return CROSS_OUTSIDE;
	if ( cross_p.x<L2_p1.x-1e-10 && cross_p.x<L2_p2.x-1e-10 )
		return CROSS_OUTSIDE;
	if ( cross_p.y>L2_p1.y && cross_p.y>L2_p2.y )
		return CROSS_OUTSIDE;
	if ( cross_p.y<L2_p1.y && cross_p.y<L2_p2.y )
		return CROSS_OUTSIDE;

	return CROSS_INSIDE;
}

bool UtilityFunctions::PointTriangleInside(cv::Point3d &w, cv::Point2d &pos, cv::Point2d &tp0, cv::Point2d &tp1, cv::Point2d &tp2)
{
	GetBaryCoord( w, pos, tp0, tp1, tp2 );
	if ( w.x<-PT_TRI_INSIDE_ERROR || w.x-1>PT_TRI_INSIDE_ERROR )
		return false;
	if ( w.y<-PT_TRI_INSIDE_ERROR || w.y-1>PT_TRI_INSIDE_ERROR )
		return false;
	if ( w.z<-PT_TRI_INSIDE_ERROR || w.z-1>PT_TRI_INSIDE_ERROR )
		return false;
	return true;
}

bool UtilityFunctions::PointTriangleInside(cv::Point3d &w, cv::Point2d pos, cv::Point2d tp[3])
{
	return PointTriangleInside(w, pos, tp[0], tp[1], tp[2]);
}

bool UtilityFunctions::CubicToLineApx( cv::Point2d cp[4] )
{
	double dist = Dis2f(cp[3], cp[0]);
	double d1 = PointLineDistance( cp[1], cp[0], cp[3] );
	double d2 = PointLineDistance( cp[2], cp[0], cp[3] );
	const double error = 1e-2;
	if ( dist < 2 )	//point
		return true;
	else if ( d1/dist < error && d2/dist < error )
		return true;

	return false;
}

double UtilityFunctions::PointLineDistance( cv::Point2d& p, cv::Point2d& L_p1, cv::Point2d& L_p2 )
{
	double dx = L_p1.x - L_p2.x;
	double dy = L_p1.y - L_p2.y;

	return fabs( dy*p.x - dx*p.y + L_p1.x*L_p2.y - L_p1.y*L_p2.x ) / sqrt( dx*dx + dy*dy );
}

cv::Point2d UtilityFunctions::GetCubicPt(cv::Point2d p1, cv::Point2d p2, cv::Point2d p3, cv::Point2d p4, double t )
{
	return (1-t)*(1-t)*(1-t)*p1 + 3*(1-t)*(1-t)*t*p2 + 3*(1-t)*t*t*p3 + t*t*t*p4;
}

cv::Point2d UtilityFunctions::GetCubicDerivative( cv::Point2d p1, cv::Point2d p2, cv::Point2d p3, cv::Point2d p4, double t )
{
	cv::Point2d g = (p2-p1)*3;
	cv::Point2d b = (p3-p2)*3 - g;
	cv::Point2d a = p4 - p1 - b - g;
	return 3*t*t*a + 2*t*b + g;
}

void UtilityFunctions::CubicRoots( double t[3], double p[4] )
{
	double A = p[1]/p[0];
	double B = p[2]/p[0];
	double C = p[3]/p[0];

	double Q = (3.0*B - A*A) / 9.0;
	double R = (9.0*A*B - 27.0*C - 2.0*A*A*A) / 54.0;
	double D = Q*Q*Q + R*R;	// polynomial discriminant

	if ( D > 0)	//complex or duplicate roots
	{
		double D_sqr = sqrt(D);
		double S = signd( R + D_sqr )*pow(abs(R+D_sqr), 1.0/3.0);
		double T = signd( R - D_sqr )*pow(abs(R-D_sqr), 1.0/3.0);

		t[0] = -A/3.0 + (S + T);	// real root
		t[1] = -A/3.0 - (S+T)/2.0;	// real part of complex root
		t[2] = t[1];			// real part of complex root
		double Im = fabs( sqrt(3.0)*(S-T)/2.0 );	/// complex part of root pair  

		//discard complex roots
		if (Im!=0)
		{
			t[1] = -1;
			t[2] = -1;
		}
	}
	else
	{
		double th = acos( R / sqrt( -Q*Q*Q ) );
		t[0] = 2*sqrt(-Q)*cos(th/3) - A/3;
		t[1] = 2*sqrt(-Q)*cos( (th + 2*M_PI)/3 ) - A/3;
		t[2] = 2*sqrt(-Q)*cos( (th + 4*M_PI)/3 ) - A/3;
		//printf("t2 = %lf\n", t[2]);
	}
}

bool UtilityFunctions::PointLineInside( cv::Point2d point, cv::Point2d Line_p1, cv::Point2d Line_p2 )
{
	double para = GetLinePara(point, Line_p1, Line_p2);
	if ( para >=0 && para <= 1 ) return true;

	return false;
}

double UtilityFunctions::GetLinePara( cv::Point2d p, cv::Point2d Lp1, cv::Point2d Lp2 )
{
	if ( Lp1.x != Lp2.x )
		return ( p.x - Lp1.x ) / ( Lp2.x - Lp1.x );
	else if ( Lp1.y != Lp2.y )
		return ( p.y - Lp1.y ) / ( Lp2.y - Lp1.y );

	return 0;
}

bool UtilityFunctions::AabbIntersect( const cv::Point2d p11, const cv::Point2d p12, const cv::Point2d p21, const cv::Point2d p22 )
{
	cv::Point2d max1 = p11;
	cv::Point2d min1 = p11;

	cv::Point2d max2 = p21;
	cv::Point2d min2 = p21;

	if ( p12.x > max1.x ) max1.x = p12.x;
	if ( p12.y > max1.y ) max1.y = p12.y;
	if ( p12.x < min1.x ) min1.x = p12.x;
	if ( p12.y < min1.y ) min1.y = p12.y;

	if ( p22.x > max2.x ) max2.x = p22.x;
	if ( p22.y > max2.y ) max2.y = p22.y;
	if ( p22.x < min2.x ) min2.x = p22.x;
	if ( p22.y < min2.y ) min2.y = p22.y;

	if ( max1.x <  min2.x)
		return false;
	if ( max1.y < min2.y )
		return false;

	if ( min1.x >  max2.x)
		return false;
	if ( min1.y > max2.y )
		return false;

	return true;
}