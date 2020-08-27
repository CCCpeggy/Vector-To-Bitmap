
/**
*
* Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
* Version: 2014
*
*/

#include "../../Include/Common.h"


double CVSystem::CurveRDP::PerpendicularDistance(CVSystem::MyPoint p, CVSystem::MyPoint p1, CVSystem::MyPoint p2) 
{
	// if start and end point are on the same x the distance is the difference in X.
	double result;
	if (abs(p1.x - p2.x) < M_EPS) { result = abs(p.x - p1.x); }
	else
	{
		double slope = (p2.y - p1.y) / (p2.x - p1.x); // a
		double intercept = p1.y - (slope * p1.x); // b
		result = abs(slope * p.x - p.y + intercept) / sqrt(pow(slope, 2) + 1);
	}

	return result;
}

void CVSystem::CurveRDP::RDP(std::vector<bool>& flags, std::vector<CVSystem::MyPoint> points, double epsilon, int startIndex, int endIndex, double rdp_point_min)
{			
	//if (endIndex - startIndex < rdp_point_min)
	//	{ return; }

	MyPoint firstPoint = points[startIndex];
	MyPoint lastPoint=points[endIndex];

	int index = -1;
	double dist = DBL_MIN;
	for (int i = startIndex + 1; i < endIndex; i++)
	{
		double cDist = PerpendicularDistance(points[i], firstPoint, lastPoint);
		if (cDist > dist)
		{
			dist = cDist;
			index = i;
		}
	}

	if (dist > epsilon)
	{
		if(index - startIndex >= rdp_point_min && endIndex - index >= rdp_point_min)
		{
			// Todo: wrong recursive sequence
			flags[index] = true;
			RDP(flags, points, epsilon, startIndex, index, rdp_point_min);
			RDP(flags, points, epsilon, index, endIndex  , rdp_point_min);
		}

		/*flags[index] = true;
		RDP(flags, points, epsilon, startIndex, index, rdp_point_min);
		RDP(flags, points, epsilon, index, endIndex  , rdp_point_min);*/
	}
}

void CVSystem::CurveRDP::RDPst(std::vector<std::vector<MyPoint>> polys, std::vector<bool>& flags, std::vector<CVSystem::MyPoint> &points, double epsilon, int startIndex, int endIndex, double rdp_point_min)
{	
	MyPoint mindisP;
	double dis = 10;
	int minIndex = -1;

	int index = -1;
	for (int i = startIndex + 1; i < endIndex; i++)
	{

		MyPoint tp = points[i];
		for(int ps = 0; ps < polys.size(); ps++)
		{
			for(int j = 0; j < polys[ps].size(); j++)
			{
				double tempDis = tp.Distance(polys[ps][j]);
				if(tempDis < dis)
				{
					mindisP = polys[ps][j];
					dis = tempDis;
					minIndex = i;
				}
			}
		}
	}
	if(dis < SystemParams::t_STDis)
	{
		points[minIndex] = mindisP;
		index = minIndex;
		flags[index] = true;
		RDPst(polys, flags, points, epsilon, startIndex, index, rdp_point_min);
		RDPst(polys, flags, points, epsilon, index, endIndex  , rdp_point_min);
	}

}

void CVSystem::CurveRDP::RDPst2(std::vector<bool>& flags, std::vector<CVSystem::MyPoint> &points, double epsilon, int startIndex, int endIndex, double rdp_point_min)
{	
	std::vector<int> stflags;
	for(int i = 0; i < flags.size() - 1; i++)
	{
		if(flags[i])
			stflags.push_back(i);
	}
	for(int i = 0; i < stflags.size() - 1; i++)
	{
		RDP(flags, points, epsilon, stflags[i], stflags[i + 1], rdp_point_min);
	}
}

void CVSystem::CurveRDP::SimplifyRDP(std::vector<CVSystem::MyPoint>& oldCurves, std::vector<CVSystem::MyPoint>& newCurves, double epsilon)
{
	newCurves.clear();
	newCurves.push_back(oldCurves[0]);
	SimplifyRDPRecursive(oldCurves, newCurves, epsilon, 0, oldCurves.size() - 1);
	newCurves.push_back(oldCurves[oldCurves.size() - 1]);
}
void CVSystem::CurveRDP::SimplifyRDPRecursive(std::vector<CVSystem::MyPoint>& oldCurves, std::vector<CVSystem::MyPoint>& newCurves, double epsilon, int startIndex, int endIndex)
{
	MyPoint firstPoint = oldCurves[startIndex];
	MyPoint lastPoint = oldCurves[endIndex];

	int index = -1;
	double dist = DBL_MIN;
	for (int i = startIndex + 1; i < endIndex; i++)
	{
		double cDist = PerpendicularDistance(oldCurves[i], firstPoint, lastPoint);
		if (cDist > dist)
		{
			dist = cDist;
			index = i;
		}
	}

	if (index != -1 && dist > epsilon)
	{
		SimplifyRDPRecursive(oldCurves, newCurves, epsilon, startIndex, index);
		newCurves.push_back(oldCurves[index]);		
		SimplifyRDPRecursive(oldCurves, newCurves, epsilon, index, endIndex);
	}
}