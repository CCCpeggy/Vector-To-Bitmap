
/*
 * Sorting lines to be oriented counter-clockwise and clockwise
 *
 * Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 */

#ifndef __Lines_Sorter__
#define __Lines_Sorter__

#include "../../Include/Common.h"

namespace CVSystem
{
	// rebuild orientation of boundary curves
	class LinesSorter
	{
	public:
		// old version. can cause infinite loop
		static std::vector<MyIndexedBezierCurves> BuildBezierCurves1(std::vector<MyIndexedLine> inputLines, 
															        std::vector<MyPoint> vertexList);

		// new version
		static std::vector<MyIndexedBezierCurves> BuildBezierCurves2(std::vector<MyIndexedLine> inputLines, 
																	std::vector<MyPoint> vertexList);
	};
}
#endif