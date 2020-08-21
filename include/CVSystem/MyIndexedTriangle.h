
/**
 *
 * Triangles with three vertices
 * Each vertices doesn't hold the real coordinate
 *
 * Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
 * Version: 2014
 *
 *
 */

#ifndef __My_Indexed_Triangles__
#define __My_Indexed_Triangles__

#include <opencv2/core/core.hpp>
#include "MyPoint.h"
#include "TriangleType.h"

namespace CVSystem
{
	
	struct MyIndexedTriangle
	{	
	public:
		// first index
		int idx0;

		// second index
		int idx1;

		// third index
		int idx2;

		// type of the triangle (black, white, and screentone)
		TriangleType tri_type;

		// Constructor
		MyIndexedTriangle(int idx0, int idx1, int idx2, TriangleType tri_type)
		{
			this->idx0 = idx0;
			this->idx1 = idx1;
			this->idx2 = idx2;
			this->tri_type = tri_type;
		}
	};

	class KLMTriangle : public MyIndexedTriangle
	{
	public:

		cv::Vec3d KLM[3];

		unsigned short int KLMTriType;//0B 1W 2KLM

		cv::Vec4d texCoord[3];

		KLMTriangle()
			:MyIndexedTriangle(0, 0, 0, TRIANGLE_NOT_DEFINED)
		{
			for(int i = 0; i < 3; i++)
				this->KLM[i] = 0;

			this->KLMTriType = 0;

			for(int i = 0; i < 3; i++)
				this->texCoord[i] = cv::Vec4d(0, 0, 0, 0);
		}

		KLMTriangle(int idx0, int idx1, int idx2, int KTT, cv::Vec3d KLM[3])
			:MyIndexedTriangle(idx0, idx1, idx2, TRIANGLE_NOT_DEFINED)
		{
			for(int i = 0; i < 3; i++)
				this->KLM[i] = KLM[i];

			this->KLMTriType = KTT;

			if(KTT != 2)
			{
				for(int i = 0; i < 3; i++)
					this->texCoord[i] = cv::Vec4d(KTT, KTT, KTT, KTT);
			}
		}

		KLMTriangle(int idx0, int idx1, int idx2, int KTT, cv::Vec3d KLM[3], cv::Vec4d texCoord[3])
			:MyIndexedTriangle(idx0, idx1, idx2, TRIANGLE_NOT_DEFINED)
		{
			for(int i = 0; i < 3; i++)
				this->KLM[i] = KLM[i];

			this->KLMTriType = KTT;

			for(int i = 0; i < 3; i++)
				this->texCoord[i] = texCoord[i];
		}
	};
}
#endif