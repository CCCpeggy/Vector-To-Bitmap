
/**
* Screentone segmentation class with two different approach:
*	- Level set and gabor wavelet (Chan-Vese algorithm)
*	- Cartoon+texture filter
*
* Author: Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
* Version: 2014
*
*
*/

#ifndef __Screentone_Segmentation__
#define __Screentone_Segmentation__

#include "../../Include/Common.h"

//Tesseract OCR
//#include <baseapi.h>
//#include <resultiterator.h>
//
//#include "MyLine.h"
//
//#include "GaborManjunathMa.h"
//#include "RunningStat.h"


namespace CVSystem
{
	class ScreentoneSegmentation
	{
		
	public:
		
		// Constructor
		ScreentoneSegmentation();

		// Destructor
		~ScreentoneSegmentation();

		// Cartoon+texture segmentation
		//void ComputeCTSegmentation();

		//void ComputeCTSegmentation(std::string fileName, bool ShowFrame);

		//// level set segmentation
		//void ComputeLSSegmentation(std::vector<std::vector<CVSystem::MyLine>> sLines);	

		//// DIstance map (;evel set)
		//void ComputeDistanceMap(std::vector<std::vector<CVSystem::MyLine>> sLines);

		//// init local variables
		//void Init(std::string filename, int width, int height);
		//void Init(cv::Mat image);

		//// Get label map
		//int* GetLabelMap() { return labelMap; }
		//int* GetLabelMap(std::string fileName);
		//void FreeLabelMap();

		//// get dilated label map
		//int* GetDilatedLabelMap() { return dilatedLabelMap; }
		//int* GetDilatedLabelMap(std::string fileName);
		//void FreeDilatedLabelMap();

		////Otsu
		int Threshold(int *hist);
		void ComputeOtsuGaussian();  
		//cv::Mat ComputeOtsuGaussian(std::string fileName, bool isExport);
		//void SetOutputBinaryStr(std::string filePath) { binary_output = filePath; }	//回傳檔案位置
		//std::string GetOutputBinaryStr() { return this->binary_output; }	//回傳檔案位置

		////OCR
		//void ComputeTesseractOCR(std::string fileName);

		//void ExportImage(std::string fileName, cv::Mat image);

		static cv::Mat Compute_Segmentation(cv::Mat image);

		//void MargeCFR(cv::Mat oriBImage, std::vector<cv::Mat> cfrList);
		
	private:
		//// refine segmentation (cartoon+texture filter)
		//void RemoveSmallArea1(cv::Mat& segm);

		//// refine segmentation (cartoon+texture filter)
		//void RemoveSmallArea2(cv::Mat& segm);

		//// gabor wavelet feature (level set)
		//std::vector<double> GetStrokeFeature(std::vector<MyPoint> sPts);

		//// gabor wavelet feature (level set)
		//std::vector<double> GetPixelFeature(MyPoint p);

		//// get distance of two features (level set)
		//double GetDistance(std::vector<double> f1, std::vector<double> f2);

		//// get distance of two features (level set)
		//double GetDistance(std::vector<double> fv, int x, int y);

		//// distance map (level set)
		//void   GetDistanceMap(std::vector<double> fv);		

		//void   BuildFeature();

		//// standard deviation to construct gabor wavelet feature (level set)
		//void FeatureStdDev(cv::Mat fv, std::vector<double>& stds);	// one pass implementation		

	private:
		int* labelMap;
		int* dilatedLabelMap;

		cv::Mat featureMatrix;

		std::vector<CVSystem::MyPoint> sPts;
		std::vector<double> stds;
		cv::Mat gbrImg;

		
		std::string filename;
		cv::Mat inpImg;
		cv::Mat oriInpImg;
		cv::Mat binImg;

		cv::Mat img[5];

		cv::Mat local;
		cv::Mat ori_image;
		std::string binary_output;
		int resize_time;

		int	img_ori_width;
		int img_ori_height;
		int	img_width_scale;
		int img_height_scale;

		bool shouldCalc;
	};
}

#endif
