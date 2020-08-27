#include "../../Include/Common.h"

std::string MakeFileNameWithFlag(std::string fileName, int dpi, std::string flag)
{
	std::stringstream temp_dpi;
	temp_dpi << dpi;

	std::string resize_dpi;
	temp_dpi >> resize_dpi;

	resize_dpi = "_" + resize_dpi;

	int find_number = fileName.find_first_of(".");
	std::string mainName = fileName.substr(0, find_number);
	std::string fileFlag = resize_dpi + flag;

	return mainName + fileFlag + ".png";
}

int CVSystem::ScreentoneSegmentation::Threshold(int* hist) //compute the threshold[OTSU algorithm]
{
	float u0, u1;
	float w0, w1;
	int count0;
	int t, maxT;
	float devi, maxDevi = 0; //��t�γ̤j��t
	int i;
	int sum = 0;
	for (i = 0; i < 256; i++)
	{
		sum = sum + hist[i];
	}
	for (t = 0; t < 255; t++)
	{
		u0 = 0;
		count0 = 0; //�H�Ȭ�t�ɡAc0�ժ����Ȥβ��ͪ����v
		for (i = 0; i <= t; i++)
		{
			u0 += i * hist[i];
			count0 += hist[i];
		}
		u0 = u0 / count0;
		w0 = (float)count0 / sum; //�H�Ȭ�t�ɡAc1�ժ����Ȥβ��ͪ����v
		u1 = 0;
		for (i = t + 1; i < 256; i++)
		{
			u1 += i * hist[i];
		}
		u1 = u1 / (sum - count0);
		w1 = 1 - w0; //��������t
		devi = w0 * w1 * (u1 - u0) * (u1 - u0); //�O���̤j����t�γ̨Φ�m
		if (devi > maxDevi)
		{
			maxDevi = devi;
			maxT = t;
		}
	}
	return maxT;
}

void CVSystem::ScreentoneSegmentation::ComputeOtsuGaussian()
{
	std::cout << "�i�hOtsu_Gaussian" << std::endl;

	//��Ƕ�	
	//cv::Mat gray_image = this->oriInpImg.clone();
	cv::Mat gray_image;
	cv::cvtColor(this->oriInpImg, gray_image, CV_BGR2GRAY);

	std::cout << "�i�hcvtColor" << std::endl;

	//�έp�����
	int hist[256] = { 0 };
	for (int i = 0; i < gray_image.cols; i++)
	{
		for (int j = 0; j < gray_image.rows; j++)
		{
			int xx = (int)gray_image.at<uchar>(j + i * gray_image.rows);
			hist[xx] += 1;
		}
	}

	//�p��֭�
	int t = Threshold(hist);
	if (t == 0)
	{
		std::cout << "done" << std::endl;
		//�g�X�G�ȤƹϤ�
		std::string ori_fileName = this->filename.substr(this->filename.find_last_of("/") + 1);
		std::string FileName = MakeFileNameWithFlag(ori_fileName, 1200, "_B");
		binary_output = SystemParams::str_Resources_Binarization + "/" + FileName;
		cv::imwrite(binary_output, gray_image);
		return;
	}
	std::cout << "�H�Ȭ��G" << t << std::endl;

	//*----------------------------�ƭȲέp---------------------------------*
	int not_zero = 0;
	for (int i = 0; i < t + 1; i++)
	{
		if (hist[i] != 0)
		{
			not_zero = i;
			break;
		}
	}
	std::cout << "�}�Y���G" << not_zero << std::endl;

	int all_size = 0;
	double all_temp = 0;
	for (int i = 0; i < 256; i++)
	{
		all_temp += i * hist[i];
		all_size += hist[i];
	}
	double all_mean = all_temp / all_size;
	double all_temp1 = 0;
	for (int i = 0; i < 256; i++)
	{
		int x = (i - all_mean) * (i - all_mean) * hist[i];
		all_temp1 += x;
	}
	double all_temp2 = all_temp1 / all_size;
	double all_std = sqrt((double)all_temp2);

	std::cout << "all_mean���G" << all_mean << " all_std���G" << all_std << std::endl;


	int white_max = 0;
	int temp0 = hist[0];
	int white_sum = 0;
	int sum0 = 0;
	for (int i = 0; i < t; i++)
	{
		if (hist[i] > temp0)
		{
			temp0 = hist[i];
			white_max = i;
		}
		sum0 += hist[i];
		white_sum += i * hist[i];
	}
	int white_mean = white_sum / sum0;
	double w_sum0 = 0;
	for (int i = 0; i < t + 1; i++)
	{
		int x = (i - white_mean) * (i - white_mean) * hist[i];
		w_sum0 += x;
	}
	double w_temp = w_sum0 / sum0;

	double white_std = sqrt((double)w_temp);
	std::cout << "white_max���G" << white_max << " �ƶq���G" << temp0 << std::endl;
	std::cout << "white_sum���G" << white_sum << " �ƶq���G" << sum0 << " mean:" << white_mean << " std:" << white_std << std::endl;

	//��White max �H�e��mean �� std
	int w_m_size = 0;
	int w_m_sum = 0;
	for (int i = 0; i < white_max + 1; i++)
	{
		w_m_size += hist[i];
		w_m_sum += i * hist[i];
	}
	int w_m_mean = w_m_sum / w_m_size;
	std::cout << "white_max���e�������Ȭ��G" << w_m_mean << std::endl;
	int w_m_sum0 = 0;
	for (int i = 0; i < white_max + 1; i++)
	{
		int x = (i - w_m_mean) * (i - w_m_mean) * hist[i];
		w_m_sum0 += x;
	}
	int w_m_temp0 = w_m_sum0 / w_m_size;
	double w_m_std = sqrt((double)w_m_temp0);
	std::cout << "white_max���e���зǮt���G" << w_m_std << std::endl;

	int black_max = 0;
	int temp1 = hist[t];
	double black_sum = 0;
	int sum1 = 0;
	for (int i = t; i < 256; i++)
	{
		if (hist[i] > temp1)
		{
			temp1 = hist[i];
			black_max = i;
		}
		sum1 += hist[i];
		black_sum += i * hist[i];
	}
	double black_mean = black_sum / sum1;
	double b_sum0 = 0;
	for (int i = t; i < 256; i++)
	{
		int x = (i - black_mean) * (i - black_mean) * hist[i];
		b_sum0 += x;
	}
	double b_temp = b_sum0 / sum1;
	double black_std = sqrt((double)b_temp);
	std::cout << "black_max���G" << black_max << " �ƶq���G" << temp1 << std::endl;
	std::cout << "black_sum���G" << black_sum << " �ƶq���G" << sum1 << " mean:" << black_mean << " std:" << black_std << std::endl;

	//��black max �H�e��mean �� std
	int b_m_size = 0;
	int b_m_sum = 0;
	for (int i = black_max; i < 256; i++)
	{
		b_m_size += hist[i];
		b_m_sum += i * hist[i];
	}
	int b_m_mean = b_m_sum / b_m_size;
	std::cout << "black_max���᪺�����Ȭ��G" << b_m_mean << std::endl;
	int b_m_sum0 = 0;
	for (int i = black_max; i < 256; i++)
	{
		int x = (i - b_m_mean) * (i - b_m_mean) * hist[i];
		b_m_sum0 += x;
	}
	int b_m_temp0 = b_m_sum0 / b_m_size;
	double b_m_std = sqrt((double)b_m_temp0);
	std::cout << "black_max���᪺�зǮt���G" << b_m_std << std::endl;


	//�S��

	int black_max_sigma = (int)white_max + (int)white_std + SystemParams::OG_sigma;
	int white_max_sigma = (int)black_max - (int)white_std + SystemParams::OG_sigma;
	std::cout << black_max_sigma << " : " << white_max_sigma << std::endl;
	if (SystemParams::OG_Bsigma != 0)
		black_max_sigma = SystemParams::OG_Bsigma;
	if (SystemParams::OG_Wsigma != 0)
		white_max_sigma = SystemParams::OG_Wsigma;

	for (int i = 0; i < gray_image.cols; i++)
	{
		for (int j = 0; j < gray_image.rows; j++)
		{
			int xx = (int)gray_image.at<uchar>(j, i);

			if (xx < black_max_sigma)
			{
				gray_image.at<uchar>(j, i) = 0;
			}
			if (xx > white_max_sigma)
			{
				gray_image.at<uchar>(j, i) = 255;
			}
		}
	}

	//�����G�Ȥ�
	//�����local����
	cv::adaptiveThreshold(gray_image, local, 255, CV_ADAPTIVE_THRESH_GAUSSIAN_C, CV_THRESH_BINARY, SystemParams::G_block_size, SystemParams::G_params);

	std::cout << "done" << std::endl;
	//�g�X�G�ȤƹϤ�
	std::string ori_fileName = this->filename.substr(this->filename.find_last_of("/") + 1);
	std::string fileName = MakeFileNameWithFlag(ori_fileName, 1200, "_B");
	binary_output = SystemParams::str_Resources_Binarization + "/" + fileName;
	cv::imwrite(binary_output, local);
}

cv::Mat CVSystem::ScreentoneSegmentation::Compute_Segmentation(cv::Mat image)
{
	cv::Mat outSegm;
	image.convertTo(outSegm, CV_8UC1);
	if (SystemParams::TU_maskSize >= 1)
	{
		cv::Mat elemE = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(SystemParams::TU_maskSize, SystemParams::TU_maskSize));
		cv::morphologyEx(outSegm, outSegm, cv::MORPH_CLOSE, elemE);
		cv::morphologyEx(outSegm, outSegm, cv::MORPH_ERODE, elemE);
		outSegm.convertTo(outSegm, CV_8UC1);
	}
	return outSegm.clone();
}